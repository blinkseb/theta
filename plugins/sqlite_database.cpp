#include "plugins/sqlite_database.hpp"
#include "interface/plugin.hpp"
#include "interface/histogram.hpp"

#include <sstream>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/scoped_array.hpp>

using namespace std;
using namespace theta;

sqlite_database::sqlite_database(const Configuration & cfg) :
    db(0), transaction_active(false), save_all_products(true){
    std::string filename = cfg.setting["filename"];
    if(cfg.setting.exists("products_data")){
        if(cfg.setting["products_data"].get_type()==Setting::TypeString){
            string s = cfg.setting["products_data"];
            if(s=="*")save_all_products = true;
            else throw ConfigurationException("products_data setting is a string but not '*'");
        }
        else{
            save_all_products = false;
            size_t n = cfg.setting["products_data"].size();
            for(size_t i=0; i<n; ++i){
                string column_name = cfg.setting["products_data"][i];
                products_data.insert(column_name);
                if(column_name=="*"){
                    save_all_products = true;
                    products_data.clear();
                    break;
                }
            }
            //if anything is written at all, also write runid and eventid:
            if(products_data.size()){
                products_data.insert("runid");
                products_data.insert("eventid");
            }
        }
    }
    if (boost::filesystem::exists(filename)) {
        boost::filesystem::remove(filename);
    }
    int res = sqlite3_open(filename.c_str(), &db);
    if (res != SQLITE_OK) {
        stringstream ss;
        ss << "sqlite_database constructor failed (filename=\"" << filename << "\") with SQLITE code " << res;
        error(ss.str());//throws.
    }
    beginTransaction();
}

void sqlite_database::close() {
    if (!db)
        return;
    endTransaction();
    //finalize all statements associated with the connection:
    sqlite3_stmt *pStmt;
    while ((pStmt = sqlite3_next_stmt(db, 0)) != 0) {
        sqlite3_finalize(pStmt);
    }
    int res = sqlite3_close(db);
    //set db to zero even in case of failure as it is unlikely
    // that the error is recoverable ...
    db = 0;
    if (res != 0) {
        error("sqlite_database::close()");
    }
}

sqlite_database::~sqlite_database() {
    //Close, but do not throw on failure, just print it:
    try {
        close();
    } catch (Exception & e) {
        cerr << "Exception in sqlite_database destructor: " << e.message << endl << "Ingoring." << endl;
    }
}

void sqlite_database::beginTransaction() {
    if(!transaction_active){
        exec("BEGIN;");
        transaction_active = true;
    }
}

void sqlite_database::endTransaction() {
    if (transaction_active)
        exec("END;");
    transaction_active = false;
}

void sqlite_database::exec(const string & query) {
    char * err = 0;
    sqlite3_exec(db, query.c_str(), 0, 0, &err);
    if (err != 0) {
        stringstream ss;
        ss << "sqlite_database::exec(\"" << query << "\") returned error: " << err;
        sqlite3_free(err);
        //database errors should not happen at all. If they do, we cannot use the log table, so write the error
        // to stderr, so the user knows what has happened:
        throw DatabaseException(ss.str());
    }
}

sqlite3_stmt* sqlite_database::prepare(const string & sql) {
    sqlite3_stmt * statement = 0;
    int ret = sqlite3_prepare_v2(db, sql.c_str(), sql.size() + 1, &statement, 0);
    if (ret != SQLITE_OK) {
        if (statement != 0) {
            sqlite3_finalize(statement);
        }
        try {
            error( __FUNCTION__);
        } catch (DatabaseException & ex) {
            stringstream ss;
            ss << ex.message << " (return code " << ret << "); SQL was " << sql;
            throw DatabaseException(ss.str());
        }
    }
    return statement;
}

void sqlite_database::error(const string & functionName) {
    stringstream ss;
    ss << "Error in function " << functionName << ": sqlite said '" << sqlite3_errmsg(db) << "'";
    throw DatabaseException(ss.str());
}

std::auto_ptr<Table> sqlite_database::create_table(const string & table_name){
    check_name(table_name);
    sqlite_table * result = new sqlite_table(table_name, boost::dynamic_pointer_cast<sqlite_database>(shared_from_this()));
    if(table_name == "products"){
        result->save_all_columns = save_all_products;
        result->save_columns = products_data;
    }
    return std::auto_ptr<Table>(result);
}


sqlite_database::sqlite_table::sqlite_table(const string & name_, const boost::shared_ptr<sqlite_database> & db_) : Table(db_),
    name(name_), table_created(false), next_insert_index(1), insert_statement(0), db(db_),
    save_all_columns(true) {
}

Column sqlite_database::sqlite_table::add_column(const std::string & name, const data_type & type){
    if(table_created) throw invalid_argument("sqlite_table::add_column called after table already created (via call to set_column / add_row).");
    if(!save_all_columns && save_columns.find(name) == save_columns.end()) return Column(-1);
    if(column_definitions.str().size() > 0)
        column_definitions << ", ";
    column_definitions << "'" << name << "' ";
    switch(type){
        case typeDouble:
            column_definitions << "DOUBLE";
            break;
        case typeInt: column_definitions << "INTEGER(4)"; break;
        case typeString: column_definitions << "TEXT"; break;
        case typeHisto: column_definitions << "BLOB"; break;
        default:
            throw invalid_argument("Table::add_column: invalid type parameter given.");
    };
    if(ss_insert_statement.str().size() > 0)
        ss_insert_statement << ", ";
    ss_insert_statement << "'" << name << "'";
    Column result(next_insert_index++);
    column_infos[result] = column_info(name, type);
    return result;
}


void sqlite_database::sqlite_table::create_table(){
    stringstream ss;
    string col_def = column_definitions.str();
    ss << "CREATE TABLE '" << name << "' (" << col_def << ");";
    db->exec(ss.str());
    
    ss.str("");
    ss << "INSERT INTO '" << name << "'(" << ss_insert_statement.str() << ") VALUES (";
    for(int i=1; i < next_insert_index; ++i){
        if(i==1)ss << "?";
        else ss << ", ?";
    }
    ss << ");";
    insert_statement = db->prepare(ss.str());
    table_created = true;
}

//create the table if it is empty to ensure that all tables have been created
// even if there are no entries
sqlite_database::sqlite_table::~sqlite_table(){
    if(not table_created) create_table();
}

//NOTE: If coding for multiple threads accessing a single sqlite database,
// we would have to do explicit locking between the insert statement execution
// and this function call to ensure that we get the correct id here.

void sqlite_database::sqlite_table::add_row(const Row & row){
    if(not table_created) create_table();
    for(map<Column, column_info>::const_iterator it=column_infos.begin(); it!=column_infos.end(); ++it){
        int index = it->first.get_id();
        if(index <= 0) continue;
        if(it->second.type == typeDouble){
            sqlite3_bind_double(insert_statement, index, row.get_column_double(it->first));
        }
        else if(it->second.type == typeInt){
            sqlite3_bind_int(insert_statement, index, row.get_column_int(it->first));
        }
        else if(it->second.type == typeString){
            const string & s = row.get_column_string(it->first);
            sqlite3_bind_text(insert_statement, index, s.c_str(), s.size(), SQLITE_TRANSIENT);
        }
        else if(it->second.type == typeHisto){
            const Histogram1D & h = row.get_column_histogram(it->first);
            boost::scoped_array<double> blob_data(new double[h.get_nbins()+4]);
            blob_data[0] = h.get_xmin();
            blob_data[1] = h.get_xmax();
            //set underflow to 0.0:
            blob_data[2] = 0.0;
            std::copy(h.get_data(), h.get_data() + h.size(), &blob_data[3]);
            //set overflow to 0.0:
            blob_data[h.get_nbins()+3] = 0.0;
            size_t nbytes = sizeof(double) * (h.get_nbins() + 4);
            sqlite3_bind_blob(insert_statement, index, &blob_data[0], nbytes, SQLITE_TRANSIENT);
        }
    }
    int res = sqlite3_step(insert_statement);
    sqlite3_reset(insert_statement);
    //reset all to NULL
    sqlite3_clear_bindings(insert_statement);
    if (res != 101) {
        db->error(__FUNCTION__);
    }
}

REGISTER_PLUGIN(sqlite_database)
