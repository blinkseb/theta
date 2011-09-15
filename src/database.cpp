#include <string>
#include <sstream>

#include <boost/date_time/local_time/local_time.hpp>

#include "interface/database.hpp"
#include "interface/histogram.hpp"

using namespace std;
using namespace theta;

REGISTER_PLUGIN_BASETYPE(theta::Database);
REGISTER_PLUGIN_BASETYPE(theta::DatabaseInput);

void Database::check_name(const string & column_name) {
    if (column_name.size() == 0)
        throw DatabaseException("Database::check_name: name was empty.");
    if (not ((column_name[0] >= 'A' && column_name[0] <= 'Z') or (column_name[0] >= 'a' && column_name[0] <= 'z'))) {
        stringstream ss;
        ss << "Database::check_name: '" << column_name << "' does not start with a letter as is required.";
        throw DatabaseException(ss.str());
    }
    for (size_t i = 1; i < column_name.size(); i++) {
        if (column_name[i] >= 'A' && column_name[i] <= 'Z')
            continue;
        if (column_name[i] >= 'a' && column_name[i] <= 'z')
            continue;
        if (column_name[i] >= '0' && column_name[i] <= '9')
            continue;
        if (column_name[i] == '_')
            continue;
        stringstream ss;
        ss << "Database::check_name: name '" << column_name << "' invalid at letter " << i;
        throw DatabaseException(ss.str());
    }
}
// place destructors here, not in header. Otherwise some gcc versions are missing the vtable which is placed
// in the module defining the first non-inline function ...
Table::~Table(){}

Column::~Column(){}

/* ProductsTable */
ProductsTable::ProductsTable(std::auto_ptr<Table> & table_): table(table_){
    c_runid = table->add_column("runid", typeInt);
    c_eventid = table->add_column("eventid", typeInt);
}

std::auto_ptr<Column> ProductsTable::declare_product(const ProductsSource & source, const std::string & column_name, const data_type & type){
    std::string new_name = source.getName() + "__" + column_name;
    return table->add_column(new_name, type);
}

void ProductsTable::add_row(int runid, int eventid){
    table->set_column(*c_runid, runid);
    table->set_column(*c_eventid, eventid);
    table->add_row();
}

/* LogTable */
LogTable::LogTable(std::auto_ptr<Table> & table_): level(info), table(table_){
   for(int i=0; i<4; ++i){
       n_messages[i]=0;
   }
   c_runid = table->add_column("runid", typeInt);
   c_eventid = table->add_column("eventid", typeInt);
   c_severity = table->add_column("severity", typeInt);
   c_message = table->add_column("message", typeString);
   c_time = table->add_column("time", typeDouble);
}

const int* LogTable::get_n_messages() const{
    return n_messages;
}
    

void LogTable::set_loglevel(e_severity s){
    level = s;
}

LogTable::e_severity LogTable::get_loglevel() const{
    return level;
}

void LogTable::really_append(int runid, int eventid, e_severity s, const string & message) {
    n_messages[s]++;
    table->set_column(*c_runid, runid);
    table->set_column(*c_eventid, eventid);
    table->set_column(*c_severity, s);
    table->set_column(*c_message, message);
    using namespace boost::posix_time;
    using namespace boost::gregorian;
    ptime t(microsec_clock::universal_time());
    time_duration td = t - ptime(date(1970, 1, 1));
    double time = td.total_microseconds() / 1000000.0;
    table->set_column(*c_time, time);
    table->add_row();
}

//RndInfoTable
RndInfoTable::RndInfoTable(std::auto_ptr<Table> & table_): table(table_){
    c_runid = table->add_column("runid", typeInt);
    c_name = table->add_column("name", typeString);
    c_seed = table->add_column("seed", typeInt);
}

void RndInfoTable::append(int runid, const string & name, int seed){
    table->set_column(*c_runid, runid);
    table->set_column(*c_name, name);
    table->set_column(*c_seed, seed);
    table->add_row();
}

