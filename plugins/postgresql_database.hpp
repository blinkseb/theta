#ifndef PLUGIN_POSTGRESQL_DATABASE_HPP
#define PLUGIN_POSTGRESQL_DATABASE_HPP

#include "interface/decls.hpp"
#include "interface/database.hpp"

#include "libpq-fe.h"
#include <memory>
#include <string>


/** \brief Database storing all data in a postgresql database
 *
 * WARNING: This plugin has not been extensively tested. While it should work in general, the
 * Table::typeHisto column type has not been tested extensively.
 *
 * This plugin is only included if you explicitely enable postgresql within cmake.
 *
 * Configured via a setting group like
 * \code
 * my_db = {
 *   type = "postgresql_database";
 *   conninfo = "host=some.host port=10815 user=abc password=def dbname=theta-1";
 *   table_prefix = "theta5"; //optional. Default is ""
 *   mode = "append"; // optional. Default is "recreate"
 * };
 * \endcode
 *
 * \c type must always be "postgresql_database" in order to select this plugin
 *
 * \c conninfo is the connection information used to create the connection to the postgresql database.
 *    See the <a href="http://www.postgresql.org/docs/8.4/static/libpq-connect.html">postgresql documentation
 *    for PQconnectdb</a> for details.
 *
 * \c table_prefix is a string which is prepended to the table names in order to avoid name colissions.
 *      This setting is optional; as default, no prefix is prepended.
 *
 * \c mode is either "append" or "recreate". In "recreate" mode, any pre-existing table of the same name 
 *     will be dropped before the tables are created. In case of "append", data is appended to a pre-existing table.
 *     Note that this can cause errors if the table format of the pre-existing table is not compatible to the one to be created.
 *
 * The types Table::typeDouble, Table::typeInt and Table::typeString are translated directly
 * to their SQL counterparts \c DOUBLE \c PRECISION, \c INT(4) and \c TEXT, respectively. For Table::typeHisto,
 * an SQL BLOB is saved which contains the lower and upper border of the histogram and the raw histogram data,
 * including underflow and overflow bin (see theta::Histogram::getData).
 */
class postgresql_database: public theta::Database{
public:
    
    /** \brief Constructor for the plugin system
     *
     * See class documentation for a description of the parsed COnfiguration settings.
     */
    postgresql_database(const theta::plugin::Configuration & cfg);
    
    
    virtual ~postgresql_database();
    
    /** \brief See documentation of Database::create_table
     */
    virtual std::auto_ptr<theta::Table> create_table(const std::string & table_name);
        
private:
    /** Execute the sql string \c query.
     */
    void exec(const std::string & query);
    
    /** Prepare the sql statement \c query.
     */
    void prepare(const std::string & query, int ncol, const std::string & statement_name);
    
    /** Start a Transaction. In case of an error, a DatabaseException is thrown.
     */
    void beginTransaction();
    
    /** End a transaction previously started with \c beginTransaction. Calling
     *  endTransaction without a call to beginTransaction is valid (it is a no-op).
     * */
    void endTransaction();

    /** Throw a DatabaseException using the last database error message.
     * 
     * This method is intended for Table objects which want to propagate
     * an error via Exceptions in a consistent manner.
     * 
     * \param functionName: name of the function (in a \c Table object) in which the database error ocurred.
     */
    void error(const std::string & functionName);
    
    void close();
    
    PGconn * conn;
    bool transaction_active;
    std::string table_prefix;
    std::string mode;
    
    
    //declare privately(!) the sqlite_table class:
    class postgresql_table: public theta::Table {
    friend class postgresql_database;

        // destructor; creates the table if empty
        virtual ~postgresql_table();
        
        std::auto_ptr<theta::Column> add_column(const std::string & name, const data_type & type);
        virtual void set_column(const theta::Column & c, double d);
        virtual void set_column(const theta::Column & c, int i);
        virtual void set_column(const theta::Column & c, const std::string & s);
        virtual void set_column(const theta::Column & c, const theta::Histogram & h);
        virtual void set_autoinc_column(const std::string & name);
        virtual int add_row();

    private:
        
        postgresql_table(const std::string & name_, postgresql_database * db_);
        bool check_existing_table();
        
        std::string name;
        bool have_autoinc;
        std::string autoinc_name;
        std::stringstream column_definitions; // use by the add_column method
        std::stringstream ss_insert_statement;
        
        //for comparison with already existing tables:
        std::vector<std::string> column_names;
        bool table_created;
        
        int next_column; // next free column to return by add_column
        
        std::string insert_statement;
        postgresql_database * db;
        
        std::vector<char *> column_content; // set by set_column
        std::vector<int> param_lengths;
        std::vector<int> param_formats;
        
        void create_table();
        
        class postgresql_column: public theta::Column{
        public:
            //starting at 1
            int column_index;
            postgresql_column(int i):column_index(i){}
            virtual ~postgresql_column(){}
        };
    };
};

#endif
