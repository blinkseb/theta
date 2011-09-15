#ifndef PLUGINS_SQLITE_DATABASE_IN_HPP
#define PLUGINS_SQLITE_DATABASE_IN_HPP


#include "interface/decls.hpp"
#include "interface/database.hpp"

#include <sqlite3.h>


/** \brief Database class to read data from an sqlite file
 *
 * In particular suitable for reading data previously written with
 * \link sqlite_database \endlink.
 *
 * Configured via a setting like
 * \code
 *  input_database = {
 *     type = "sqlite_database_in";
 *     filename = "in.db";
 *     //alternatively: filenames = ("in1.db", "in2.db", "in3.db");
 *  };
 * \endcode
 *
 * Exactly one of the settings \c filename and \c filenames must be given.
 */
class sqlite_database_in: public theta::DatabaseInput{
public:
    /// Constructor used by the plugin system to build an instance given the configuration
    sqlite_database_in(const theta::plugin::Configuration & cfg);
    
    ~sqlite_database_in();
    
    virtual std::auto_ptr<ResultIterator> query(const std::string & table_name, const std::vector<std::string> & column_names);
    
private:
   class SqliteResultIterator: public DatabaseInput::ResultIterator{
   private:
       sqlite3 * db;
       sqlite3_stmt * statement;
       bool has_data_;
   public:
       SqliteResultIterator(sqlite3_stmt * st, sqlite3 * db_): db(db_), statement(st){
           operator++();
       }
       ~SqliteResultIterator(){
           sqlite3_finalize(statement);
       }
       bool has_data(){
           return has_data_;
       }
       virtual void operator++();
       double get_double(size_t icol);
       int get_int(size_t icol);
       std::string get_string(size_t icol);
       theta::Histogram get_histogram(size_t icol);
   };

   sqlite3 * db;
   size_t n_files;
};



#endif

