#ifndef PLUGIN_BLACKHOLE_DATABASE_HPP
#define PLUGIN_BLACKHOLE_DATABASE_HPP

#include "interface/decls.hpp"
#include "interface/database.hpp"

#include <memory>
#include <string>


/** \brief Database which discards all information
 *
 * This is only useful for testing and benchmarking purposes.
 *
 * Configured via a setting group like
 * \code
 * {
 *   type = "blackhole_database";
 * }
 * \endcode
 *
 * \c type must always be "blackhole_database" in order to select this plugin
 */
class blackhole_database: public theta::Database{
public:
    
    /** \brief Constructor for the plugin system
     */
    blackhole_database(const theta::plugin::Configuration & cfg);
    
    virtual ~blackhole_database();
    
    /** \brief See documentation of Database::create_table
     */
    virtual std::auto_ptr<theta::Table> create_table(const std::string & table_name);
        
private:
    //declare privately(!) the sqlite_table class:
    class blackhole_table: public theta::Table {
        friend class blackhole_database;
        blackhole_table(const boost::shared_ptr<Database> & db): Table(db){}
        
        virtual ~blackhole_table();
        
        std::auto_ptr<theta::Column> add_column(const std::string & name, const theta::data_type & type);
        virtual void set_autoinc_column(const std::string & name);
        virtual void set_column(const theta::Column & c, double d);
        virtual void set_column(const theta::Column & c, int i);
        virtual void set_column(const theta::Column & c, const std::string & s);
        virtual void set_column(const theta::Column & c, const theta::Histogram & h);
        virtual int add_row();
    
        class blackhole_column: public theta::Column{
        public:
            virtual ~blackhole_column(){}
        };
    };
};

#endif
