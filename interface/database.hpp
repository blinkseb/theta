#ifndef DATABASE_HPP
#define DATABASE_HPP

#include <string>
#include <memory>

#include <boost/utility.hpp>
#include <boost/enable_shared_from_this.hpp>

#include "interface/decls.hpp"
#include "interface/plugin.hpp"
#include "interface/data_type.hpp"
#include "interface/producer.hpp"

namespace theta {


/** \brief Abstract class for data input
 *
 * A database in theta is a collection of tables. Tables are identified by a unique name
 * within the database and can contain an arbitrary number of columns (identified by name) of type double, int, string and
 * Histogram.
 *
 * This is an abstract class. Concrete instances can be retrieved via the plugin system using
 * this class a template argument to PluginManager.
 */
class DatabaseInput: private boost::noncopyable{
public:

    /// Required for the plugin system
    typedef DatabaseInput base_type;
    
    /** \brief Iterator-like class to step through the result of a Database query
     */
    class ResultIterator: private boost::noncopyable{
    public:
        /** \brief Retrieve next row in the table
         *
         * If there is no next row, a subsequent call to has_data() will return false.
         *
         * In case of an error, a DatabaseException is thrown.
         */
        virtual void operator++() = 0;
        
        /** \brief Returns true iff this iterator points to a valid column containing data
         *
         * Returns false if there are no more result columns
         */
        virtual bool has_data() = 0;
        
        //@{
        /** \brief Retrieve the column values of the current row
         *
         * Before calling any of these, make sure has_data() returns true.
         *
         * The column index argument is zero-based and refers to the column_names vector
         * passed to DatabaseInput::query.
         *
         * In case of type mismatch, a DatabaseException is thrown.
         */
        virtual double get_double(size_t icol) = 0;
        virtual int get_int(size_t icol) = 0;
        virtual theta::Histogram get_histogram(size_t icol) = 0;
        virtual std::string get_string(size_t icol) = 0;
        //@}
        
        virtual ~ResultIterator(){}
    };

    /** \brief Select some columns from a table
     *
     * the returned ResultIterator points to the first result row (if any), i.e.,
     * call has_data and get_* on it first to get the result of the first row and *then* ResultIterator::operator++().
     *
     * In case of an error, a DatabaseException is thrown.
     */
    virtual std::auto_ptr<ResultIterator> query(const std::string & table_name, const std::vector<std::string> & column_names) = 0;
    virtual ~DatabaseInput(){}

};


/** \brief Abstract database class for data output
 *
 * A database in theta is a collection of tables. Tables are identified by a unique name
 * within the database and can contain an arbitrary number of columns (identified by name) of type double, int, string and
 * Histogram.
 *
 * This is an abstract class. Concrete instances can be retrieved via the plugin system using
 * this class a template argument to PluginManager.
 */
class Database: private boost::noncopyable, public boost::enable_shared_from_this<Database> {
public:
    
    /// Required for the plugin system
    typedef Database base_type;

    /** \brief Virtual Destructor
     *
     * Virtual, as polymrophic access to derived classes will happen.
     *
     * Should do any cleanup work (like closing the database file for file-based databases,
     * closing the network connection for network-based, etc.).
     */
    virtual ~Database(){}
    
    
    /** \brief Create a new table instance within this database
     *
     * Table names must start with a letter and otherwise consist only of letters,
     * digits and underscores. Throws an InvalidArgumentException if the name does not meet the requirements.
     *
     * The returned Table can be used only as long as this Database instance is not destroyed.
     * Using Tables created from a destroyed Database yields undefined behaviour.
     */
    virtual std::auto_ptr<Table> create_table(const std::string & table_name) = 0;

protected:
    
    /** \brief Checks the table name requirements
     *
     * Derived classes can call this at the beginning of their implementation of create_table
     * in order to meet the specification: this method will throw an InvalidArgumentException
     * if \c name violates the specification.
     */
    void check_name(const std::string & table_name);
};



/** \brief Abstract class for a table in a Database
 *
 * Tables are always constrcuted via a Datase instance. Once created, it can be used by:
 * <ol>
 *   <li>calling add_column one or more times at initialization and saving the result
 *       in some data member of the derived class</li>
 *   <li>To write out a row, call the set_column methods on all defined columns,
 *       followed by a call to add_row</li>
 * </ol>
 *
 * After the first call of set_column or add_row, to more calls to add_column are allowed.
 *
 * To meet this specification, derived classes will usually defer the actual table creation in the underlying
 * implementation until the first call of set_column or add_row. Note that even if no rows are added via add_row,
 * the table must still be created.
 */
class Table: private boost::noncopyable {
public:

    /// destructor; creates the table if empty
    virtual ~Table();
    
    /** \brief Add a column to the table
     *
     * The column will have name \c name and a data type scpecified by \c type.
     *
     * The result can be used as first argument to the set_column() methods. Other than saving it
     * and using it in \c set_column, you should not do anything with it.
     *
     * This method should only be called at the beginning, after the construction of the object and
     * before using it to actually write any data. If this method is called after a call
     * to either set_column or add_row, an IllegalStateException will be thrown.
     */
    virtual std::auto_ptr<Column> add_column(const std::string & name, const data_type & type) = 0;
    
    /** \brief Set an autoincrement column to this table
     *
     * Can be called at most once per table, i.e., a table can have at most one
     * autoinc column. If called more than once, an InvalidArgumentException is thrown.
     *
     * It is only valid to call this method if it is valid to call add_column.
     */
    virtual void set_autoinc_column(const std::string & name) = 0;
    
    //@{
    /** \brief Set column contents
     *
     * The first argument must be a column as returned by add_column.
     *
     * After setting all columns, add_row() has to be called to actually write
     * a row to the result table.
     */
    virtual void set_column(const Column & c, double d) = 0;
    virtual void set_column(const Column & c, int i) = 0;
    virtual void set_column(const Column & c, const std::string & s) = 0;
    virtual void set_column(const Column & c, const theta::Histogram & histo) = 0;
    //@}
    
    
    /** \brief Add a row to the table
     *
     * This uses the column values previously set with the set_column methods.
     * If set_column was not called for all columns in the table, the contents of these columns
     * will have a NULL-like value which depends on the particular implementation.
     *
     * If the table has an autoinc_column, i.e., if set_autoinc_column has been called,
     * the added id is returned. Otherwise, the result will always be 0.
     */
    virtual int add_row() = 0;
protected:
    //the tables always hold a shared ptr to the database to prevent
    // the database being destroyed earlier than all its tables (!)
    boost::shared_ptr<Database> db;
    Table(const boost::shared_ptr<Database> & db_): db(db_){}
private:
    Table(); //not implemented
};

/** \brief Base class for columns managed by Tables
 */
class Column{
public:
    virtual ~Column() = 0;
};


/** \brief A Table to store products, i.e., per-event output
 *
 * Per \link theta::Run Run \endlink, there is exactly one ProductsTable. Products are produced by
 * <ul>
 *   <li>Producer instances to save the (per-event) result of their computation</li>
 *   <li>DataSource to save some per-event information about data production</li>
 * </ul>
 *
 * An ProductsTable will have an integer column named "runid" and an integer column named "eventid". Additionally,
 * any columns defined by the Producer or DataSource, with column names as described in ProductsTable::add_column.
 *
 * Clients use ProductsTable very similarly to a Table. Differences are (i) the
 * signature of the add_column method which takes an additional \c name argument and (ii)
 * the add_row column, which takes a Run instance as argument here.
 *
 * The actual write is done by the \link Run \endlink instance; the Producer / DataSource must not
 * call add_row directly.
 */
class ProductsTable: public ProductsSink {
public:
        //@{
        /** \brief Forwards to Table::set_column
         */
        virtual void set_product(const Column & c, double d){
            table->set_column(c, d);
        }
        
        virtual void set_product(const Column & c, int i){
            table->set_column(c, i);
        }
        
        virtual void set_product(const Column & c, const std::string & s){
            table->set_column(c, s);
        }
        
        virtual void set_product(const Column & c, const theta::Histogram & histo){
            table->set_column(c, histo);
        }
        
        //@}
        
        /** \brief Add a column to this table
         *
         * Similar to Table::add_column, but expects an additional name of the caller. The caller class
         * must derive from ProducerTableWriter and pass themselves as first argument.
         *
         * The actual column name used in the table will be
         * \code
         *   tw.getName() + "__" + column_name
         * \endcode
         */
        virtual std::auto_ptr<Column> declare_product(const theta::ProductsSource & source, const std::string & product_name, const data_type & type);
        
        /** \brief Construct a new ProductsTable based on the given table
         *
         * Ownership of object held by table will be transferred.
         */
        ProductsTable(std::auto_ptr<Table> & table);
        
        /** \brief Add a row to the table, given the current run
         *
         * This is called by theta::Run after a producer has executed.
         */
        void add_row(int runid, int eventid);
        
private:
    std::auto_ptr<Table> table;
    std::auto_ptr<Column> c_runid, c_eventid;
};



/** \brief Table to store all logging information
 *
 * The corresponding table has following columns:
 * <ol>
 * <li>runid (typeInt)</li>
 * <li>eventid (typeInt)</li>
 * <li>severity (typeInt)</li>
 * <li>message (typeString)</li>
 * <li>time (typeDouble)</li>
 * </ol>
 * 
 * \c runid and \c eventid are the run and event id, respectively, the log entry
 * is associated to. \c eventid is set to 0, if no particular event (but the run as a whole)
 * is referred to.
 * 
 * \c severity is one level from LogTable::e_severity. See there for a description of the meaning
 * of these levels.
 * 
 * \c message is the human-readable log message.
 * 
 * \c time is the number of seconds since the unix epoch (1970-01-01 UTC),
 *    with sub-second accuracy.
 **/
class LogTable: private boost::noncopyable {
public:
    /** \brief Severity levels for log messages
     * 
     * <ol>
     * <li>\c error should be used if a serious condition is reported which will likely
     *     affect the whole result.</li>
     * <li>\c warning should be used if a problem is not as serious as \c error, but might affect
     *     the result in an undesired way. </li>
     * <li>\c info is used purely informational, without the implicit action request of
     *    \c error and \c warning. It covers such events like the start and end of a pseudo experiment. </li>
     * <li>\c debug is used like \c info, but for very detailed reporting which is usually not required.</li>
     * </ol>
     **/
    enum e_severity {
        error = 0, warning, info, debug
    };
    
    /** \brief Construct a new logTable based on the given Table
     *
     * The default loglevel is warning.
     *
     * Ownership of table will be transferred.
     */
    LogTable(std::auto_ptr<Table> & table);
    
    /** \brief Set the current log level.
     *
     * Future calls to append() will only write the message to the table if the messages
     * has a severity which is euqal to or exceeds the level given here.
     *
     * Note that is it not possible to disable logging of error messages.
     */
    void set_loglevel(e_severity s);
    
    /** \brief Get the currently active log level
     */
    e_severity get_loglevel() const;
    
    /** \brief Append message to log table, if severity is larger than currently configured level
     * 
     *
     * \c s is the severity level of the log message
     *
     * \c message is the message, in human-readable english
     */
    void append(int runid, int eventid, e_severity s, const std::string & message){
        //define inline as hint to the compiler; keep this function as short as possible to
        // encourage the compiler actually inlining it and to have high performance benefits in
        // case the user disables logging.
        if(s <= level) really_append(runid, eventid, s, message);
    }
    
    /** \brief Returns the number of messages
     *
     * Use severity::e_severity converted to int to access the number of messages written to the log.
     * Messages suppressed by the set_loglevel() method <em>not</em> counted.
     */
    const int* get_n_messages() const;

private:

    /// really append the log message. Called from append() in case severity is large enough
    void really_append(int runid, int eventid, e_severity s, const std::string & message);
    e_severity level;
    int n_messages[4];
    std::auto_ptr<Column> c_runid, c_eventid, c_severity, c_message, c_time;
    std::auto_ptr<Table> table;
};


/** \brief Table to store information about the random number seeds.
 *
 * There is a one-to-one relationship between this class and \link theta::Run Run \endlink, i.e.,
 * to each RndInfoTable instance, there is one Run instance and vice versa.
 *
 * The corresponding table has following columns:
 * <ol>
 * <li>runid (typeInt): the run id the entry refers to.</li>
 * <li>seed (typeInt): the seed of the random number generator used in the run.</li>
 * </ol>
 */
class RndInfoTable: private boost::noncopyable {
public:
    /** \brief Construct with an underlying Table
     *
     * Ownership of table will be transferred.
     */
    RndInfoTable(std::auto_ptr<Table> & table);

    /** \brief append an entry to the RndInfoTable
     *
     * \c runid is the current runid
     * \c name is the name of the module of the seed, according to \link theta::ProductsTableWriter ProductsTableWriter \endlink . <br />
     * \c seed is the seed of the random number generator used for this module
     */
    void append(int runid, const std::string & name, int seed);
private:
    std::auto_ptr<Column> c_runid, c_name, c_seed;
    std::auto_ptr<Table> table;
};

}

#endif

