#include "plugins/blackhole_database.hpp"
#include "interface/plugin.hpp"
#include "interface/histogram.hpp"

#include <sstream>

#include <boost/filesystem.hpp>
#include <boost/scoped_array.hpp>

using namespace std;
using namespace theta;

blackhole_database::blackhole_database(const plugin::Configuration & cfg){}

blackhole_database::~blackhole_database() {}

std::auto_ptr<Table> blackhole_database::create_table(const string & table_name){
    check_name(table_name);
    return std::auto_ptr<Table>(new blackhole_table(boost::dynamic_pointer_cast<blackhole_database>(shared_from_this())));
}

blackhole_database::blackhole_table::~blackhole_table(){}

std::auto_ptr<Column> blackhole_database::blackhole_table::add_column(const std::string & name,
                                                                      const data_type & type){
    return std::auto_ptr<Column>(new blackhole_column());
}

void blackhole_database::blackhole_table::set_autoinc_column(const std::string &){
}

void blackhole_database::blackhole_table::set_column(const Column & c, double d){}

void blackhole_database::blackhole_table::set_column(const Column & c, int i){}

void blackhole_database::blackhole_table::set_column(const Column & c, const std::string & s){}

void blackhole_database::blackhole_table::set_column(const Column & c, const theta::Histogram & h){}

int blackhole_database::blackhole_table::add_row(){ return 0;}

REGISTER_PLUGIN(blackhole_database)
