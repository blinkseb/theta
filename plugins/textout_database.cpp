#include "plugins/textout_database.hpp"
#include "interface/plugin.hpp"
#include "interface/redirect_stdio.hpp"
#include "interface/histogram.hpp"

using namespace std;
using namespace theta;

textout_database::textout_database(const plugin::Configuration & cfg) : save_all_products(true){
    if(cfg.setting.exists("products_data")){
      try{
         string s = cfg.setting["products_data"];
         if(s=="*")save_all_products = true;
         else throw ConfigurationException("products_data setting is a string but not '*'");
      }
      catch(libconfig::SettingTypeException & e){
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
}

textout_database::~textout_database() {
}


std::auto_ptr<Table> textout_database::create_table(const string & table_name){
    check_name(table_name);
    textout_table * result = new textout_table(table_name, boost::dynamic_pointer_cast<textout_database>(shared_from_this()));
    if(table_name == "products"){
        result->save_all_columns = save_all_products;
        result->save_columns = products_data;
    }
    return std::auto_ptr<Table>(result);
}


textout_database::textout_table::textout_table(const string & name_, const boost::shared_ptr<textout_database> & db_) : Table(db_), irow(0),
    name(name_), have_autoinc(false), next_autoinc_value(1), db(db_), save_all_columns(true) {
}

std::auto_ptr<Column> textout_database::textout_table::add_column(const std::string & name, const data_type & type){
    if(!save_all_columns && save_columns.find(name) == save_columns.end()) return std::auto_ptr<Column>(new textout_column(-1));
    column_names.push_back(name);
    current_col_values.resize(current_col_values.size()+1);
    return std::auto_ptr<Column>(new textout_column(current_col_values.size() - 1));
}


void textout_database::textout_table::set_autoinc_column(const std::string & name){
    if(have_autoinc)
         throw InvalidArgumentException("textout_database::textout_table::set_autoinc_column: already have an autoinc column");
    have_autoinc = true;
}

void textout_database::textout_table::set_column(const Column & c, double d){
    int index = static_cast<const textout_column&>(c).index;
    if(index >= 0)
        current_col_values[index] = d;
}

void textout_database::textout_table::set_column(const Column & c, int i){
    int index = static_cast<const textout_column&>(c).index;
    if(index >= 0)
        current_col_values[index] = i;
}

void textout_database::textout_table::set_column(const Column & c, const std::string & s){
    int index = static_cast<const textout_column&>(c).index;
    if(index >= 0)
        current_col_values[index] = s;
}

void textout_database::textout_table::set_column(const Column & c, const theta::Histogram & h){
    int index = static_cast<const textout_column&>(c).index;
    if(index >= 0)
        current_col_values[index] = h;
}

ostream & operator<<(ostream & out, const theta::Histogram & h){
    out << "histogram(nbins=" << h.get_nbins() << ",xmin=" << h.get_xmin() << ",xmax=" << h.get_xmax() << ",data=[";
    for(size_t i=0; i<=h.get_nbins(); ++i){
        out << h.get(i) << ",";
    }
    return out << h.get(h.get_nbins()+1) << "])";
}

class printer: public boost::static_visitor<>{
private:
    ostream & out;
public:
    printer(ostream & out_): out(out_){}
    
    template <typename T>
    void operator()(const T& t){
        out << t;
    }
};

int textout_database::textout_table::add_row(){
    theta::cout << endl << "Table '" << name << "', row " << ++irow << ":" << endl;
    printer p(theta::cout);
    for(size_t i=0; i<column_names.size(); ++i){
        theta::cout << column_names[i] << "=";
        apply_visitor(p, current_col_values[i]);
        theta::cout << endl;
        //<< current_col_values[i] << endl;
    }
    if(have_autoinc){
        return next_autoinc_value++;
    }
    else{
        return 0;
    }
}

REGISTER_PLUGIN(textout_database)
