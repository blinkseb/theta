#include "root/rootfile_database.hpp"
#include "interface/plugin.hpp"
#include "interface/histogram.hpp"

#include <sstream>

#include <boost/filesystem.hpp>
#include <boost/scoped_array.hpp>

using namespace std;
using namespace theta;

rootfile_database::rootfile_database(const plugin::Configuration & cfg): file(0){
   std::string filename = cfg.replace_theta_dir(cfg.setting["filename"]);
   file = new TFile(filename.c_str(), "recreate");
   if(not file->IsOpen()){
       delete file;
       throw ConfigurationException("error opening output root file '" + filename + "'");
   }
   
   if(cfg.setting.exists("products_data")){
      try{
         string s = cfg.setting["products_data"];
         if(s=="*")save_all_products = true;
         else throw ConfigurationException("products_data setting is a string but not '*'");
      }
      catch(libconfig::SettingTypeException & e){
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
   else{
      save_all_products = true;
   }
   
   if(cfg.setting.exists("products_histograms")){
      size_t n = cfg.setting["products_histograms"].size();
      for(size_t i=0; i<n; ++i){
          int nbins = cfg.setting["products_histograms"][i]["nbins"];
          double xmin = cfg.setting["products_histograms"][i]["range"][0];
          double xmax = cfg.setting["products_histograms"][i]["range"][1];
          hist_infos.push_back(hist_info());
          hist_infos.back().h.reset(nbins, xmin, xmax);
          hist_infos.back().name = static_cast<string>(cfg.setting["products_histograms"][i]["name"]);
          hist_infos.back().column_name = static_cast<string>(cfg.setting["products_histograms"][i]["column"]);
      }
   }
}

rootfile_database::~rootfile_database() {
    if(file){
       //write those root histos:
       TDirectory * dir = file->mkdir("products_histograms");
       for(size_t i=0; i<hist_infos.size(); ++i){
           TH1D * root_histo = new TH1D(hist_infos[i].name.c_str(), hist_infos[i].name.c_str(), hist_infos[i].h.get_nbins(),
                                        hist_infos[i].h.get_xmin(), hist_infos[i].h.get_xmax());
           root_histo->SetContent(hist_infos[i].h.getData());
           root_histo->SetDirectory(dir);
       }
       file->Write();
       file->Close();
       delete file;
       file = 0;
    }
}

std::auto_ptr<Table> rootfile_database::create_table(const string & table_name){
    check_name(table_name);
    rootfile_table * result = new rootfile_table(table_name, boost::dynamic_pointer_cast<rootfile_database>(shared_from_this()));
    if(table_name == "products"){
        result->save_all_columns = save_all_products;
        result->save_columns = products_data;
    }
    return std::auto_ptr<Table>(result);
}


rootfile_database::rootfile_table::rootfile_table(const std::string & tablename, const boost::shared_ptr<rootfile_database> & db_):Table(db_),
   db(db_), save_all_columns(true){
       tree = new TTree(tablename.c_str(), tablename.c_str());
       tree->SetDirectory(db->file);
       //ugly switch, but should work:
       products_table = tablename == "products";
}

rootfile_database::rootfile_table::~rootfile_table(){}

std::auto_ptr<Column> rootfile_database::rootfile_table::add_column(const std::string & name, const data_type & type){
    std::auto_ptr<Column> result;
    bool make_branch = false;
    if(save_all_columns || save_columns.find(name)!=save_columns.end()) make_branch = true;
    switch(type){
        case theta::typeDouble:
            result.reset(new rootfile_column_double(name));
            if(make_branch)
                tree->Branch(name.c_str(), &static_cast<rootfile_column_double*>(result.get())->d, "data/D");
            break;
        case theta::typeInt:
            result.reset(new rootfile_column_int(name));
            if(make_branch)
               tree->Branch(name.c_str(), &static_cast<rootfile_column_int*>(result.get())->i, "data/I");
            break;
        case theta::typeString:
            result.reset(new rootfile_column_string(name));
            if(make_branch)
               tree->Branch(name.c_str(), "TString", &static_cast<rootfile_column_string*>(result.get())->s);
            break;
        case theta::typeHisto:
            result.reset(new rootfile_column_histo(name));
            if(make_branch)
               tree->Branch(name.c_str(), "TH1D", &static_cast<rootfile_column_histo*>(result.get())->h);
            break;
    }
    return result;
}

void rootfile_database::rootfile_table::set_autoinc_column(const std::string & s){
    std::auto_ptr<Column> col = add_column(s, theta::typeInt);
    autoinc_column.reset(static_cast<rootfile_column_int*>(col.release()));
}

void rootfile_database::rootfile_table::set_column(const Column & c_, double d){
   const rootfile_column_double& c = static_cast<const rootfile_column_double&> (c_);
   c.d = d;
   if(products_table){
       for(size_t i=0; i< db->hist_infos.size(); ++i){
           if(db->hist_infos[i].column_name == c.name){
               db->hist_infos[i].h.fill(d, 1.0);
           }
       }
   }
}

void rootfile_database::rootfile_table::set_column(const Column & c, int i){
   (static_cast<const rootfile_column_int&> (c)).i = i;
}

void rootfile_database::rootfile_table::set_column(const Column & c, const std::string & s){
   *((static_cast<const rootfile_column_string&> (c)).s) = s.c_str();
}

void rootfile_database::rootfile_table::set_column(const Column & c, const theta::Histogram & h){
   TH1D * root_histo = (static_cast<const rootfile_column_histo&> (c)).h;
   root_histo->SetBins(h.get_nbins(), h.get_xmin(), h.get_xmax());
   root_histo->SetContent(h.getData());
}

int rootfile_database::rootfile_table::add_row(){
    int result = 0;
    if(autoinc_column.get()){
        result = ++autoinc_column->i;
    }
    tree->Fill();
    return result;

}

REGISTER_PLUGIN(rootfile_database)


