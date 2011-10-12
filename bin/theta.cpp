#include "interface/cfg-utils.hpp"
#include "interface/plugin.hpp"
#include "interface/histogram.hpp"
#include "interface/variables-utils.hpp"
#include "interface/variables.hpp"
#include "interface/main.hpp"
#include "interface/redirect_stdio.hpp"

#include "libconfig/libconfig.h++"

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include "boost/date_time/posix_time/posix_time_types.hpp"

#include <termios.h>

#include <fstream>

using namespace std;
using namespace theta;
using namespace theta::utils;

using namespace libconfig;

namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace btime = boost::posix_time;

class MyProgressListener: public ProgressListener{
public:
    virtual void progress(int done_, int total_, int n_error_){
        done = done_;
        total = total_;
        n_error = n_error_;
        print();
    }
    
    void print(){
        if(!is_tty) return;
        btime::ptime now = btime::microsec_clock::local_time();
        if(now < next_update && done < total) return;
        //move back to beginning of terminal line:
        theta::cout << "\033[" << chars_written << "D";
        chars_written = 0;  
        char c[200];
        double error_fraction = 100.0 * n_error / done;
        const char * color_start;
        const char * color_end = "\033[0m";
        if(error_fraction > 30){ //red:
           color_start = "\033[31m";
        }
        else if(error_fraction > 5){ //yellow:
           color_start = "\033[1;33m";
        }
        else{ //green:
           color_start = "\033[32m";
        }
        if(total > 0){
            double progress_fraction = 100.0 * done / total;
            chars_written += snprintf(c, 200, "progress: %6d / %-6d [%5.1f%%]   errors: %s%6d [%5.1f%%]%s", done, total, progress_fraction, color_start, n_error, error_fraction, color_end);
        }
        else{
            chars_written += snprintf(c, 200, "progress: %6d   errors: %s%6d [%5.1f%%]%s", done, color_start, n_error, error_fraction, color_end);
        }
        theta::cout << c;
        theta::cout.flush();
        next_update = now + btime::milliseconds(50);
    }

    MyProgressListener(): stdout_fd(theta::cout_fd), done(0), total(0), n_error(0), is_tty(isatty(stdout_fd)),
      chars_written(0), next_update(btime::microsec_clock::local_time()) {
        if(!is_tty) return;
        //disable terminal echoing; we don't expect any input.
        termios settings;
        if (tcgetattr (stdout_fd, &settings) < 0) return; //ignore error
        settings.c_lflag &= ~ECHO;
        tcsetattr (stdout_fd, TCSANOW, &settings);
    }
    
    ~MyProgressListener(){
        if(!is_tty) return;
        //enable terminal echoing again; don't be evil
        termios settings;
        if (tcgetattr (stdout_fd, &settings) < 0) return; //ignore error
        settings.c_lflag |= ECHO;
        tcsetattr (stdout_fd, TCSANOW, &settings);
        if(chars_written > 0)
            theta::cout << endl;
    }
    
private:
    int stdout_fd;
    int done, total, n_error;
    bool is_tty;
    int chars_written;
    btime::ptime next_update;
};



namespace{
    //note: these functions are useful to have compatibility
    // with both V2 and V3 of boost::filesystem
    // as in V2, path::filename returns a string whereas in
    // V3, path::filename returns a path.
    std::string to_string(const std::string & s){
        return s;
    }
    std::string to_string(const fs::path & p){
        return p.string();
    }
}


boost::shared_ptr<Main> build_main(string cfg_filename, bool nowarn){
    Config cfg;
    boost::shared_ptr<SettingUsageRecorder> rec(new SettingUsageRecorder());
    boost::shared_ptr<Main> main;
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    bool init_complete = false;
    try {
        try {
            //as includes in config files should always be resolved relative to the config file's location:
            string old_path = fs::current_path().string();
            //convert any failure to a FileIOException:
            try{
                 if(fs::path(cfg_filename).has_parent_path()){
                    fs::current_path(fs::path(cfg_filename).parent_path());
                 }
                 cfg_filename = to_string(fs::path(cfg_filename).filename());
            }
            catch(fs::filesystem_error & ex){
                 throw FileIOException();
            }
            cfg.readFile(cfg_filename.c_str());
            fs::current_path(old_path);
        } catch (FileIOException & f) {
            stringstream s;
            s << "Configuration file " << cfg_filename << " could not be read";
            throw ConfigurationException(s.str());
        } catch (ParseException & p) {
            stringstream s;
            s << "Error parsing configuration file: " << p.getError() << " in line " << p.getLine() << ", file " << p.getFile();
            throw ConfigurationException(s.str());
        }
        
        SettingWrapper root(cfg.getRoot(), cfg.getRoot(), rec);
        Configuration config(root);
        config.pm->set("default", vm);
        
        //process options:
        Configuration cfg_options(config, config.setting["options"]);
        PluginLoader::execute(cfg_options);
        
        //populate VarIdManager from config:
        apply_vm_settings(config);
        //build run:
        main = PluginManager<Main>::build(Configuration(config, root["main"]));
        init_complete = true;
    }
    catch (SettingNotFoundException & ex) {
        theta::cerr << "Error: the required setting " << ex.getPath() << " was not found." << endl;
    } catch (SettingTypeException & ex) {
        theta::cerr << "Error: the setting " << ex.getPath() << " has the wrong type." << endl;
    } catch (Exception & e) {
        theta::cerr << "Error: " << e.message << endl;
    }
    if(not init_complete){
        main.reset();
        return main;
    }
    
    if(not nowarn){
        vector<string> unused;
        rec->get_unused(unused, cfg.getRoot());
        if (unused.size() > 0) {
            theta::cout << "WARNING: following setting paths in the configuration file have not been used: " << endl;
            for (size_t i = 0; i < unused.size(); ++i) {
                theta::cout << "  " << (i+1) << ". " << unused[i] << endl;
            }
            theta::cout << "Comment out these settings to get rid of this message." << endl;
        }
    }
    return main;
}


int main(int argc, char** argv) {
    bool redirect_io;
    
    po::options_description desc("Options");
    desc.add_options()("help,h", "show help message")
    ("quiet,q", "quiet mode (suppress progress message)")
    ("nowarn", "do not warn about unused configuration file statements")
    ("redirect-io", po::value<bool>(&redirect_io)->default_value(true), "redirect stio/stderr of libraries to /dev/null");

    po::options_description hidden("Hidden options");

    hidden.add_options()
    ("cfg-file", po::value<vector<string> >(), "configuration file");

    po::positional_options_description p;
    p.add("cfg-file", -1);

    po::options_description cmdline_options;
    cmdline_options.add(desc).add(hidden);

    po::variables_map cmdline_vars;
    try{
        po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), cmdline_vars);
    }
    catch(std::exception & ex){
        theta::cerr << "Error parsing command line options: " << ex.what() << endl;
        return 1;
    }
    po::notify(cmdline_vars);
    
    if(cmdline_vars.count("help")){
        theta::cout << desc << endl;
        return 0;
    }
    
    if(cmdline_vars.count("cfg-file")==0){
        theta::cerr << "Error: you have to specify a configuration file" << endl;
        return 1;
    }
    
    if(redirect_io) theta::redirect_stdio();
    
    vector<string> cfg_filenames = cmdline_vars["cfg-file"].as<vector<string> >();
    bool quiet = cmdline_vars.count("quiet");
    bool nowarn = cmdline_vars.count("nowarn");
    
    //determine theta_dir (for config file replacements with $THETA_DIR
    fill_theta_dir(argv);
    if(theta_dir==""){
        theta::cerr << "WARNING: could not determine theta_dir, leaving empty" << endl;
    }
    
    try {
        for(size_t i=0; i<cfg_filenames.size(); ++i){
            if(!quiet and cfg_filenames.size() > 1){
                theta::cout << "processing file " << (i+1) << " of " << cfg_filenames.size() << ", " << cfg_filenames[i] << endl;
            }
            boost::shared_ptr<Main> main = build_main(cfg_filenames[i], nowarn);
            if(!main) return 1;
            if(not quiet){
                boost::shared_ptr<ProgressListener> l(new MyProgressListener());
                main->set_progress_listener(l);
            }

            //install signal handler now, not much earlier. Otherwise, plugin loading in build_run()
            // might change it ...
            install_sigint_handler();
            main->run();
            if(stop_execution) break;
        }
    }
    catch(ExitException & ex){
       theta::cerr << "Exit requested: " << ex.message << endl;
       return 1;
    }
    catch (Exception & ex) {
        theta::cerr << "An error ocurred in Main::run: " << ex.what() << endl;
        return 1;
    }
    catch(FatalException & ex){
        theta::cerr << "A fatal error ocurred in Main::run: " << ex.message << endl;
        return 1;
    }
    if(theta::stop_execution){
        theta::cout << "(exiting on SIGINT)" << endl;
    }
    return 0;
}
