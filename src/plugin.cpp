#include "interface/plugin.hpp"

using namespace theta::plugin;
using namespace std;

int theta::plugin::plugin_build_depth=0;

void PluginLoader::execute(const Configuration & cfg) {
    SettingWrapper files = cfg.setting["plugin_files"];
    size_t n = files.size();
    for (size_t i = 0; i < n; i++) {
        std::string filename = cfg.replace_theta_dir(files[i]);
        load(filename);
    }
}

void PluginLoader::load(const std::string & soname) {
    void* handle = 0;
    try {
        handle = dlopen(soname.c_str(), RTLD_NOW | RTLD_GLOBAL);
    } catch (Exception & ex) {
        std::stringstream ss;
        ss << ex.message << " (in PluginLoader::load while loading plugin file '" << soname << "')";
        ex.message = ss.str();
        throw;
    }
    if (handle == 0) {
        std::stringstream s;
        const char * error = dlerror();
        if (error == 0) error = "0";
        s << "PluginLoader::load: error loading plugin file '" << soname << "': " << error << std::endl;
        throw InvalidArgumentException(s.str());
    }
}
