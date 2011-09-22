#include "interface/variables.hpp"
#include "interface/codegen.hpp"
#include "interface/pm.hpp"

#include <boost/shared_ptr.hpp>

#include <cstdlib>
#include <cstdio>
#include <dlfcn.h>

using namespace std;
using namespace theta;

size_t codegen::get_id(const ParId & id){
   return id.id;
}

ParId codegen::create_pid(size_t id){
   return ParId(id);
}


size_t codegen::get_id(const ObsId & id){
   return id.id;
}

ObsId codegen::create_oid(size_t id){
   return ObsId(id);
}

void codegen::header(std::ostream & out, const theta::PropertyMap & pm, const theta::ParIds& pids, const ObsIds & oids){
    out << "#include \"interface/variables.hpp\"" << endl;
    out << "#include \"interface/codegen.hpp\"" << endl;
    out << "#include \"interface/utils.hpp\"" << endl;
    out << "#include \"interface/phys.hpp\"" << endl << endl;
    out << "#include <string.h>" << endl;
    out << "using namespace theta::utils;" << endl;
    out << "using namespace theta;" << endl << endl;
    out << "namespace {" << endl;
    out << "ParIds get_all_par_ids(){" << endl;
    out << "    ParIds result;" << endl;
    for(theta::ParIds::const_iterator it=pids.begin(); it!=pids.end(); ++it){
        out << "    result.insert(codegen::create_pid(" << codegen::get_id(*it) << "));" << endl;
    }
    out << "    return result;" << endl;
    out << "}" << endl << endl;
    out << "const ParIds all_par_ids(get_all_par_ids());" << endl;
    size_t i=0;
    boost::shared_ptr<theta::VarIdManager> vm = pm.get<theta::VarIdManager>();
    for(theta::ParIds::const_iterator it=pids.begin(); it!=pids.end(); ++it, ++i){
        out << "const size_t pindex_" << vm->getName(*it) << " = " << i << ";" << endl;
    }
    out << "const size_t N_parameters = " << pids.size() << ";" << endl << endl;
    
    out << "ObsIds get_all_obs_ids(){" << endl;
    out << "    ObsIds result;" << endl;
    for(theta::ObsIds::const_iterator it=oids.begin(); it!=oids.end(); ++it){
        out << "    result.insert(codegen::create_oid(" << codegen::get_id(*it) << "));" << endl;
    }
    out << "    return result;" << endl;
    out << "}" << endl << endl;
    out << "const ObsIds all_obs_ids(get_all_obs_ids());" << endl;
    out << "const size_t N_observables = " << oids.size() << ";" << endl;
    out << "const size_t obs_offset[N_observables + 1] = {0";
    size_t offset = 0;
    for(ObsIds::const_iterator it=oids.begin(); it!=oids.end(); ++it){
        offset += vm->get_nbins(*it);
        if(offset % 2) ++offset;
        out << ", " << offset;
    }
    out << "};" << endl;
    
    out << "}" << endl << endl; // namespace {
    out << "extern \"C\" {" << endl;
}


void codegen::footer(std::ostream & out){
   out << "}" << endl;
}

string codegen::dtos(double d){
    char buffer[100];
    snprintf(buffer, 100, "%.25e", d);
    return buffer;
}

void* codegen::compile_load_so(const std::string & cpp){
    size_t p = cpp.find_last_of('/');
    if(p==string::npos) p=0;
    else ++p;
    std::string name = cpp.substr(p);
    stringstream ss;
    ss << utils::theta_dir  << "/lib/codegen_" << name << "_" << time(0) << ".so";
    string command = "g++ -rdynamic -shared -fPIC -O3 -I" + utils::theta_dir + " " + cpp + " -o " + ss.str();
    cout << command << endl;
    int ret = system(command.c_str());
    if(ret!=0){
        throw Exception("compilation failed (command was '" + command + "')");
    }
    void * result = dlopen(ss.str().c_str(), RTLD_NOW | RTLD_LOCAL);
    if(!result){
       stringstream ss_error;
       ss_error << "dlopen failed: " << dlerror();
       throw Exception(ss_error.str());
    }
    return result;
}


