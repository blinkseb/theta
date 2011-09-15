#include "interface/plugin.hpp"
#include "interface/phys.hpp"
#include "interface/utils.hpp"

#include "test/utils.hpp"

#include <boost/test/unit_test.hpp>

using namespace theta;
using namespace theta::plugin;
using namespace std;


BOOST_AUTO_TEST_SUITE(sysrate_tests)

BOOST_AUTO_TEST_CASE(sysrate0){
    BOOST_CHECKPOINT("sysrate0 entry");
    load_core_plugins();
    
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    ParId delta1 = vm->createParId("delta1");
    ParId delta2 = vm->createParId("delta2");
    ParId beta1 = vm->createParId("beta1");
    ParId beta2 = vm->createParId("beta2");
    
    //boost::ptr_vector<Function> coeffs;
    BOOST_CHECKPOINT("parsing config");
    
    ConfigCreator cc("f = {type = \"sys_rate_function\";"
         "factors = (\"beta1\", \"beta2\");\n"
         "sys_rates = ((\"delta1\", -0.2, 0.17), (\"delta2\", -0.9, 0.0));\n"
         "};", vm);
    
    BOOST_CHECKPOINT("config parsed");
    
    const theta::plugin::Configuration & cfg = cc.get();
    BOOST_CHECKPOINT("building function");
    std::auto_ptr<Function> f;
    try{
        f = PluginManager<Function>::instance().build(Configuration(cfg, cfg.setting["f"]));
    }
    catch(Exception & ex){
        cerr << ex.message << endl;
    }
    catch(libconfig::SettingNotFoundException & ex){
        cerr << ex.getPath() << " not found" << endl;
    }
    catch(libconfig::SettingTypeException & ex){
        cerr << ex.getPath() << " has wrong type" << endl;
    }
    BOOST_REQUIRE(f.get()!=0);
    
    ParValues values;
    //check behaviour for delta=0, i.e. "no systematics":
    values.set(beta1, 1.0);
    values.set(beta2, 1.0);
    values.set(delta1, 0.0);
    values.set(delta2, 0.0);
    BOOST_REQUIRE(utils::close_to_relative((*f)(values), 1.0));
    values.set(beta1, 1.7);
    values.set(beta2, 2.1);
    BOOST_REQUIRE(utils::close_to_relative((*f)(values), 1.7*2.1));
    //no truncation must occur here:
    values.set(beta1, 1.5e-7);
    values.set(beta2, 3.15e-8);
    BOOST_REQUIRE(utils::close_to_relative((*f)(values), 1.5e-7 * 3.15e-8));
    //vary one delta:
    values.set(beta1, 1.3).set(beta2, 0.88).set(delta1, 1.0);
    BOOST_REQUIRE(utils::close_to_relative((*f)(values), values.get(beta1) * values.get(beta2) * (1 + 0.17)));
    values.set(delta1, -1.0);
    BOOST_REQUIRE(utils::close_to_relative((*f)(values), values.get(beta1) * values.get(beta2) * (1 - 0.2)));
    values.set(delta1, -0.5);
    BOOST_REQUIRE(utils::close_to_relative((*f)(values), values.get(beta1) * values.get(beta2) * (1 - 0.5 * 0.2)));
    values.set(delta1, 2.5);
    BOOST_REQUIRE(utils::close_to_relative((*f)(values), values.get(beta1) * values.get(beta2) * (1 + 2.5 * 0.17)));
    //vary both deltas:
    values.set(delta2, 0.5);
    BOOST_REQUIRE(utils::close_to_relative((*f)(values), values.get(beta1) * values.get(beta2) * (1 + 2.5 * 0.17) * (1 + 0.5 * 0.0)));
    values.set(delta2, -0.88);
    BOOST_REQUIRE(utils::close_to_relative((*f)(values), values.get(beta1) * values.get(beta2) * (1 + 2.5 * 0.17) * (1 + 0.88 * -.9)));
    //set to an extreme value to test truncation:
    values.set(delta2, -2.0);
    BOOST_REQUIRE(utils::close_to((*f)(values), 0.0, 1.0));
}

BOOST_AUTO_TEST_SUITE_END()
