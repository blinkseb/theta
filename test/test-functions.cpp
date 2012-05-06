#include "interface/phys.hpp"
#include "test/utils.hpp"

#include <iostream>
#include <boost/test/unit_test.hpp>

using namespace theta;
using namespace std;

BOOST_AUTO_TEST_SUITE(functions)

BOOST_AUTO_TEST_CASE(add){
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    ParId p1 = vm->create_par_id("p1");
    ParId p2 = vm->create_par_id("p2");
    ParId p3 = vm->create_par_id("p3");
    load_core_plugins();
    ConfigCreator cc("f = {type = \"add\"; addends = (\"p1\", \"p2\", 1.7, 3.2, \"@f2\"); };\n"
                     "f2 = {type = \"add\"; addends = (\"p3\"); };"
            , vm);
    const theta::Configuration & cfg = cc.get();
    std::auto_ptr<Function> f_ptr;
    try{
        f_ptr = PluginManager<Function>::build(Configuration(cfg, cfg.setting["f"]));
    }
    catch(ConfigurationException & ex){
        cout << ex.message << endl;
        BOOST_CHECK(false);
        return;
    }
    const Function & f = *f_ptr;
    ParValues values;
    values.set(p1, 0.0);
    values.set(p2, 0.0);
    values.set(p3, 0.0);
    BOOST_REQUIRE(values.contains_all(f.get_parameters()));
    double expected = 1.7 + 3.2;
    BOOST_CHECK(close_to_relative(f(values),expected));
    values.set(p1, 24.7);
    values.set(p2, -25.2);
    values.set(p3, 42.0);
    expected = (1.7 + 3.2) + 24.7 - 25.2 + 42.0;
    BOOST_CHECK(close_to_relative(f(values), expected));
}

BOOST_AUTO_TEST_SUITE_END()

