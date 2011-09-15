#include "interface/histogram-function.hpp"
#include "interface/plugin.hpp"
#include "interface/histogram.hpp"
#include "interface/variables.hpp"
#include "test/utils.hpp"
#include <vector>

#include <boost/test/unit_test.hpp>

using namespace theta;
using namespace theta::plugin;
using namespace std;

BOOST_AUTO_TEST_SUITE(histogram_function_tests)

BOOST_AUTO_TEST_CASE(histofunction){
    VarIdManager vm;
    size_t nbins = 100;
    Histogram h(nbins, -1, 1);
    for(size_t i = 1; i<=nbins; i++){
        double x = 2.0 / 100 * i - 1;
        h.set(i, exp(-0.5*x*x));
    }
    std::auto_ptr<HistogramFunction> hp(new ConstantHistogramFunction(h));
    Histogram hh = (*hp)(ParValues());
    for(size_t i = 1; i<=nbins; i++){
        BOOST_REQUIRE(h.get(i) == hh.get(i));
    }
}

BOOST_AUTO_TEST_CASE(root_histogram_range){
    bool loaded =  load_root_plugins();
    if(!loaded){
       cout << "In test root_histogram_range: root plugin not loaded, not executing test" << endl;
       return;
    }
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    ConfigCreator cc(
            "root-histo1 = {type = \"root_histogram\"; filename=\"testhistos.root\"; histoname = \"histo1d\";};\n" // no range
            "root-histo2 = {type = \"root_histogram\"; filename=\"testhistos.root\"; histoname = \"histo1d\"; range = (-4.0, 20.0); };\n" // full range
            "root-histo3 = {type = \"root_histogram\"; filename=\"testhistos.root\"; histoname = \"histo1d\"; range = (-2.0, 2.0); };\n" // sub range
            "root-histo4 = {type = \"root_histogram\"; filename=\"testhistos.root\"; histoname = \"histo1d\"; range = (-4.1, 20.1); };\n" // range excluding underflow / overflow
            "root-histo5 = {type = \"root_histogram\"; filename=\"testhistos.root\"; histoname = \"histo1d\"; range = (-3.9, 19.9); };\n" // invalid range
            , vm);
    const theta::plugin::Configuration & cfg = cc.get();
    BOOST_CHECKPOINT("building hf");
    std::auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::instance().build(Configuration(cfg, cfg.setting["root-histo1"]));
    ParValues pv;
    Histogram h1 = (*hf)(pv);
    BOOST_REQUIRE(h1.get_nbins()==24);
    hf = PluginManager<HistogramFunction>::instance().build(Configuration(cfg, cfg.setting["root-histo2"]));
    Histogram h2 = (*hf)(pv);
    BOOST_REQUIRE(h2.get_nbins()==24);
    for(int i=0; i<=25; ++i){
        BOOST_ASSERT(h1.get(i) == h2.get(i));
    }
    hf = PluginManager<HistogramFunction>::instance().build(Configuration(cfg, cfg.setting["root-histo3"]));
    Histogram h3 = (*hf)(pv);
    BOOST_REQUIRE(h3.get_nbins()==4);
    for(int i=1; i<=4; ++i){
       BOOST_ASSERT(h3.get(i) == i+3.0);
    }
    hf = PluginManager<HistogramFunction>::instance().build(Configuration(cfg, cfg.setting["root-histo4"]));
    Histogram h4 = (*hf)(pv);
    BOOST_REQUIRE(h4.get_nbins()==26);
    for(int i=1; i<=26; ++i){
       BOOST_ASSERT(h4.get(i) == i);
    }
    
    bool except = false;
    try{
        hf = PluginManager<HistogramFunction>::instance().build(Configuration(cfg, cfg.setting["root-histo5"]));
    }
    catch(ConfigurationException & ex){
       except = true;
    }
    BOOST_ASSERT(except);
}


BOOST_AUTO_TEST_CASE(root_histogram){
    bool loaded =  load_root_plugins();
    if(!loaded){
       cout << "In test root_histogram: root plugin not loaded, not executing test" << endl;
       return;
    }
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    ConfigCreator cc(
            "root-histo1 = {type = \"root_histogram\"; filename=\"testhistos.root\"; histoname = \"histo1d\";};\n"
            "root-histo2 = {type = \"root_histogram\"; filename=\"testhistos.root\"; histoname = \"histo2d\";};\n"
            "root-histo3 = {type = \"root_histogram\"; filename=\"testhistos.root\"; histoname = \"histo3d\";};\n"
            , vm);
    const theta::plugin::Configuration & cfg = cc.get();
    BOOST_CHECKPOINT("building hf");
    std::auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::instance().build(Configuration(cfg, cfg.setting["root-histo1"]));
    ParValues pv;
    Histogram h = (*hf)(pv);
    BOOST_REQUIRE(h.get_nbins()==24);
    for(size_t i=1; i<=h.get_nbins(); ++i){
        BOOST_ASSERT(h.get(i)==i+1);
    }
    //overflow and underflow are not included:
    BOOST_ASSERT(h.get(0)==0);
    BOOST_ASSERT(h.get(h.get_nbins()+1)==0);
    //2D histogram:
    hf = PluginManager<HistogramFunction>::instance().build(Configuration(cfg, cfg.setting["root-histo2"]));
    h = (*hf)(pv);
    BOOST_REQUIRE(h.get_nbins()==10*11);
    //calculate the expected integral (excluding overflow / underflow!!);
    // if this matches, we are satisfied and believe the rest is Ok as well:
    double expected_integral = 0.0;
    for(int i=1; i<=10; ++i){
       for(int j=1; j<=11; ++j){
          expected_integral += (i + 0.78) * (j + 3.02);
       }
    }
    BOOST_ASSERT(utils::close_to_relative(h.get_sum_of_bincontents(),expected_integral));
    //3D histogram, same as 2D:
    hf = PluginManager<HistogramFunction>::instance().build(Configuration(cfg, cfg.setting["root-histo3"]));
    h = (*hf)(pv);
    BOOST_REQUIRE(h.get_nbins()==10*11*12);
    //calculate the expected integral (excluding overflow / underflow!!);
    // if this matches, we are satisfied and believe the rest is Ok as well:
    expected_integral = 0.0;
    for(int i=1; i<=10; ++i){
       for(int j=1; j<=11; ++j){
          for(int k=1; k<=12; ++k){
             expected_integral += (i + 0.12) * (j + 1.34) * (k + 5.67);
          }
       }
    }
    BOOST_ASSERT(utils::close_to_relative(h.get_sum_of_bincontents(),expected_integral));
}


BOOST_AUTO_TEST_CASE(cubiclinear_histomorph){
    load_core_plugins();
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    
    const size_t nbins = 1;
    ParId delta = vm->createParId("delta");
    ObsId obs = vm->createObsId("obs", nbins, -1, 1);
    BOOST_CHECKPOINT("parsing config");
    ConfigCreator cc("flat-histo0 = {type = \"fixed_poly\"; observable=\"obs\"; coefficients = [1.0]; normalize_to = 1.0;};\n"
            "flat-histo1 = {type = \"fixed_poly\"; observable=\"obs\"; coefficients = [1.0]; normalize_to = 1.12;};\n"
            "flat-histo-1 = {type = \"fixed_poly\"; observable=\"obs\"; coefficients = [1.0]; normalize_to = 0.83;};\n"
            "histo = { type = \"cubiclinear_histomorph\"; parameters = (\"delta\"); nominal-histogram = \"@flat-histo0\";\n"
            "      delta-plus-histogram = \"@flat-histo1\"; delta-minus-histogram = \"@flat-histo-1\";\n"
            "};\n"
            , vm);
            
    const theta::plugin::Configuration & cfg = cc.get();
    BOOST_CHECKPOINT("building hf");
    std::auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::instance().build(Configuration(cfg, cfg.setting["histo"]));
    BOOST_CHECKPOINT("hf built");
    ParValues pv;
    //check central and +- 1 sigma values:
    pv.set(delta, 0.0);
    Histogram h = (*hf)(pv);
    BOOST_ASSERT(utils::close_to_relative(h.get(1), 1.0));
    pv.set(delta, 1.0);
    h = (*hf)(pv);
    BOOST_ASSERT(utils::close_to_relative(h.get(1), 1.12));
    pv.set(delta, -1.0);
    h = (*hf)(pv);
    BOOST_ASSERT(utils::close_to_relative(h.get(1), 0.83));
    //+- 2 sigma values, should be interpolated linearly:
    pv.set(delta, 2.0);
    h = (*hf)(pv);
    BOOST_ASSERT(utils::close_to_relative(h.get(1), 1.24));
    pv.set(delta, -2.0);
    h = (*hf)(pv);
    BOOST_ASSERT(utils::close_to_relative(h.get(1), 0.66));
    //cutoff at zero:
    pv.set(delta, -10.0);
    h = (*hf)(pv);
    BOOST_ASSERT(h.get(1) == 0.0);
    //derivative at zero should be smooth:
    pv.set(delta, 1e-8);
    h = (*hf)(pv);
    double eps = h.get(1);
    pv.set(delta, -1e-8);
    h = (*hf)(pv);
    double eps_minus = h.get(1);
    BOOST_ASSERT(utils::close_to(eps - 1, 1 - eps_minus, 1000.0));
    //derivative at zero should be (0.12 + 0.17) / 2.
    double der = (eps - eps_minus) / (2e-8);
    BOOST_ASSERT(fabs(der - (0.12 + 0.17)/2) < 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
