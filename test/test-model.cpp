#include "interface/phys.hpp"
#include "interface/histogram.hpp"
#include "interface/random.hpp"
#include "interface/plugin.hpp"
#include "interface/histogram-function.hpp"
#include "interface/model.hpp"

#include "test/utils.hpp"

#include "libconfig/libconfig.h++"

#include <boost/test/unit_test.hpp>

using namespace theta;
using namespace theta::plugin;
using namespace std;


BOOST_AUTO_TEST_SUITE(model_tests)

BOOST_AUTO_TEST_CASE(model0){
    BOOST_CHECKPOINT("model0 entry");
    load_core_plugins();
    
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    ParIds pars;
    ObsIds obs;
    const size_t nbins = 100;
    ParId beta1 = vm->createParId("beta1");
    ParId beta2 = vm->createParId("beta2");
    ObsId obs0 = vm->createObsId("obs0", nbins, -1, 1);
    pars.insert(beta1);
    pars.insert(beta2);
    obs.insert(obs0);
    
    //boost::ptr_vector<Function> coeffs;
    BOOST_CHECKPOINT("parsing config");
    
    ConfigCreator cc("flat-histo = {type = \"fixed_poly\"; observable=\"obs0\"; coefficients = [1.0]; normalize_to = 1.0;};\n"
            "gauss-histo = {type = \"fixed_gauss\"; observable=\"obs0\"; width = 0.5; mean = 0.5; normalize_to = 1.0;};\n"
            "c1 = {type = \"mult\"; parameters=(\"beta1\");};\n"
            "c2 = {type = \"mult\"; parameters=(\"beta2\");};\n"
            "dist-flat = {\n"
            "       type = \"flat_distribution\";\n"
            "       beta1 = { range = (\"-inf\", \"inf\"); fix-sample-value = 1.0; }; \n"
            "       beta2 = { range = (\"-inf\", \"inf\"); fix-sample-value = 1.0; };\n"
            " };\n"
            "m = {\n"
            "  obs0 = {\n"
            "       signal = {\n"
            "          coefficient-function = \"@c1\";\n"
            "          histogram = \"@gauss-histo\";\n"
            "       };\n"
            "       background = {\n"
            "           coefficient-function = \"@c2\";\n"
            "           histogram = \"@flat-histo\";\n"
            "       };\n"
            "   };\n"
            "  parameter-distribution = \"@dist-flat\";\n"
            "};\n"
            , vm);
    
    BOOST_CHECKPOINT("config parsed");
    
    const theta::plugin::Configuration & cfg = cc.get();
    BOOST_CHECKPOINT("building model");
    std::auto_ptr<Model> m;
    try{
        m = PluginManager<Model>::instance().build(Configuration(cfg, cfg.setting["m"]));
    }
    catch(Exception & ex){
        cerr << ex.message << endl;
    }
    catch(libconfig::SettingNotFoundException & ex){
        cerr << ex.getPath() << " not found" << endl;
    }
    
    BOOST_REQUIRE(m.get()!=0);
    
    BOOST_CHECKPOINT("building signal histo");
    std::auto_ptr<HistogramFunction> f_signal_histo = PluginManager<HistogramFunction>::instance().build(Configuration(cfg, cfg.setting["gauss-histo"]));
    BOOST_CHECKPOINT("building bkg histo");
    std::auto_ptr<HistogramFunction> f_bkg_histo = PluginManager<HistogramFunction>::instance().build(Configuration(cfg, cfg.setting["flat-histo"]));
    
    ParValues values;
    Histogram signal = (*f_signal_histo)(values);
    Histogram background = (*f_bkg_histo)(values);
    
    values.set(beta1, 1.0);
    values.set(beta2, 0.0);
    Histogram s;
    Data pred;
    BOOST_CHECKPOINT("");
    m->get_prediction(pred, values);
    BOOST_CHECKPOINT("");
    s = pred[obs0];
    BOOST_CHECKPOINT("");
    //s should be signal only:
    for(size_t i = 1; i<=nbins; i++){
        BOOST_REQUIRE(signal.get(i)==s.get(i));
    }
    //background only:
    values.set(beta1, 0.0);
    values.set(beta2, 1.0);
    m->get_prediction(pred, values);
    s = pred[obs0];
    for(size_t i = 1; i<=nbins; i++){
        BOOST_REQUIRE(background.get(i)==s.get(i));
    }
    //zero prediction:
    values.set(beta1, 0.0);
    values.set(beta2, 0.0);
    m->get_prediction(pred, values);
    s = pred[obs0];
    for(size_t i = 1; i<=nbins; i++){
        BOOST_REQUIRE(0.0==s.get(i));
    }

    //The likelihood, take double background. Use average as data:
    values.set(beta1, 1.0);
    values.set(beta2, 2.0);
    m->get_prediction(pred, values);
    s = pred[obs0];
    Data data;
    BOOST_CHECKPOINT("check");
    data[obs0] = s;
    BOOST_CHECKPOINT("check2");
    std::auto_ptr<NLLikelihood> nll = m->getNLLikelihood(data);
    BOOST_CHECKPOINT("check3");
    double x[2];
    x[0] = 0.9;
    x[1] = 1.9;
    double nll09 = (*nll)(x);
    x[0] = 1.0;
    x[1] = 2.0;
    double nll10 = (*nll)(x);
    x[0] = 1.1;
    x[1] = 2.1;
    double nll11 = (*nll)(x);
    BOOST_CHECKPOINT("check4");
    BOOST_CHECK(nll10 < nll11);
    BOOST_CHECK(nll10 < nll09);
}



BOOST_AUTO_TEST_SUITE_END()
