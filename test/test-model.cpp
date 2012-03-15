#include "interface/phys.hpp"
#include "interface/histogram.hpp"
#include "interface/random.hpp"
#include "interface/plugin.hpp"
#include "interface/histogram-function.hpp"
#include "interface/model.hpp"

#include "test/utils.hpp"

#include "libconfig/libconfig.h++"

#include <iostream>

#include <boost/test/unit_test.hpp>

using namespace theta;
using namespace std;

namespace{
    Histogram1DWithUncertainties apply(const HistogramFunction & hf, const ParValues & values){
        Histogram1DWithUncertainties result;
        hf.apply_functor(copy_to<Histogram1DWithUncertainties>(result), values);
        return result;
    }
}


BOOST_AUTO_TEST_SUITE(model_tests)

BOOST_AUTO_TEST_CASE(model0){
    BOOST_CHECKPOINT("model0 entry");
    load_core_plugins();
    utils::fill_theta_dir(0);
    
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    ParIds pars;
    ObsIds obs;
    const size_t nbins = 100;
    ParId beta1 = vm->create_par_id("beta1");
    ParId beta2 = vm->create_par_id("beta2");
    ObsId obs0 = vm->create_obs_id("obs0", nbins, -1, 1);
    pars.insert(beta1);
    pars.insert(beta2);
    obs.insert(obs0);
    
    //boost::ptr_vector<Function> coeffs;
    BOOST_CHECKPOINT("parsing config");
    
    ConfigCreator cc("flat-histo = {type = \"fixed_poly\"; observable=\"obs0\"; coefficients = [1.0]; normalize_to = 1.0;};\n"
            "gauss-histo = {type = \"fixed_gauss\"; observable=\"obs0\"; width = 0.5; mean = 0.5; normalize_to = 1.0;};\n"
            "c1 = {type = \"multiply\"; factors=(\"beta1\");};\n"
            "c2 = {type = \"multiply\"; factors=(\"beta2\");};\n"
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
    
    const theta::Configuration & cfg = cc.get();
    BOOST_CHECKPOINT("building model");
    std::auto_ptr<Model> m;
    try{
        m = PluginManager<Model>::build(Configuration(cfg, cfg.setting["m"]));
    }
    catch(Exception & ex){
        std::cerr << ex.message << endl;
    }
    catch(libconfig::SettingNotFoundException & ex){
        std::cerr << ex.getPath() << " not found" << endl;
    }
    
    BOOST_REQUIRE(m.get()!=0);
    
    BOOST_CHECKPOINT("building signal histo");
    std::auto_ptr<HistogramFunction> f_signal_histo = PluginManager<HistogramFunction>::build(Configuration(cfg, cfg.setting["gauss-histo"]));
    BOOST_CHECKPOINT("building bkg histo");
    std::auto_ptr<HistogramFunction> f_bkg_histo = PluginManager<HistogramFunction>::build(Configuration(cfg, cfg.setting["flat-histo"]));
    
    ParValues values;
    Histogram1DWithUncertainties signal = apply(*f_signal_histo,values);
    Histogram1DWithUncertainties background = apply(*f_bkg_histo, values);
    
    values.set(beta1, 1.0);
    values.set(beta2, 0.0);
    Histogram1DWithUncertainties s;
    DataWithUncertainties pred;
    BOOST_CHECKPOINT("");
    m->get_prediction(pred, values);
    BOOST_CHECKPOINT("");
    s = pred[obs0];
    BOOST_CHECKPOINT("");
    //s should be signal only:
    for(size_t i = 0; i<nbins; i++){
        BOOST_REQUIRE(signal.get_value(i)==s.get_value(i));
    }
    //background only:
    values.set(beta1, 0.0);
    values.set(beta2, 1.0);
    m->get_prediction(pred, values);
    s = pred[obs0];
    for(size_t i = 0; i<nbins; i++){
        BOOST_REQUIRE(background.get_value(i)==s.get_value(i));
    }
    //zero prediction:
    values.set(beta1, 0.0);
    values.set(beta2, 0.0);
    m->get_prediction(pred, values);
    s = pred[obs0];
    for(size_t i = 0; i<nbins; i++){
        BOOST_REQUIRE(0.0==s.get_value(i));
    }

    //The likelihood, take double background. Use average as data:
    values.set(beta1, 1.0);
    values.set(beta2, 2.0);
    m->get_prediction(pred, values);
    s = pred[obs0];
    Data data;
    BOOST_CHECKPOINT("check");
    data[obs0] = s.get_values_histogram();
    BOOST_CHECKPOINT("check2");
    std::auto_ptr<NLLikelihood> nll = m->get_nllikelihood(data);
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

BOOST_AUTO_TEST_CASE(model_unc){
    load_core_plugins();
    utils::fill_theta_dir(0);
    
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    ParIds pars;
    ObsIds obs;
    const size_t nbins = 4;
    ParId beta1 = vm->create_par_id("beta1");
    ParId beta2 = vm->create_par_id("beta2");
    ParId beta3 = vm->create_par_id("beta3");
    ObsId obs0 = vm->create_obs_id("obs0", nbins, -1, 1);
    pars.insert(beta1);
    pars.insert(beta2);
    obs.insert(obs0);
    
    //boost::ptr_vector<Function> coeffs;
    BOOST_CHECKPOINT("parsing config");
    
    ConfigCreator cc("flat-histo = {type = \"direct_data_histo\"; nbins = 4; range = (-1.0, 1.0); data = (100.0, 100.0, 50.0, 50.0); uncertainties = (10.0, 10.0, 5.0, 5.0); };\n"
            "c1 = {type = \"multiply\"; factors=(\"beta1\");};\n"
            "c2 = {type = \"multiply\"; factors=(\"beta2\");};\n"
            "c3 = {type = \"multiply\"; factors=(\"beta3\");};\n"
            "dist-flat = {\n"
            "       type = \"flat_distribution\";\n"
            "       beta1 = { range = (\"-inf\", \"inf\"); fix-sample-value = 1.0; };\n"
            "       beta2 = { range = (\"-inf\", \"inf\"); fix-sample-value = 1.0; };\n"
            "       beta3 = { range = (\"-inf\", \"inf\"); fix-sample-value = 1.0; };\n"
            " };\n"
            "m = {\n"
            "  obs0 = {\n"
            "       signal = {\n"
            "          coefficient-function = \"@c1\";\n"
            "          histogram = \"@flat-histo\";\n"
            "       };\n"
            "       background = {\n"
            "           coefficient-function = \"@c2\";\n"
            "           histogram = \"@flat-histo\";\n"
            "       };\n"
            "       background2 = {\n"
            "           coefficient-function = \"@c3\";\n"
            "           histogram = \"@flat-histo\";\n"
            "       };\n"
            "   };\n"
            "  parameter-distribution = \"@dist-flat\";\n"
            "  bb_uncertainties = true;\n"
            "};\n"
            , vm);
    
    BOOST_CHECKPOINT("config parsed");
    
    const theta::Configuration & cfg = cc.get();
    BOOST_CHECKPOINT("building model");
    std::auto_ptr<Model> m;
    try{
        m = PluginManager<Model>::build(Configuration(cfg, cfg.setting["m"]));
    }
    catch(Exception & ex){
        std::cerr << ex.message << endl;
    }
    catch(libconfig::SettingNotFoundException & ex){
        std::cerr << ex.getPath() << " not found" << endl;
    }
    
    ParValues values;
    
    values.set(beta1, 1.0);
    values.set(beta2, 1.0);
    values.set(beta3, 1.0);
    Histogram1DWithUncertainties s;
    DataWithUncertainties pred;
    BOOST_CHECKPOINT("");
    m->get_prediction(pred, values);
    BOOST_CHECKPOINT("");
    s = pred[obs0];
    BOOST_CHECKPOINT("");
    //s should be signal only:
    BOOST_CHECK(close_to_relative(s.get(0), 300.0));
    BOOST_CHECK(close_to_relative(s.get(1), 300.0));
    BOOST_CHECK(close_to_relative(s.get(2), 150.0));
    BOOST_CHECK(close_to_relative(s.get(3), 150.0));
    
    BOOST_CHECK(close_to_relative(s.get_uncertainty(0), sqrt(3) * 10.0));
    BOOST_CHECK(close_to_relative(s.get_uncertainty(1), sqrt(3) * 10.0));
    BOOST_CHECK(close_to_relative(s.get_uncertainty(2), sqrt(3) * 5.0));
    BOOST_CHECK(close_to_relative(s.get_uncertainty(3), sqrt(3) * 5.0));
    
    values.set(beta1, 2.0);
    m->get_prediction(pred, values);
    s = pred[obs0];
    BOOST_CHECK(close_to_relative(s.get(0), 400.0));
    BOOST_CHECK(close_to_relative(s.get(1), 400.0));
    BOOST_CHECK(close_to_relative(s.get(2), 200.0));
    BOOST_CHECK(close_to_relative(s.get(3), 200.0));
    
    BOOST_CHECK(close_to_relative(s.get_uncertainty(0), sqrt(20*20 + 2 * 10*10)));
    BOOST_CHECK(close_to_relative(s.get_uncertainty(1), sqrt(20*20 + 2 * 10*10)));
    BOOST_CHECK(close_to_relative(s.get_uncertainty(2), sqrt(10*10 + 2 * 5*5)));
    BOOST_CHECK(close_to_relative(s.get_uncertainty(3), sqrt(10*10 + 2 * 5*5)));
}

BOOST_AUTO_TEST_SUITE_END()
