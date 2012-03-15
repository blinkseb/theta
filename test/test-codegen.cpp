#include "interface/phys.hpp"
#include "interface/model.hpp"
#include "interface/random.hpp"
#include "interface/histogram-function.hpp"
#include "test/utils.hpp"

#include "libconfig/libconfig.h++"

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <fstream>

using namespace std;
using namespace theta;

namespace{
    
bool operator==(const Histogram1DWithUncertainties& h0, const Histogram1DWithUncertainties & h1){
    if(h0.get_nbins()!=h1.get_nbins()) return false;
    if(h0.get_xmin()!=h1.get_xmin()) return false;
    if(h0.get_xmax()!=h1.get_xmax()) return false;
    for(size_t i=0; i<h0.get_nbins(); ++i){
        if(h0.get_value(i)!=h1.get_value(i)) return false;
        if(h0.get_uncertainty(i)!=h1.get_uncertainty(i)) return false;
    }
    return true;
}

bool operator==(const Histogram1D& h0, const Histogram1D & h1){
    if(h0.get_nbins()!=h1.get_nbins()) return false;
    if(h0.get_xmin()!=h1.get_xmin()) return false;
    if(h0.get_xmax()!=h1.get_xmax()) return false;
    for(size_t i=0; i<h0.get_nbins(); ++i){
        if(h0.get(i)!=h1.get(i)) return false;
    }
    return true;
}


}

BOOST_AUTO_TEST_SUITE(codegen)

BOOST_AUTO_TEST_CASE(llvm_function_wrapper){
    if(!load_llvm_plugins()){
        std::cout << "skipping llvm test" << endl;
        return;
    }
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    utils::fill_theta_dir(0);
    PropertyMap pm;
    pm.set("default", vm);
    vm->create_par_id("p0");
    ParId p1 = vm->create_par_id("p1");
    vm->create_par_id("p2");
    ConfigCreator cc(
            "f0 = {type = \"exp_function\"; parameters = (\"p1\"); lambdas_plus = (0.2); lambdas_minus = (0.1); };\n"
            "llvm_f0 = {type = \"llvm_enable_function\"; function = \"@f0\";};\n"
            , vm);
    const theta::Configuration & cfg = cc.get();
    BOOST_CHECKPOINT("building function");
    std::auto_ptr<Function> f0 = PluginManager<Function>::build(Configuration(cfg, cfg.setting["f0"]));
    std::auto_ptr<Function> llvm_f0;
    try{
         llvm_f0 = PluginManager<Function>::build(Configuration(cfg, cfg.setting["llvm_f0"]));
    }
    catch(Exception & ex){
        std::cout << ex.message << endl;
        throw;
    }
    ParValues vals;
    vals.set(p1, 1.7);
    BOOST_CHECK((*f0)(vals) == (*llvm_f0)(vals));
    vals.set(p1, -1.8);
    BOOST_CHECK((*f0)(vals) == (*llvm_f0)(vals));
}


BOOST_AUTO_TEST_CASE(llvm_multiply){
    if(!load_llvm_plugins()){
        return;
    }
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    utils::fill_theta_dir(0);
    PropertyMap pm;
    pm.set("default", vm);
    vm->create_par_id("p0");
    ParId p1 = vm->create_par_id("p1");
    ParId p2 = vm->create_par_id("p2");
    ConfigCreator cc(
            "f0 = {type = \"exp_function\"; parameters = (\"p1\"); lambdas_plus = (0.2); lambdas_minus = (0.1); };\n"
            "mul = {type = \"multiply\"; factors = (\"@f0\", \"p2\", 1.8); };\n"
            "llvm_mul = {type = \"llvm_multiply\"; factors = \"@mul.factors\";};\n"
            , vm);
    const theta::Configuration & cfg = cc.get();
    BOOST_CHECKPOINT("building function");
    std::auto_ptr<Function> mul = PluginManager<Function>::build(Configuration(cfg, cfg.setting["mul"]));
    std::auto_ptr<Function> llvm_mul = PluginManager<Function>::build(Configuration(cfg, cfg.setting["llvm_mul"]));
    ParValues vals;
    std::auto_ptr<RandomSource> rnd_src(new RandomSourceTaus());
    Random rnd(rnd_src);//no setSeed to keep it the same every time ...
    
    for(size_t i=0; i<30; ++i){
        vals.set(p1, rnd.uniform() * 10 - 5);
        vals.set(p2, rnd.uniform() * 10 - 5);
        BOOST_CHECK((*mul)(vals) == (*llvm_mul)(vals));
    }
}


BOOST_AUTO_TEST_CASE(llvm_exp_function){
    if(!load_llvm_plugins()){
        return;
    }
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    utils::fill_theta_dir(0);
    PropertyMap pm;
    pm.set("default", vm);
    vm->create_par_id("p0");
    ParId p1 = vm->create_par_id("p1");
    ParId p2 = vm->create_par_id("p2");
    ConfigCreator cc(
            "f = {type = \"exp_function\"; parameters = (\"p1\", \"p2\"); lambdas_plus = (0.1, 0.2); lambdas_minus = (0.1, 0.18); };\n"
            "llvm_f = {type = \"llvm_exp_function\"; parameters = \"@f.parameters\"; lambdas_plus = \"@f.lambdas_plus\"; lambdas_minus = \"@f.lambdas_minus\"; };\n"
            , vm);
    const theta::Configuration & cfg = cc.get();
    BOOST_CHECKPOINT("building function mul");
    std::auto_ptr<Function> mul = PluginManager<Function>::build(Configuration(cfg, cfg.setting["f"]));
    BOOST_CHECKPOINT("building function llvm_mul");
    std::auto_ptr<Function> llvm_mul = PluginManager<Function>::build(Configuration(cfg, cfg.setting["llvm_f"]));
    ParValues vals;
    
    std::auto_ptr<RandomSource> rnd_src(new RandomSourceTaus());
    Random rnd(rnd_src);//no setSeed to keep it the same every time ...
    
    for(size_t i=0; i<30; ++i){
        vals.set(p1, rnd.uniform() * 10 - 5);
        vals.set(p2, rnd.uniform() * 10 - 5);
        BOOST_CHECK((*mul)(vals) == (*llvm_mul)(vals));
    }
}


BOOST_AUTO_TEST_CASE(constant_histo){
    if(!load_llvm_plugins()){
        return;
    }
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    utils::fill_theta_dir(0);
    PropertyMap pm;
    pm.set("default", vm);
    ConfigCreator cc(
            "histo = {type = \"direct_data_histo\"; range = (0.0, 3.0); nbins = 3; data = (17.0, 23.0, 47.7); };\n"
            "llvm_h = {type = \"llvm_enable_histogram_function\"; histogram_function = \"@histo\"; debug = true; };\n"
            , vm);
    const theta::Configuration & cfg = cc.get();
    BOOST_CHECKPOINT("building histo");
    std::auto_ptr<HistogramFunction> histo = PluginManager<HistogramFunction>::build(Configuration(cfg, cfg.setting["histo"]));
    BOOST_CHECKPOINT("building llvm_h");
    std::auto_ptr<HistogramFunction> llvm_h = PluginManager<HistogramFunction>::build(Configuration(cfg, cfg.setting["llvm_h"]));
    BOOST_CHECKPOINT("evaluating");
    Histogram1DWithUncertainties h;
    histo->apply_functor(copy_to<Histogram1DWithUncertainties>(h), ParValues());
    Histogram1DWithUncertainties h_l;
    llvm_h->apply_functor(copy_to<Histogram1DWithUncertainties>(h_l), ParValues());
    BOOST_REQUIRE(h.get_nbins() == h_l.get_nbins());
    for(size_t i=0; i<h_l.get_nbins(); ++i){
        BOOST_REQUIRE(h.get_value(i) == h_l.get_value(i));
    }
}




BOOST_AUTO_TEST_CASE(model){
    if(!load_llvm_plugins()){
        return;
    }
    utils::fill_theta_dir(0);
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    PropertyMap pm;
    pm.set("default", vm);
    ParIds pars;
    ObsIds obs;
    const size_t nbins0 = 101;
    const size_t nbins1 = 50;
    ParId beta1 = vm->create_par_id("beta1");
    ParId beta2 = vm->create_par_id("beta2");
    ParId delta0 = vm->create_par_id("delta0");
    ParId delta1 = vm->create_par_id("delta1");
    ObsId obs0 = vm->create_obs_id("obs0", nbins0, -1, 1);
    ObsId obs1 = vm->create_obs_id("obs1", nbins1, -1, 1);
    pars.insert(beta1);
    pars.insert(beta2);
    obs.insert(obs0);
    
    //boost::ptr_vector<Function> coeffs;
    BOOST_CHECKPOINT("parsing config");
    
    ConfigCreator cc(
            "flat-histo0 = {type = \"fixed_poly\"; observable=\"obs0\"; coefficients = [1.0]; normalize_to = 1.0;};\n"
            "gauss-histo0 = {type = \"fixed_gauss\"; observable=\"obs0\"; width = 0.5; mean = 0.5; normalize_to = 1.0;};\n"
            "flat-histo1 = {type = \"fixed_poly\"; observable=\"obs1\"; coefficients = [1.0]; normalize_to = 1.7;};\n"
            "gauss-histo1 = {type = \"fixed_gauss\"; observable=\"obs1\"; width = 0.23; mean = 0.43; normalize_to = 1.12;};\n"
            "c-gauss = {type = \"multiply\"; factors=(1.1, \"beta1\", {type = \"exp_function\"; parameters = (\"delta0\", \"delta1\"); lambdas_plus = (0.1, 0.12); lambdas_minus = (0.1, 0.13);});};\n"
            "c-flat = {type = \"multiply\"; factors=(\"beta2\");};\n"
            "dist-flat = {\n"
            "       type = \"flat_distribution\";\n"
            "       beta1 = { range = (\"-inf\", \"inf\"); fix-sample-value = 1.0; };\n"
            "       beta2 = { range = (\"-inf\", \"inf\"); fix-sample-value = 1.0; };\n"
            "       delta0 = { range = (\"-inf\", \"inf\"); fix-sample-value = 0.0; };\n"
            "       delta1 = { range = (\"-inf\", \"inf\"); fix-sample-value = 0.0; };\n"
            " };\n"
            "m = {\n"
            "  obs0 = {\n"
            "       signal = {\n"
            "          coefficient-function = \"@c-gauss\";\n"
            "          histogram = \"@gauss-histo0\";\n"
            "       };\n"
            "       background = {\n"
            "           coefficient-function = \"@c-flat\";\n"
            "           histogram = \"@flat-histo0\";\n"
            "       };\n"
            "  };"
            "  obs1 = {\n"
            "       signal = {\n"
            "          coefficient-function = \"@c-gauss\";\n"
            "          histogram = \"@gauss-histo1\";\n"
            "       };\n"
            "       background = {\n"
            "           coefficient-function = \"@c-flat\";\n"
            "           histogram = \"@flat-histo1\";\n"
            "       };\n"
            "   };\n"
            "  parameter-distribution = \"@dist-flat\";\n"
            "};\n"
            "\n"
            "llvm_m = {\n"
            "   type = \"llvm_model\";\n"
            "   llvm_always = true;\n"
            "   obs0 = \"@m.obs0\";\n"
            "   obs1 = \"@m.obs1\";\n"
            "   parameter-distribution = \"@m.parameter-distribution\";\n"
            "};"
            , vm);
    
    BOOST_CHECKPOINT("config parsed");
    
    const theta::Configuration & cfg = cc.get();
    BOOST_CHECKPOINT("building model m");
    std::auto_ptr<Model> m;
    try{
       m = PluginManager<Model>::build(Configuration(cfg, cfg.setting["m"]));
    }
    catch(libconfig::SettingNotFoundException & ex){
        std::cout << ex.getPath() << endl;
        throw;
    }
    BOOST_CHECKPOINT("building model llvm_m");
    std::auto_ptr<Model> llvm_m = PluginManager<Model>::build(Configuration(cfg, cfg.setting["llvm_m"]));
    
    BOOST_REQUIRE(m->get_parameters() == llvm_m->get_parameters());
    BOOST_REQUIRE(m->get_observables() == llvm_m->get_observables());
        
    Data d;
    Data llvm_d;
    
    // calculate some predictions:
    ParValues vals;
    std::auto_ptr<RandomSource> rnd_src(new RandomSourceTaus());
    Random rnd(rnd_src);//no setSeed to keep it the same every time ...
    
    for(size_t i=0; i<50; ++i){
        vals.set(beta1, rnd.uniform() * 5);
        vals.set(beta2, rnd.uniform() * 5);
        vals.set(delta0, rnd.uniform() * 10 - 5);
        vals.set(delta1, rnd.uniform() * 10 - 5);
        m->get_prediction(d, vals);
        llvm_m->get_prediction(llvm_d, vals);
        BOOST_REQUIRE(d.get_observables().contains(obs0) && d.get_observables().contains(obs1));
        BOOST_REQUIRE(d.get_observables() == llvm_d.get_observables());
        BOOST_CHECK(d[obs0] == llvm_d[obs0]);
        BOOST_CHECK(d[obs1] == llvm_d[obs1]);
    }
    
    
    // calculate some nll values:
    vals.set(beta1, 1.17);
    vals.set(beta2, 0.83);
    vals.set(delta0, 0.12);
    vals.set(delta1, -0.18);    
    std::auto_ptr<NLLikelihood> nll = m->get_nllikelihood(d);
    std::auto_ptr<NLLikelihood> llvm_nll = llvm_m->get_nllikelihood(d);
    for(size_t i=0; i<50; ++i){
        vals.set(beta1, rnd.uniform() * 5);
        vals.set(beta2, rnd.uniform() * 5);
        vals.set(delta0, rnd.uniform() * 10 - 5);
        vals.set(delta1, rnd.uniform() * 10 - 5);
        double nll_value = (*nll)(vals);
        double llvm_nll_value = (*llvm_nll)(vals);
        //cout << i << " " << nll_value << " " << llvm_nll_value << " diff = " << fabs(nll_value - llvm_nll_value) << endl;
        // note: the values are not completely equal, as llvm_nll calculates the template likelihood over ALL bins in one
        // go and nll calculates it for each observable seperately. The differences can only be at most N_obs roundings, though
        BOOST_CHECK(close_to_relative(nll_value, llvm_nll_value));
    }
}


BOOST_AUTO_TEST_SUITE_END()

