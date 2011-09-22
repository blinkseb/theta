#include "interface/phys.hpp"
#include "interface/model.hpp"
#include "test/utils.hpp"

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <fstream>

using namespace std;
using namespace theta;

BOOST_AUTO_TEST_SUITE(codegen)

BOOST_AUTO_TEST_CASE(function){
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    utils::fill_theta_dir(0);
    PropertyMap pm;
    pm.set("default", vm);
    vm->createParId("p0");
    vm->createParId("p1");
    vm->createParId("p2");
    ConfigCreator cc(
            "f0 = {type = \"exp_function\"; parameters = (\"p1\"); lambdas_plus = (0.2); lambdas_minus = (0.1); };"
            , vm);
    const theta::Configuration & cfg = cc.get();
    BOOST_CHECKPOINT("building function");
    std::auto_ptr<Function> f0 = PluginManager<Function>::instance().build(Configuration(cfg, cfg.setting["f0"]));
    ofstream out("testgen.cpp");
    codegen::header(out, pm, f0->getParameters(), ObsIds());
    f0->codegen(out, "f0", pm);
    codegen::footer(out);
    void * h = compile_load_so("testgen.cpp");
    codegen::t_function_evaluate f0_generated = reinterpret_cast<codegen::t_function_evaluate>(dlsym(h, "f0_evaluate"));
    double par_values[1];
    par_values[0] = 1.7;
    double res = f0_generated(par_values);
    BOOST_CHECK(res == exp(0.2 * 1.7));
    par_values[0] = -0.3;
    res = f0_generated(par_values);
    BOOST_CHECK(res == exp( -0.1 * 0.3));
}

BOOST_AUTO_TEST_CASE(model){
    load_core_plugins();
    utils::fill_theta_dir(0);
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    PropertyMap pm;
    pm.set("default", vm);
    ParIds pars;
    ObsIds obs;
    const size_t nbins = 100;
    ParId beta1 = vm->createParId("beta1");
    ParId beta2 = vm->createParId("beta2");
    ParId delta0 = vm->createParId("delta0");
    ParId delta1 = vm->createParId("delta1");
    ObsId obs0 = vm->createObsId("obs0", nbins, -1, 1);
    pars.insert(beta1);
    pars.insert(beta2);
    obs.insert(obs0);
    
    //boost::ptr_vector<Function> coeffs;
    BOOST_CHECKPOINT("parsing config");
    
    ConfigCreator cc("flat-histo = {type = \"fixed_poly\"; observable=\"obs0\"; coefficients = [1.0]; normalize_to = 1.0;};\n"
            "gauss-histo = {type = \"fixed_gauss\"; observable=\"obs0\"; width = 0.5; mean = 0.5; normalize_to = 1.0;};\n"
            "c1 = {type = \"multiply\"; factors=(1.1, \"beta1\", {type = \"exp_function\"; parameters = (\"delta0\", \"delta1\"); lambdas_plus = (0.1, 0.12); lambdas_minus = (0.1, 0.13);});};\n"
            "c2 = {type = \"multiply\"; factors=(\"beta2\");};\n"
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
    m = PluginManager<Model>::instance().build(Configuration(cfg, cfg.setting["m"]));
    BOOST_CHECKPOINT("model built, generating code");
    ofstream out("testgen_model.cpp");
    codegen::header(out, pm, m->getParameters(), m->getObservables());
    m->codegen(out, "m", pm);
    codegen::footer(out);
    BOOST_CHECKPOINT("generated code; compiling and loading");
    void * h = compile_load_so("testgen_model.cpp");
    codegen::t_model_get_prediction m_pred_generated = reinterpret_cast<codegen::t_model_get_prediction>(dlsym(h, "m_get_prediction"));
    
    // get generated prediction:
    double par_values[4] = {0.8, 1.7, 0.1, -0.1};
    Histogram1D pred(100);
    m_pred_generated(par_values, pred.getData());
    
    // get original prediction:
    ParValues values;
    values.set(beta1, 0.8);
    values.set(beta2, 1.7);
    values.set(delta0, 0.1);
    values.set(delta1, -0.1);
    Data d;
    m->get_prediction(d, values);
    
    //check if they are equal:
    BOOST_REQUIRE(d[obs0].get_nbins() == nbins);
    for(size_t i=0; i<nbins; ++i){
        BOOST_ASSERT(d[obs0].get(i) == pred.get(i));
    }
    
    // use the prediction to calculate the likelihood functions: compiled and "normal":
    {
        std::auto_ptr<NLLikelihood> nll = m->getNLLikelihood(d);
        BOOST_CHECKPOINT("getting compiled nll");
        try{
            std::auto_ptr<NLLikelihood> nll_gen = dynamic_cast<default_model&>(*m).getCompiledNLLikelihood(d);
            BOOST_ASSERT((*nll)(values) == (*nll_gen)(values));
        }
        catch(Exception & ex){
           cout << ex.message << endl;
           BOOST_ASSERT(false);
        }
    }
    
    values.set(beta1, 1.3);
    values.set(delta0, 1.5);
    values.set(delta1, -1.5);
    m->get_prediction(d, values);
    {
        std::auto_ptr<NLLikelihood> nll = m->getNLLikelihood(d);
        BOOST_CHECKPOINT("getting compiled nll");
        std::auto_ptr<NLLikelihood> nll_gen = dynamic_cast<default_model&>(*m).getCompiledNLLikelihood(d);
        BOOST_ASSERT((*nll)(values) == (*nll_gen)(values));
        values.set(beta2, 0.888);
        BOOST_ASSERT((*nll)(values) == (*nll_gen)(values));
        values.set(delta1, 1.1);
        BOOST_ASSERT((*nll)(values) == (*nll_gen)(values));
    }
}


BOOST_AUTO_TEST_SUITE_END()

