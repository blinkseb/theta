#ifndef LLVM_MODEL_HPP
#define LLVM_MODEL_HPP

#include "interface/decls.hpp"
#include "interface/model.hpp"
#include "llvm/llvm_interface.hpp"

#include <vector>
#include <string>
#include <set>
#include <map>


typedef double (*t_model_get_prediction)(const double *, double *);

class llvm_model_nll;


// like default_model, with same configuration.
// One additional parameter: llvm_always = true;  (default is false) will use the
// compiled prediction function not only for likelihood construction but also for get_prediction.
class llvm_model: public theta::Model{
friend class llvm_model_nll;
private:
    typedef boost::ptr_map<theta::ObsId, boost::ptr_vector<theta::HistogramFunction> > histos_type;
    typedef boost::ptr_map<theta::ObsId, boost::ptr_vector<theta::Function> > coeffs_type;
    histos_type histos;
    coeffs_type coeffs;
    std::auto_ptr<theta::Distribution> parameter_distribution;
    bool llvm_always;
    
    mutable t_model_get_prediction model_get_prediction;
    mutable std::auto_ptr<llvm_module> module;
    
    void set_prediction(const theta::ObsId & obs_id, boost::ptr_vector<theta::Function> & coeffs, boost::ptr_vector<theta::HistogramFunction> & histos);
    // generate and compile llvm, fill mode_get_prediction function pointer.
    void generate_llvm() const;
    
 public:
    llvm_model(const theta::Configuration & cfg);
    //the pure virtual functions:
    virtual void get_prediction_randomized(theta::Random & rnd, theta::Data & result, const theta::ParValues & parameters) const;
    virtual void get_prediction(theta::Data & result, const theta::ParValues & parameters) const;
    virtual std::auto_ptr<theta::NLLikelihood> getNLLikelihood(const theta::Data & data) const;
    
    virtual const theta::Distribution & get_parameter_distribution() const {
       return *parameter_distribution;
    }
    virtual ~llvm_model();
    
};
 
class llvm_model_nll: public theta::NLLikelihood{
friend class llvm_model;
public:
    using theta::Function::operator();
    virtual double operator()(const theta::ParValues & values) const;
    
    virtual void set_additional_term(const boost::shared_ptr<theta::Function> & term);
    virtual void set_override_distribution(const boost::shared_ptr<theta::Distribution> & d);
    virtual const theta::Distribution & get_parameter_distribution() const{
        if(override_distribution) return *override_distribution;
        else return model.get_parameter_distribution();
    }        
private:
    const llvm_model & model;
    
    boost::shared_ptr<theta::Function> additional_term;
    boost::shared_ptr<theta::Distribution> override_distribution;

    theta::Histogram1D data_concatenated;
    t_model_get_prediction model_get_prediction;
    //cached predictions:
    mutable theta::Histogram1D pred_concatenated;
    mutable std::vector<double> parameter_values;
    
    llvm_model_nll(const llvm_model & m, const theta::Data & data, t_model_get_prediction model_get_prediction);
 };

#endif

