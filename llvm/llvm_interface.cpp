#include "llvm/llvm_interface.hpp"

#include "llvm/DerivedTypes.h"
#include "llvm/LLVMContext.h"
#include "llvm/Module.h"
#include "llvm/Support/IRBuilder.h"
#include "llvm/Support/StandardPasses.h"
#include "llvm/Target/TargetSelect.h"
#include "llvm/Target/TargetData.h"

#include "llvm/Analysis/Verifier.h"
#include "llvm/Attributes.h"

#include "llvm/PassManager.h"
#include "llvm/Transforms/Scalar.h"
#include "llvm/Transforms/IPO.h"
#include "llvm/ExecutionEngine/ExecutionEngine.h"
#include "llvm/ExecutionEngine/JIT.h"

#include "interface/plugin.hpp"
#include <dlfcn.h>


using namespace llvm;
using namespace theta;
using namespace std;

namespace {
    struct llvm_initializer{
        llvm_initializer(){
            InitializeNativeTarget();
        }
    };
    llvm_initializer init;

    std::map<theta::ParId, size_t> create_pid_to_index(const ParIds & pids){
        std::map<theta::ParId, size_t> result;
        size_t i=0;
        for(ParIds::const_iterator it=pids.begin(); it!=pids.end(); ++it){
            result[*it] = i++;
        }
        return result;
    }

}

// create the module function
// inline void add_with_coeff(double coeff, double * data, const double * rhs_data, int n);
// which calculates data += coeff * rhs_data where data and rhs_data are 16-byte aligned double-vectors of
// length n  (or n + (n%2)).
void emit_add_with_coeff(llvm::Module * mod){
    LLVMContext & context = mod->getContext();
    const Type * double_t = Type::getDoubleTy(context);
    VectorType* double2_t = VectorType::get(double_t, 2);
    const Type * void_t = Type::getVoidTy(context);
    const Type * i32_t = Type::getInt32Ty(context);
    std::vector<const Type*> arg_types(4);
    arg_types[0] = double_t;
    arg_types[1] = arg_types[2] = double_t->getPointerTo();
    arg_types[3] = i32_t;
    FunctionType * FT = FunctionType::get(void_t, arg_types, false);
    llvm::Function * F = llvm::Function::Create(FT, llvm::Function::ExternalLinkage, "add_with_coeff", mod);
    F->addFnAttr(Attribute::InlineHint);
    llvm::Function::arg_iterator iter = F->arg_begin();
    Argument * coeff = iter++;
    Argument * data_ = iter++;
    Argument * rhs_data_ = iter++;
    Argument * n_ = iter;
    data_->addAttr(Attribute::constructAlignmentFromInt(16));
    rhs_data_->addAttr(Attribute::constructAlignmentFromInt(16));
    
    IRBuilder<> Builder(context);
    BasicBlock * BB = BasicBlock::Create(context, "entry", F);
    Builder.SetInsertPoint(BB);
    //coeff2 = Vector<double, 2>(coeff, coeff):
    std::vector<Constant *> zeros(2);
    zeros[0] = zeros[1] = ConstantFP::get(double_t, 0.0);
    Value * coeff2tmp = Builder.CreateInsertElement(ConstantVector::get(zeros), coeff, ConstantInt::get(i32_t, 0));
    Value * coeff2 = Builder.CreateInsertElement(coeff2tmp, coeff, ConstantInt::get(i32_t, 1), "coeff2");
    // n = n_ + n_%2  (n_%2 = n_ % 1)
    Value * n_mod_2 = Builder.CreateAnd(n_, ConstantInt::get(i32_t, 1), "n_mod_2");
    Value * n = Builder.CreateAdd(n_, n_mod_2, "n_aligned");
    Value * n_half = Builder.CreateUDiv(n, ConstantInt::get(i32_t, 2), "n2");
    // create data and rhs_data to point to data_ / rhs_data_ but of type pointer to double2_t:
    Value * data = Builder.CreateBitCast(data_, double2_t->getPointerTo(), "data");
    Value * rhs_data = Builder.CreateBitCast(rhs_data_, double2_t->getPointerTo(), "rhs_data");
    // make a loop i = 0..n_half
    BasicBlock * BB_loop_header = BasicBlock::Create(context, "loop_header", F); // here compare i with n and jump to body or exit
    BasicBlock * BB_loop_body = BasicBlock::Create(context, "loop_body", F);
    BasicBlock * BB_loop_exit = BasicBlock::Create(context, "loop_exit", F);
    Builder.CreateBr(BB_loop_header);
    Builder.SetInsertPoint(BB_loop_header);
    PHINode * i_phi = Builder.CreatePHI(i32_t, "i");
    theta_assert(i_phi != 0);
    i_phi->addIncoming(ConstantInt::get(i32_t, 0), BB);
    // add the i_phi incoming value for next loop value later
    Value * i_eq_n = Builder.CreateICmpEQ(i_phi, n_half);
    Builder.CreateCondBr(i_eq_n, BB_loop_exit, BB_loop_body);
    Builder.SetInsertPoint(BB_loop_body);
    // get data[i] and rhs_data[i] as double2_t
    Value * p_data_i = GetElementPtrInst::Create(data, i_phi, "p_data_i", BB_loop_body);
    Value * p_rhs_data_i = GetElementPtrInst::Create(rhs_data, i_phi, "p_rhs_data_i", BB_loop_body);
    Value * rhs_data_i = Builder.CreateLoad(p_rhs_data_i, "rhs_data_i");
    Value * multmp = Builder.CreateFMul(rhs_data_i, coeff2);
    Value * data_i = Builder.CreateLoad(p_data_i, "data_i");
    Value * res = Builder.CreateFAdd(multmp, data_i, "res");
    Builder.CreateStore(res, p_data_i);
    Value * next_i = Builder.CreateAdd(i_phi, ConstantInt::get(i32_t, 1), "next_i");
    i_phi->addIncoming(next_i, BB_loop_body);
    Builder.CreateBr(BB_loop_header);
    Builder.SetInsertPoint(BB_loop_exit);
    Builder.CreateRetVoid();
    //llvm_verify(F, "add_with_coeff");
}


llvm_module::llvm_module(const ParIds & pids_): pids(pids_), pid_to_index(create_pid_to_index(pids)){
   LLVMContext & context = getGlobalContext();
   module =  new Module("llvm_module", context);
   std::string err;
   ee =  ExecutionEngine::createJIT(module, &err, 0, CodeGenOpt::Aggressive);
   if(!ee){
       throw Exception("llvm_module: could not create llvm::ExecutionEngine: " + err);
   }
   //define function prototype for the codegen_*_evaluate functions:
   //double codegen_f_evaluate(llvm_module* mod, void* fptr, const double * par_value);
   const Type * double_t = Type::getDoubleTy(context);
   const Type * void_t = Type::getVoidTy(context);
   const Type * p_char_t = Type::getInt8Ty(context)->getPointerTo();
   std::vector<const Type*> arg_types(3, p_char_t);
   arg_types[2] = double_t->getPointerTo();
   FunctionType * FT = FunctionType::get(double_t, arg_types, false);
   llvm::Function::Create(FT, llvm::Function::ExternalLinkage, "codegen_f_evaluate", module);
   //void codegen_hf_add_with_coeff(llvm_module* mod, void* hfptr, double coeff, const double * par_values, double * data)
   arg_types.resize(5);
   arg_types[2] = double_t;
   arg_types[3] = arg_types[4] = double_t->getPointerTo();
   FunctionType * FT2 = FunctionType::get(void_t, arg_types, false);
   llvm::Function::Create(FT2, llvm::Function::ExternalLinkage, "codegen_hf_add_with_coeff", module);
   
   // create an alias 'exp_function' which either redirects to amd_exp (if it exists) or to exp:
   arg_types.resize(1);
   arg_types[0] = double_t;
   FunctionType * FT_exp = FunctionType::get(double_t, arg_types, false);
   void * handle = dlopen(0, 0);
   bool amd_exp_available = handle && dlsym(handle, "amd_exp");
   llvm::Function * f_exp = llvm::Function::Create(FT_exp, llvm::Function::ExternalLinkage, amd_exp_available?"amd_exp":"exp", module);
   new GlobalAlias(FT_exp, GlobalValue::PrivateLinkage, "exp_function", f_exp, module);
   emit_add_with_coeff(module);
   //module->dump();
}

size_t llvm_module::get_index(const theta::ParId & pid) const{
    map<theta::ParId, size_t>::const_iterator it = pid_to_index.find(pid);
    if(it==pid_to_index.end()) throw invalid_argument("llvm_module::get_index");
    return it->second;
}

void* llvm_module::getFunctionPointer(llvm::Function * function){
    return ee->getPointerToFunction(function);
}

llvm_module::~llvm_module(){
    delete ee;
}

void llvm_module::optimize(){
    PassManager pm;    
    // this:
    //createStandardFunctionPasses(&pm, 3);
    // does not work because it need a FunctionPassManager, so add it by hand:
    pm.add(createCFGSimplificationPass());
    pm.add(createScalarReplAggregatesPass());
    pm.add(createInstructionCombiningPass());
    createStandardModulePasses(&pm, 3, /* OptimizeSize = */ false, /* UnitAtATime = */ true, /* UnrollLoops = */ true,
        /* SimplifyLibCalls= */ false, /* HaveExceptions = */ false, createFunctionInliningPass());
    pm.run(*module);
}

void set_private_inline(llvm::Function * f){
    f->addFnAttr(Attribute::AlwaysInline);
    f->setLinkage(GlobalValue::PrivateLinkage);
}

// get llvm function types for the functions created by HistogramFunction and Function code
// generators:
llvm::FunctionType * get_ft_hf_add_with_coeff(llvm_module & mod){
    LLVMContext & context = mod.module->getContext();
    const Type * double_t = Type::getDoubleTy(context);
    const Type * void_t = Type::getVoidTy(context);
    std::vector<const Type*> arg_types(3);
    arg_types[0] = double_t;
    arg_types[1] = arg_types[2] = double_t->getPointerTo();
    return FunctionType::get(void_t, arg_types, false);
}

llvm::FunctionType * get_ft_function_evaluate(llvm_module & mod){
    LLVMContext & context = mod.module->getContext();
    const Type * double_t = Type::getDoubleTy(context);
    std::vector<const Type*> arg_types(1);
    arg_types[0] = double_t->getPointerTo();
    return FunctionType::get(double_t, arg_types, false);
}

void llvm_verify(llvm::Function* f, const std::string & fname){
    if(verifyFunction(*f, llvm::ReturnStatusAction)){
        f->dump();
        throw invalid_argument("llvm_verify failed for function " + fname);
    }
}

GlobalVariable * add_global_const_ddata(Module * m, const std::string & name, const double * data, size_t n){
    size_t n_orig = n;
    if(n%2)++n;
    std::vector<Constant*> elements(n);
    for(size_t i=0; i<n_orig; ++i){
        elements[i] = ConstantFP::get(m->getContext(), APFloat(data[i]));
    }
    if(n_orig % 2){
        elements[n_orig] = ConstantFP::get(m->getContext(), APFloat(0.0));
    }
    ArrayType* typ = ArrayType::get(Type::getDoubleTy(m->getContext()), n);
    Constant * ini = ConstantArray::get(typ, elements);
    GlobalVariable * gvar = new GlobalVariable(*m, typ, true, GlobalValue::PrivateLinkage, ini, name);
    gvar->setAlignment(16);
    return gvar;
}

llvm::Function * llvm_generic_codegen(const theta::Function * f, llvm_module & mod, const std::string & prefix){
    FunctionType * FT = get_ft_function_evaluate(mod);
    llvm::Function * F = llvm::Function::Create(FT, llvm::Function::ExternalLinkage, prefix + "_evaluate", mod.module);
    LLVMContext & context = mod.module->getContext();
    IRBuilder<> Builder(context);
    llvm::Function * llvm_codegen_f_evaluate = mod.module->getFunction("codegen_f_evaluate");
    const Type * p_char_t = Type::getInt8Ty(context)->getPointerTo();
    const Type * i64_t = Type::getInt64Ty(context);
    BasicBlock * BB = BasicBlock::Create(context, "entry", F);
    Builder.SetInsertPoint(BB);
    Value * par_values = &*(F->arg_begin());
    Value * mod_ = Builder.CreateBitCast(ConstantInt::get(i64_t, reinterpret_cast<unsigned long>(&mod)), p_char_t);
    Value * fptr = Builder.CreateBitCast(ConstantInt::get(i64_t, reinterpret_cast<unsigned long>(f)), p_char_t);
    Value * ret = Builder.CreateCall3(llvm_codegen_f_evaluate, mod_, fptr, par_values);
    Builder.CreateRet(ret);
    return F;
}

llvm::Function * llvm_generic_codegen(const theta::HistogramFunction * hf, llvm_module & mod, const std::string & prefix){
    FunctionType * FT = get_ft_hf_add_with_coeff(mod);
    llvm::Function * F = llvm::Function::Create(FT, llvm::Function::ExternalLinkage, prefix + "_add_with_coeff", mod.module);
    LLVMContext & context = mod.module->getContext();
    IRBuilder<> Builder(context);
    llvm::Function * llvm_codegen_hf_add_with_coeff = mod.module->getFunction("codegen_hf_add_with_coeff");
    const Type * p_char_t = Type::getInt8Ty(context)->getPointerTo();
    const Type * i64_t = Type::getInt64Ty(context);
    BasicBlock * BB = BasicBlock::Create(context, "entry", F);
    Builder.SetInsertPoint(BB);
    llvm::Function::arg_iterator iter = F->arg_begin();
    Value * coeff = iter++;
    Value * par_values = iter++;
    Value * data = iter;
    Value * mod_ = Builder.CreateBitCast(ConstantInt::get(i64_t, reinterpret_cast<unsigned long>(&mod)), p_char_t);
    Value * hfptr = Builder.CreateBitCast(ConstantInt::get(i64_t, reinterpret_cast<unsigned long>(hf)), p_char_t);
    Value *Args[] = { mod_, hfptr, coeff, par_values, data };
    Builder.Insert(CallInst::Create(llvm_codegen_hf_add_with_coeff, Args, Args+5), "");
    Builder.CreateRetVoid();
    return F;
}

llvm::Function * create_llvm_function(const theta::Function * f, llvm_module & mod, const std::string & prefix){
    llvm::Function * result;
    if(const llvm_enabled<theta::Function>* en = dynamic_cast<const llvm_enabled<theta::Function>*>(f)){
        result = en->llvm_codegen(mod, prefix);
    }
    else{
        result = llvm_generic_codegen(f, mod, prefix);
    }
    llvm_verify(result, prefix);
    return result;
}

llvm::Function * create_llvm_histogram_function(const theta::HistogramFunction * hf, llvm_module & mod, const std::string & prefix){
    llvm::Function * result;
    if(const llvm_enabled<theta::HistogramFunction>* en = dynamic_cast<const llvm_enabled<theta::HistogramFunction>*>(hf)){
        result = en->llvm_codegen(mod, prefix);
    }
    else if(hf->get_parameters().size()==0){
        Histogram1D h0;
        hf->apply_functor(copy_to<Histogram1D>(h0), ParValues());
        GlobalVariable * gv_data = add_global_const_ddata(mod.module, prefix + "_data", h0.get_data(), h0.size());
        llvm::FunctionType * FT = get_ft_hf_add_with_coeff(mod);
        result = llvm::Function::Create(FT, llvm::Function::ExternalLinkage, prefix + "_add_with_coeff", mod.module);
        llvm::Function * add_with_coeff = mod.module->getFunction("add_with_coeff");
        LLVMContext & context = mod.module->getContext();
        IRBuilder<> Builder(context);
        llvm::Function::arg_iterator iter = result->arg_begin();
        Value * coeff = iter++;
        /*Value * par_values =*/ iter++; // not used ...
        Value * data = iter;
        const Type * double_t = Type::getDoubleTy(context);
        const Type * i32_t = Type::getInt32Ty(context);
        Value * gv_data_conv = Builder.CreateBitCast(gv_data, double_t->getPointerTo());
        BasicBlock * BB = BasicBlock::Create(context, "entry", result);
        Builder.SetInsertPoint(BB);
        Builder.CreateCall4(add_with_coeff, coeff, data, gv_data_conv, ConstantInt::get(i32_t, h0.size()));
        Builder.CreateRetVoid();
        //mod.module->dump();
    }else{
        result = llvm_generic_codegen(hf, mod, prefix);
    }
    llvm_verify(result, prefix);
    return result;
}



namespace{
    
    class add_to_vdouble: public functor<Histogram1D>{
    private:
        mutable double * data;
        double coeff;
    public:
        add_to_vdouble(double * data_, double coeff_): data(data_), coeff(coeff_){}
        virtual void operator()(const Histogram1D & h) const{
            utils::add_fast_with_coeff(data, h.get_data(), coeff, h.size());
        }
    };
    
}

// the callback functions llvm will call for Functions and HistogramFunctions which are not derived from llvm_enabled<T>:
extern "C" {
double codegen_f_evaluate(llvm_module* mod, void* fptr, const double * par_value);
void codegen_hf_add_with_coeff(llvm_module* mod, void* hfptr, double coeff, const double * par_value, double * data);
}

double codegen_f_evaluate(llvm_module* mod, void* fptr, const double * par_values){
    ParValues values(par_values, mod->get_parameters());
    return static_cast<theta::Function*>(fptr)->operator()(values);
}

void codegen_hf_add_with_coeff(llvm_module* mod, void* hfptr, double coeff, const double * par_values, double * data){
    ParValues values(par_values, mod->get_parameters());
    //const Histogram1D & h = static_cast<theta::HistogramFunction*>(hfptr)->operator()(values);
    //utils::add_fast_with_coeff(data, h.get_data(), coeff, h.size());
    static_cast<theta::HistogramFunction*>(hfptr)->apply_functor(add_to_vdouble(data, coeff), values);
}


/*  llvm_enable_function */
double llvm_enable_function::operator()(const theta::ParValues & values) const{
    std::vector<double> v_values(par_ids.size());
    size_t i=0;
    for(ParIds::const_iterator it=par_ids.begin(); it!=par_ids.end(); ++it, ++i){
        v_values[i] = values.get(*it);
    }
    return llvmf(&v_values[0]);
}

llvm_enable_function::llvm_enable_function(const theta::Configuration & cfg): f(PluginManager<theta::Function>::build(Configuration(cfg, cfg.setting["function"]))),
   m(f->get_parameters()){
    par_ids = f->get_parameters();
    llvm::Function * tmp_llvmf = create_llvm_function(f.get(), m, "enablef");
    llvmf = reinterpret_cast<t_function_evaluate>(m.getFunctionPointer(tmp_llvmf));
    if(llvmf==0){
        throw ConfigurationException("could not compile llvm");
    }
}

/*  llvm_enable_histogram_function */
void llvm_enable_histogram_function::apply_functor(const functor<Histogram1D> & f, const theta::ParValues & values) const{
    std::vector<double> v_values(par_ids.size());
    size_t i=0;
    for(ParIds::const_iterator it=par_ids.begin(); it!=par_ids.end(); ++it, ++i){
        v_values[i] = values.get(*it);
    }
    h0.set_all_values(0.0);    llvmhf(1.0, &v_values[0], h0.get_data());
    f(h0);
}

void llvm_enable_histogram_function::apply_functor(const functor<Histogram1DWithUncertainties> & f, const theta::ParValues & values) const{
    Histogram1D h;
    apply_functor(copy_to<Histogram1D>(h), values);
    f(Histogram1DWithUncertainties(h));
}

llvm_enable_histogram_function::llvm_enable_histogram_function(const theta::Configuration & cfg):
  hf(PluginManager<theta::HistogramFunction>::build(Configuration(cfg, cfg.setting["histogram_function"]))), m(hf->get_parameters()){
  par_ids = hf->get_parameters();
  size_t nbins;
  double xmin, xmax;
  hf->get_histogram_dimensions(nbins, xmin, xmax);
  h0 = Histogram1D(nbins, xmin, xmax);
  llvm::Function * tmp_llvmf = create_llvm_histogram_function(hf.get(), m, "enablehf");
  llvmhf = reinterpret_cast<t_hf_add_with_coeff>(m.getFunctionPointer(tmp_llvmf));
  if(llvmhf==0){
      throw ConfigurationException("could not compile llvm");
  }
}

REGISTER_PLUGIN(llvm_enable_function)
REGISTER_PLUGIN(llvm_enable_histogram_function)

    
