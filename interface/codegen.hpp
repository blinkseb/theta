#ifndef CODEGEN_HPP
#define CODEGEN_HPP

#include "interface/decls.hpp"

#include <ostream>

namespace codegen{

typedef double (*t_function_evaluate)(const double *);
typedef double (*t_hf_add_with_coeff)(double, const double *, double *);
typedef double (*t_model_get_prediction)(const double *, double *);

void header(std::ostream & out, const theta::PropertyMap & pm, const theta::ParIds& pids, const theta::ObsIds & oids);
void footer(std::ostream & out);
size_t get_id(const theta::ParId & id);
size_t get_id(const theta::ObsId & id);
theta::ParId create_pid(size_t id);
theta::ObsId create_oid(size_t id);

// convert double to string, with guaranteed enough accuracy ...
std::string dtos(double d);

// returns full path to so file
void* compile_load_so(const std::string & cpp);

}

#endif

