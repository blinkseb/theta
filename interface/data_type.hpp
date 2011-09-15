#ifndef DATA_TYPE_HPP
#define DATA_TYPE_HPP

//note: data types are used by btoh database and produssource. In order
// not to generate cyclic dependencies, provide the definition in this separate file.

namespace theta{
    enum data_type { typeDouble, typeInt, typeString, typeHisto };
}

#endif

