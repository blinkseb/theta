#ifndef THETA_DECLS_HPP
#define THETA_DECLS_HPP

// this header is included in almost any other header which require the declaration
// of classes. Whenever possible, only this header should be included.

namespace libconfig{
    class Setting;
}

namespace theta{
    class Histogram;
    class HistogramFunction;
    class Minimizer;
    class Producer;
    class ProductsSource;
    class ProductsSource;
    class Distribution;
    
    //database.hpp
    class Database;
    class DatabaseInput;
    class Table;
    class Column;
    class ProductsTable;
        
    namespace plugin{
          class Configuration;
    }
    class Random;
    class Run;
    
    //variables.hpp
    class VarIdManager;
    template<typename tag> class VarId;
    struct ParIdTag;
    struct ObsIdTag;
    typedef VarId<ParIdTag> ParId;
    typedef VarId<ObsIdTag> ObsId;
    template<class id_type> class VarIds;
    typedef VarIds<ObsId> ObsIds;
    typedef VarIds<ParId> ParIds;
    class ParValues;
    
    //phys.hpp
    class Data;
    class Function;
    class Model;
    class DataSource;
    class NLLikelihood;
    
    class SettingWrapper;
}

#endif
