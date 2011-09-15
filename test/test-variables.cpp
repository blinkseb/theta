#include "interface/random.hpp"
#include "interface/variables.hpp"

#include <boost/test/unit_test.hpp>

using namespace theta;
using namespace std;

BOOST_AUTO_TEST_SUITE(variables_tests)

struct varidtest{
    VarIdManager vm;
    const ParId par0;
    varidtest(): par0(vm.createParId("par0")){
    }
};


BOOST_AUTO_TEST_CASE(basic){
    varidtest v;
    BOOST_CHECK(v.vm.parNameExists("par0"));
    BOOST_CHECK(v.vm.getName(v.par0)=="par0");
    BOOST_CHECK(v.vm.getParId("par0")==v.par0);
}

BOOST_AUTO_TEST_CASE(getParIds){
    varidtest v;
    ParIds ids;
    ids.insert(v.par0);
    BOOST_CHECK(ids==v.vm.getAllParIds());
    
    ParId par1 = v.vm.createParId("par1");
    BOOST_REQUIRE(par1!=v.par0);
    BOOST_CHECK(not (ids==v.vm.getAllParIds()));
    
    ids.insert(par1);
    BOOST_CHECK(ids==v.vm.getAllParIds());
}

BOOST_AUTO_TEST_CASE(par_exceptions){
    varidtest v;
    
    //request variables not there:
    bool ex = false;
    try{
        v.vm.getParId("var1");
    }
    catch(NotFoundException &){
        ex = true;
    }
    BOOST_REQUIRE(ex);
    
    //create parameter already there
    ex = false;
    try{
        v.vm.createParId("par0");
    }
    catch(InvalidArgumentException &){
        ex = true;
    }
    BOOST_REQUIRE(ex);
    
    ParId par1 = v.vm.createParId("par1");
    BOOST_REQUIRE(v.par0!=par1);
    ParIds all_ids;
    all_ids.insert(par1);
    all_ids.insert(v.par0);
    BOOST_REQUIRE(all_ids == v.vm.getAllParIds());
}

BOOST_AUTO_TEST_CASE(basic_obs){
    VarIdManager vm;
    BOOST_REQUIRE(not vm.obsNameExists("obs0"));
    ObsId obs0 = vm.createObsId("obs0", 100, -0.2, 0.8);
    BOOST_CHECK(vm.getObsId("obs0")==obs0);
    BOOST_CHECK(vm.get_nbins(obs0)==100);
    BOOST_CHECK(vm.get_range(obs0).first==-0.2);
    BOOST_CHECK(vm.get_range(obs0).second==0.8);
}

BOOST_AUTO_TEST_CASE(all_obs){
    VarIdManager vm;
    ObsIds ids;
    BOOST_REQUIRE(vm.getAllObsIds()==ids);
    ObsId obs0 = vm.createObsId("obs0", 100, -0.8, 1.2);
    BOOST_REQUIRE(not (vm.getAllObsIds()==ids));
    ids.insert(obs0);
    BOOST_REQUIRE(vm.getAllObsIds()==ids);
}

BOOST_AUTO_TEST_CASE(exceptions_obs){
    VarIdManager vm;
    ObsIds ids;
    bool ex = false;
    try{
        vm.createObsId("obs0", 0, -1, 1);
    }
    catch(InvalidArgumentException &){
        ex = true;
    }
    BOOST_CHECK(ex);
    BOOST_REQUIRE(vm.getAllObsIds()==ids);
    
    ex = false;
    try{
        vm.createObsId("obs0", 100, 1, -1);
    }
    catch(InvalidArgumentException &){
        ex = true;
    }
    BOOST_CHECK(ex);
    BOOST_REQUIRE(vm.getAllObsIds()==ids);
}


BOOST_AUTO_TEST_CASE(parvalues_basic){
    ParValues vv;
    VarIdManager vm;
    ParId v0 = vm.createParId("v0");
    BOOST_REQUIRE(not vv.contains(v0));
    vv.set(v0, -7.0);
    BOOST_REQUIRE(vv.get(v0)==-7.0);
    BOOST_REQUIRE(vv.contains(v0));
    ParId v1 = vm.createParId("v1");
    bool ex = false;
    try{
        vv.get(v1);
    }
    catch(NotFoundException &){
        ex = true;
    }
    BOOST_REQUIRE(ex);
}

//as ParValues uses an internally growing vector, test with many parameters ...
BOOST_AUTO_TEST_CASE(parvalues_many){
    ParValues vv;
    VarIdManager vm;
    vector<ParId> parameters;
    ParValues values;
    for(size_t i=0; i<100; ++i){
        stringstream name;
        name << "parameter" << i;
        parameters.push_back(vm.createParId(name.str()));
        values.set(parameters.back(), 0.1*i + i*i);
    }
    for(size_t i=0; i<100; ++i){
        BOOST_REQUIRE_EQUAL(values.get(parameters[i]), 0.1*i + i*i);
    }
    
    ParValues values_sparse;
    values_sparse.set(parameters[90], 99.0);
    for(size_t i=0; i<100; ++i){
        if(i==90){
            BOOST_CHECK(values_sparse.get(parameters[i])==99.0);
        }
        else{
            BOOST_REQUIRE(not values_sparse.contains(parameters[i]));
        }
    }
    
    values.set(values_sparse);
    for(size_t i=0; i<100; ++i){
        if(i!=90){
            BOOST_REQUIRE_EQUAL(values.get(parameters[i]), 0.1*i + i*i);
        }
        else{
            BOOST_CHECK(values.get(parameters[i]) == values_sparse.get(parameters[i]));
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

