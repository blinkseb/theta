#include "interface/histogram.hpp"
#include "interface/phys.hpp"
#include "interface/variables.hpp"

#include <boost/test/unit_test.hpp>

using namespace theta;
using namespace std;

BOOST_AUTO_TEST_SUITE(data_tests)


BOOST_AUTO_TEST_CASE(basic){
    VarIdManager vm;
    ObsId oid1 = vm.createObsId("oid1", 10, 0, 1);
    ObsId oid2 = vm.createObsId("oid2", 10, 0, 1);
    ObsId oid3 = vm.createObsId("oid3", 10, 0, 1);
    ObsId oid4 = vm.createObsId("oid4", 10, 0, 1);
    ObsId oid5 = vm.createObsId("oid5", 10, 0, 1);
    
    
    Histogram h(10, 0, 1);
    Histogram h_invalid;

    Data d;
    ObsIds oids = d.getObservables();
    BOOST_CHECK(oids == ObsIds());
    BOOST_CHECK(d[oid1].get_nbins()==0);
    BOOST_CHECK(d[oid2].get_nbins()==0);
    BOOST_CHECK(d[oid3].get_nbins()==0);
    BOOST_CHECK(oids == ObsIds());
    
    oids.insert(oid1);
    d[oid1] = h;
    BOOST_CHECK(oids == d.getObservables());
    d[oid3] = h;
    oids.insert(oid3);
    BOOST_CHECK(oids == d.getObservables());
    
    //insert an oid not been there before:
    d[oid1] = h_invalid;
    d[oid3] = h_invalid;
    oids = ObsIds();
    BOOST_CHECK(oids == d.getObservables());
    d[oid5] = h;
    oids.insert(oid5);
    BOOST_CHECK(oids == d.getObservables());
    d[oid5] = h_invalid;
    BOOST_CHECK(ObsIds() == d.getObservables());
    
    const Data & c_d = d;
    bool ex = false;
    try{
       Histogram h3 = c_d[oid5];
    }
    catch(Exception &){
      ex = true;
    }
    BOOST_CHECK(ex);
}


BOOST_AUTO_TEST_SUITE_END()

