#include "interface/histogram.hpp"
#include "interface/histogram-function.hpp"
#include "interface/phys.hpp"
#include "interface/utils.hpp"
#include "interface/random.hpp"
#include "interface/exception.hpp"


#include <boost/test/unit_test.hpp>
#include <iostream>


using namespace std;
using namespace theta;

void check_histos_equal(const Histogram & h1, const Histogram & h2){
    BOOST_REQUIRE(h1.get_nbins()==h2.get_nbins());
    BOOST_REQUIRE(h1.get_xmin()==h2.get_xmin());
    BOOST_REQUIRE(h1.get_xmax()==h2.get_xmax());
    //for the total weight, do not assume too much ...
    BOOST_REQUIRE(utils::close_to_relative(h1.get_sum_of_bincontents(),h2.get_sum_of_bincontents()));
    const size_t n = h1.get_nbins();
    for(size_t i=0; i<=n+1; i++){
        BOOST_CHECK(h1.get(i)==h2.get(i));
    }
}

// general note: use odd bin numbers to test SSE implementation for which
// an odd number of bins is a special case.

BOOST_AUTO_TEST_SUITE(histogram_tests)

//test construcors and copy assignment
BOOST_AUTO_TEST_CASE(ctest){
   //default construction:
   const size_t nbins = 101;
   Histogram m(nbins, -1, 1);
   BOOST_CHECK(m.get_nbins()==nbins);
   BOOST_CHECK(m.get_xmin()==-1);
   BOOST_CHECK(m.get_xmax()==1);
   BOOST_CHECK(m.get_sum_of_bincontents()==0);
   //fill a bit:
   for(size_t i=0; i<=nbins+1; i++){
       m.set(i, i*i);
   }

   //copy constructos:
   Histogram mcopy(m);
   BOOST_REQUIRE(mcopy.getData()!=m.getData());
   check_histos_equal(m, mcopy);

   //copy assignment:
   Histogram m200(2 * nbins, -1, 1);
   BOOST_CHECK(m200.get_nbins()==2*nbins);
   m200 = m;
   BOOST_CHECK(m200.get_nbins()==nbins);
   check_histos_equal(m200, m);

   //copy assignment with empty histo:
   Histogram h_empty;
   h_empty = m;
   check_histos_equal(h_empty, m);
   
   bool exception = false;
   try{
      Histogram m2(nbins, 1, 0);
   }
   catch(InvalidArgumentException & ex){
      exception = true;
   }
   BOOST_CHECK(exception);
   
   exception = false;
   try{
      Histogram m2(nbins, -1, -1);
   }
   catch(InvalidArgumentException & ex){
      exception = true;
   }
   BOOST_CHECK(exception);
}


//test get, set, fill:
BOOST_AUTO_TEST_CASE(getset){
    const size_t nbins=100;
    Histogram m(nbins, -1, 1);
    volatile double sum = 0.0;
    for(size_t i=0; i<=nbins; i++){
        volatile double a = sqrt(i+0.0);
        sum += a;
        m.set(i, a);
    }

    BOOST_CHECK(utils::close_to_relative(m.get_sum_of_bincontents(),sum));

    for(size_t i=0; i<=nbins; i++){
        volatile double a = sqrt(i+0.0);
        BOOST_CHECK(m.get(i) == a);
    }

    //fill:
    volatile double content = m.get(1);
    m.fill(-0.999, 1.7);
    content += 1.7;
    BOOST_CHECK(content==m.get(1));
    sum += 1.7;
    BOOST_CHECK(utils::close_to_relative(m.get_sum_of_bincontents(),sum));

    //fill in underflow:
    content = m.get(0);
    double delta = 10.032;
    m.fill(-1.001, delta);
    content += delta;
    BOOST_CHECK(content==m.get(0));
    sum += delta;
    BOOST_CHECK(utils::close_to_relative(m.get_sum_of_bincontents(),sum));

    //fill in overflow:
    content = m.get(nbins+1);
    delta = 7.032;
    m.fill(1.001, delta);
    content += delta;
    BOOST_CHECK(content==m.get(nbins+1));
    sum += delta;
    BOOST_CHECK(utils::close_to_relative(m.get_sum_of_bincontents(),sum));
}

//test +=
BOOST_AUTO_TEST_CASE(test_plus){
    Random rnd(new RandomSourceTaus());
    const size_t nbins = 101;
    Histogram m0(nbins, 0, 1);
    Histogram m1(m0);
    Histogram m_expected(m0);
    for(size_t i=0; i<=nbins+1; ++i){
        volatile double g0 = rnd.get();
        m0.set(i, g0);
        volatile double g1 = rnd.get();
        m1.set(i, g1);
        g0 += g1;
        m_expected.set(i, g0);
    }
    //m0+=m1 should add the histogram m1 to m0, leaving m1 untouched ...
    Histogram m1_before(m1);
    Histogram m0_before(m0);
    m0+=m1;
    //the sum of weights could be slightly off, but check_histos_equal does handle this.
    check_histos_equal(m0, m_expected);
    check_histos_equal(m1, m1_before);
    //... and it should commute:
    m1+=m0_before;
    check_histos_equal(m1, m_expected);
    //BOOST_CHECKPOINT("test_plus m1, m0");
    check_histos_equal(m1, m0);
}

//test *=
BOOST_AUTO_TEST_CASE(test_multiply){
    Random rnd(new RandomSourceTaus());
    const size_t nbins = 101;
    Histogram m0(nbins, 0, 1);
    Histogram m1(m0);
    Histogram m0m1(m0);
    Histogram m0factor_expected(m0);
    double factor = rnd.get();
    for(size_t i=0; i<=nbins+1; i++){
        double g0 = rnd.get();
        m0.set(i, g0);
        m0factor_expected.set(i, g0*factor);
        double g1 = rnd.get();
        m1.set(i, g1);
        m0m1.set(i, g0*g1);
    }
    Histogram m0_before(m0);
    m0*=factor;
    check_histos_equal(m0, m0factor_expected);
    m1*=m0_before;
    check_histos_equal(m1, m0m1);
    bool exception = false;
    //check error behaviour:
    try{
       Histogram m2(nbins + 1, 0, 1);
       m0 *= m2;
    }
    catch(InvalidArgumentException & ex){
       exception = true;
    }
    BOOST_REQUIRE(exception);
}

BOOST_AUTO_TEST_CASE(test_reset){
   Random rnd(new RandomSourceTaus());
   const size_t nbins = 101;
   Histogram m(nbins, rnd.uniform(), rnd.uniform() + 2.0);
   for(size_t i=0; i<=nbins+1; i++){
       m.set(i, 0.1*i);
   }
   Histogram m0=m;
   check_histos_equal(m0, m);
   
   //should not change nbins, xmin, xmax:   
   m.reset();
   BOOST_CHECK(m.get_xmin()==m0.get_xmin());
   BOOST_CHECK(m.get_xmax()==m0.get_xmax());
   BOOST_CHECK(m.get_nbins()==m0.get_nbins());
   BOOST_CHECK(m.get_sum_of_bincontents()==0.0);
   for(size_t i=0; i<=nbins+1; i++){
       BOOST_REQUIRE(m.get(i)==0.0);
   }
   
   bool exception = false;
   try{
      m.reset(nbins + 1);
   }
   catch(InvalidArgumentException & ex){
      exception = true;
   }
   BOOST_CHECK(exception);
}

BOOST_AUTO_TEST_SUITE_END()
