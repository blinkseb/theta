#include "interface/histogram.hpp"
#include "interface/histogram-function.hpp"
#include "interface/phys.hpp"
#include "test/utils.hpp"
#include "interface/random.hpp"
#include "interface/exception.hpp"


#include <boost/test/unit_test.hpp>
#include <iostream>


using namespace std;
using namespace theta;

bool histos_equal(const Histogram1D & h1, const Histogram1D & h2){
    if(h1.get_nbins()!=h2.get_nbins()) return false;
    if(h1.get_xmin()!=h2.get_xmin()) return false;
    if(h1.get_xmax()!=h2.get_xmax()) return false;
    const size_t n = h1.get_nbins();
    for(size_t i=0; i<n; i++){
        if(h1.get(i)!=h2.get(i)) return false;
    }
    return true;
}

// general note: use odd bin numbers to test SSE implementation for which
// an odd number of bins is a special case.

BOOST_AUTO_TEST_SUITE(histogram_tests)

//test constructors and copy assignment
BOOST_AUTO_TEST_CASE(ctest){
   //default construction:
    Histogram1D h_def;
    BOOST_CHECK(h_def.get_nbins()==0);
    BOOST_CHECK(h_def.get_data()==0);
    
   const size_t nbins = 101;
   Histogram1D m(nbins, -1, 1);
   BOOST_CHECK(m.get_nbins()==nbins);
   BOOST_CHECK(m.get_xmin()==-1);
   BOOST_CHECK(m.get_xmax()==1);
   BOOST_CHECK(m.get_sum()==0);
   //fill a bit:
   for(size_t i=0; i<nbins; i++){
       m.set(i, i*i);
   }

   //copy constructor:
   Histogram1D mcopy(m);
   BOOST_REQUIRE(mcopy.get_data()!=m.get_data());
   BOOST_CHECK(histos_equal(m, mcopy));

   //copy assignment:
   Histogram1D m200(2 * nbins, -1, 1);
   BOOST_CHECK(m200.get_nbins()==2*nbins);
   m200 = m;
   BOOST_CHECK(m200.get_nbins()==nbins);
   BOOST_CHECK(histos_equal(m200, m));

   //copy assignment with empty histo:
   Histogram1D h_empty;
   h_empty = m;
   BOOST_CHECK(histos_equal(h_empty, m));
   
   Histogram1D h_empty2;
   BOOST_CHECK(h_empty2.get_nbins()==0);
   BOOST_CHECK(h_empty2.get_data()==0);
   h_empty = h_empty2;
   BOOST_CHECK(h_empty.get_nbins()==0);
   BOOST_CHECK(h_empty.get_data()==0);
   
   bool exception = false;
   try{
      Histogram1D m2(nbins, 1, 0);
   }
   catch(invalid_argument & ex){
      exception = true;
   }
   BOOST_CHECK(exception);
   
   exception = false;
   try{
      Histogram1D m2(nbins, -1, -1);
   }
   catch(invalid_argument & ex){
      exception = true;
   }
   BOOST_CHECK(exception);
}


//test get, set, fill:
BOOST_AUTO_TEST_CASE(getset){
    const size_t nbins=100;
    Histogram1D m(nbins, -1, 1);
    volatile double sum = 0.0;
    for(size_t i=0; i<nbins; i++){
        double a = sqrt(i+0.0);
        sum += a;
        m.set(i, a);
    }

    BOOST_CHECK(close_to_relative(m.get_sum(),sum));

    for(size_t i=0; i<nbins; i++){
        double a = sqrt(i+0.0);
        BOOST_CHECK(m.get(i) == a);
    }

    //fill:
    double content = m.get(0);
    m.fill(-0.999, 1.7);
    content += 1.7;
    BOOST_CHECK(content==m.get(0));
    sum += 1.7;
    BOOST_CHECK(close_to_relative(m.get_sum(),sum));

    //fill in underflow, content should not change
    content = m.get(0);
    double delta = 10.032;
    m.fill(-1.001, delta);
    BOOST_CHECK(content==m.get(0));
    BOOST_CHECK(close_to_relative(m.get_sum(),sum));

    //fill in overflow, content should not change:
    content = m.get(nbins-1);
    delta = 7.032;
    m.fill(1.001, delta);
    BOOST_CHECK(content==m.get(nbins-1));
    BOOST_CHECK(close_to_relative(m.get_sum(),sum));
}

//test +=
BOOST_AUTO_TEST_CASE(test_plus){
    std::auto_ptr<RandomSource> rnd_src(new RandomSourceTaus());
    Random rnd(rnd_src);
    const size_t nbins = 101;
    Histogram1D m0(nbins, 0, 1);
    Histogram1D m1(m0);
    Histogram1D m_expected(m0);
    for(size_t i=0; i<nbins; ++i){
        volatile double g0 = rnd.get();
        m0.set(i, g0);
        volatile double g1 = rnd.get();
        m1.set(i, g1);
        g0 += g1;
        m_expected.set(i, g0);
    }
    //m0+=m1 should add the histogram m1 to m0, leaving m1 untouched ...
    Histogram1D m1_before(m1);
    Histogram1D m0_before(m0);
    m0+=m1;
    //the sum of weights could be slightly off, but check_histos_equal does handle this.
    BOOST_CHECK(histos_equal(m0, m_expected));
    BOOST_CHECK(histos_equal(m1, m1_before));
    //... and it should commute:
    m1+=m0_before;
    BOOST_CHECK(histos_equal(m1, m_expected));
    //BOOST_CHECKPOINT("test_plus m1, m0");
    BOOST_CHECK(histos_equal(m1, m0));
}

//test *=
BOOST_AUTO_TEST_CASE(test_multiply){
    std::auto_ptr<RandomSource> rnd_src(new RandomSourceTaus());
    Random rnd(rnd_src);
    const size_t nbins = 101;
    Histogram1D m0(nbins, 0, 1);
    Histogram1D m1(m0);
    Histogram1D m0m1(m0);
    Histogram1D m0factor_expected(m0);
    double factor = rnd.get();
    for(size_t i=0; i<nbins; i++){
        double g0 = rnd.get();
        m0.set(i, g0);
        m0factor_expected.set(i, g0*factor);
        double g1 = rnd.get();
        m1.set(i, g1);
        m0m1.set(i, g0*g1);
    }
    Histogram1D m0_before(m0);
    m0*=factor;
    BOOST_CHECK(histos_equal(m0, m0factor_expected));
    m1*=m0_before;
    BOOST_CHECK(histos_equal(m1, m0m1));
    bool exception = false;
    //check error behaviour:
    try{
       Histogram1D m2(nbins + 1, 0, 1);
       m0 *= m2;
    }
    catch(invalid_argument & ex){
       exception = true;
    }
    BOOST_REQUIRE(exception);
}


BOOST_AUTO_TEST_SUITE_END()
