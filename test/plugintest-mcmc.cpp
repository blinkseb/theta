#include "interface/utils.hpp"
#include "interface/random.hpp"
#include "plugins/mcmc-result.hpp"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(mcmc_tests)

using namespace theta;
using namespace std;


//test the numerical stability of the covariance calculation of the
// mcmcResult class.
BOOST_AUTO_TEST_CASE(mcmc_result0){
    Result res(2);
    
    //check a minimally filled result: fill in two points with 0.0.
    double vec[2] = {0,1};
    res.fill(vec, 0, 1);
    res.fill(vec, 0, 1);
    res.end();
    BOOST_CHECK(res.getCount()==2);
    BOOST_CHECK(res.getnpar()==2);
    vector<double> means = res.getMeans();
    BOOST_REQUIRE(means.size()==2);
    BOOST_CHECK(means[0] == 0.0);
    BOOST_CHECK(means[1] == 1.0);
    Matrix cov = res.getCov();
    BOOST_CHECK(cov(0,0)==0.0);
    BOOST_CHECK(cov(1,1)==0.0);
    BOOST_CHECK(cov(1,0)==0.0);
    BOOST_CHECK(cov(0,1)==0.0);
}
    
BOOST_AUTO_TEST_CASE(mcmc_result1){
    Result res(2);
    double vec[2];
    //now fill in data in which first component is half
    // of the time 7.9, half of the time 12.1, i.e. 10 +- 2.1
    // and the second component is constant 1.1.
    res.reset();
    vec[1] = 1.1;
    const int N = 100000;
    for(int i=0; i<N; ++i){
        vec[0] = i%2?(10-2.1):(10+2.1);
        res.fill(vec, 0.0, 1);
    }
    res.end();
    vector<double> means = res.getMeans();
    BOOST_CHECK(utils::close_to_relative(means[0], 10.0));
    BOOST_CHECK(utils::close_to_relative(means[1], 1.1));
    Matrix cov = res.getCov();
    BOOST_CHECK(utils::close_to(cov(1,0), 0.0, 1.0));
    BOOST_CHECK(utils::close_to(cov(0,1), 0.0, 1.0));
    BOOST_CHECK(fabs(cov(1,1)) < 1e-15);
    BOOST_CHECK(utils::close_to_relative(sqrt(cov(0,0)), 2.1));
}

BOOST_AUTO_TEST_CASE(mcmc_result2){
    //fill using non-trivial weights ...
    Result res(2), res2(2);
    double vec[2];
    vec[1] = 1.1;
    double mean0 = 17.38;
    double width0 = 2.754;
    const int N = 100000;
    int res2_weight_odd=0, res2_weight_even=0;
    for(int i=0; i<N; ++i){
        //in result, fill alternating:
        vec[0] = mean0 + (i%2?width0:-width0);
        res.fill(vec, 0.0, 1);
        // in res2, fill with different weights:
        if(i%2){
            if(res2_weight_odd==N/2) continue;
            int w = min(N/2 - res2_weight_odd, i%10 + 1);
            res2_weight_odd += w;
            res2.fill(vec, 0.0, w);
        }
        else{
            if(res2_weight_even==N/2) continue;
            int w = min(N/2 - res2_weight_even, i%10 + 1);
            res2_weight_even += w;
            res2.fill(vec, 0.0, w);
        }
    }
    Matrix cov = res.getCov();
    Matrix cov2 = res2.getCov();
    for(size_t i=0; i<2; ++i){
        for(size_t j=0; j<2; ++j){
            BOOST_CHECK(utils::close_to(cov(i,j), cov2(i,j), 10.0));
        }
    }
/*    cout << setprecision(18);
    cout << "cov:" << endl << cov(0,0) << " " << cov(0, 1) << endl << cov(1,0) << " " << cov(1,1) << endl;
    cout << "cov2:" << endl << cov2(0,0) << " " << cov2(0, 1) << endl << cov2(1,0) << " " << cov2(1,1) << endl;*/
}


BOOST_AUTO_TEST_CASE(mcmc_result3){
    //similar to mcmc_result2: fill using non-trivial weights, but use
    // a bit more "complicated" data ...
    Result res(2), res2(2);
    Random rnd(new RandomSourceTaus());
    //rnd.set_seed(time(0));
    double vec[2];
    const int N = 1000;
    for(int i=0; i<N; ++i){
        double value0 = rnd.gauss();
        double value1 = rnd.gauss();
        vec[0] = value0 + 3 * value1 + 5.0;
        vec[1] = value0  - 2*value1 - 7.0;
        int weight = i%10 + 1;
        for(int k=0; k<weight; ++k){
           res.fill(vec, 0.0, 1);
        }
        res2.fill(vec, 0.0, weight);
    }
    Matrix cov = res.getCov();
    Matrix cov2 = res2.getCov();
    for(size_t i=0; i<2; ++i){
        for(size_t j=0; j<2; ++j){
            BOOST_CHECK(utils::close_to_relative(cov(i,j), cov2(i,j)));
        }
    }
    cout << setprecision(18);
}

BOOST_AUTO_TEST_SUITE_END()


