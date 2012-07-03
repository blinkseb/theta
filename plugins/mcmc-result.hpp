#ifndef PLUGINS_MCMC_RESULT_HPP
#define PLUGINS_MCMC_RESULT_HPP

#include <algorithm>
#include <vector>
#include "interface/matrix.hpp"
#include "interface/histogram.hpp"

namespace theta {

    /** \brief Class to be used as second argument to metropolisHastings to save some information of the Markov-Chain
     *
     * Also useful as base class for more specialized classes for use with metropolisHastings.
     * These derived classes should implement fill2 which will be called by fill and
     * does nothing in this class.
     */
    class Result {
    protected:
        // number of parameters of the likelihood function
        size_t npar;
        
        // number of total points in the chain (including rejected proposal points)
        size_t count;
        
        // number of different points in the chain, i.e., not counting rejected proposal points
        size_t count_different_points;
        
        // sliding mean of the parameter values in the chain
        std::vector<double> means;
        
        // sliding covariance times count
        Matrix count_covariance;
        
        virtual void fill2(const double * p, double nll, size_t weight){}
    public:
        /** \brief Construct result with \c npar parameters
         */
        Result(size_t npar);
        
        /// Declare destructor virtual to allow polymorphic access to derived classes
        virtual ~Result(){}
        
        /** \brief fill a new chain point with the given parameter values, nll value and weight
         *
         * This method is called by the metropolisHastings routine.
         * \param p contains getnpar() parameter values of the point
         * \param nll is the negative logarithm of the likelihood / posterior
         * \param weight is the weight of the point, i.e., the number of rejected proposals to jump away from it, plus one.
         */
        void fill(const double * p, double nll, size_t weight);
        
        /// Returns the number of parameters specified in the constructor
        size_t getnpar() const;
        
        /// Returns the number of points in the chain, also counting rejected proposal points
        size_t getCount() const;
        
        /// Returns the number of different point in the chain, i.e., not including rejected proposals
        size_t getCountDifferent() const;
        
        /// Returns the mean of the parameter values in the chain
        std::vector<double> getMeans() const;
        
        /// Returns the covariance matrix of the parameter values in the chain
        Matrix getCov() const;
    };

}

#endif
