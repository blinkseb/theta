#ifndef PLUGINS_SECANT_HPP
#define PLUGINS_SECANT_HPP

#include "interface/exception.hpp"


#include <iomanip>
using namespace std;

/** \brief The secant method to find the root of a one-dimensional function
 *
 * \param x_low The lower end of the start interval
 * \param x_high The higher end of the start interval
 * \param x_accuracy If the found interval is shorter than this, the iteration will stop
 * \param f_x_low is function(x_low).Used to save one function evalutation.
 * \param f_x_high is function(x_high). Used to save one function evaluation.
 * \param f_accuracy If the absolute function value is lower than this, the iteration will stop
 * \param function The function object to use.
 *
 * The iteration will stop if either one if the x_accuracy and f_accuracy criteria is fulfilled.
 * To use only one criterion, set the other value to 0.0 to prevent it from being fulfilled.
 *
 * Note that the function values at x_low and x_high must have different sign. Otherwise,
 * an InvalidArgumentException will be thrown.
 * All x and function values mus be finite.
 */
template<typename T>
double secant(double x_low, double x_high, double x_accuracy, double f_x_low, double f_x_high, double f_accuracy, const T & function){
    assert(isfinite(x_low) && isfinite(x_high));
    assert(x_low <= x_high);
    assert(isfinite(f_x_low) && isfinite(f_x_high));
    if(f_x_low * f_x_high >= 0) throw theta::InvalidArgumentException("secant: function values have the same sign!");
    if(fabs(f_x_low) <= f_accuracy) return x_low;
    if(fabs(f_x_high) <= f_accuracy) return x_high;

    const double old_interval_length = x_high - x_low;    
    //calculate intersection point for secant method:
    double x_intersect = x_low - (x_high - x_low) / (f_x_high - f_x_low) * f_x_low;
    assert(x_intersect >= x_low);
    assert(x_intersect <= x_high);
    if(old_interval_length < x_accuracy){
        return x_intersect;
    }
    double f_x_intersect = function(x_intersect);
    double f_mult = f_x_low * f_x_intersect;
    //fall back to bisection if the new interval would not be much smaller:
    double new_interval_length = f_mult < 0 ? x_intersect - x_low : x_high - x_intersect;
    if(new_interval_length > 0.5 * old_interval_length){
        x_intersect = 0.5*(x_low + x_high);
        f_x_intersect = function(x_intersect);
        f_mult = f_x_low * f_x_intersect;
    }
    if(f_mult < 0){
        return secant(x_low, x_intersect, x_accuracy, f_x_low, f_x_intersect, f_accuracy, function);
    }
    else if(f_mult > 0.0){
        return secant(x_intersect, x_high, x_accuracy, f_x_intersect, f_x_high, f_accuracy, function);
    }
    //it can actually happen that we have 0.0. In this case, return the x value for
    // the smallest absolute function value:
    else{
        f_x_intersect = fabs(f_x_intersect);
        f_x_low = fabs(f_x_low);
        f_x_high = fabs(f_x_high);
        if(f_x_low < f_x_high && f_x_low < f_x_intersect) return x_low;
        if(f_x_high < f_x_intersect) return x_high;
        return x_intersect;
    }
}

#endif
