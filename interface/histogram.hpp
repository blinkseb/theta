#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

#include "interface/decls.hpp"
#include "interface/utils.hpp"
#include <cstring>

void get_allocs_frees(int & n_alls, int & n_frs);

namespace theta{

/** A Histogram class holding binned data.
 *
 * Only equally spaced binning is supported. Binning is therefore fully determined
 * by a range of double (\c xmin, \c xmax) and the number of bins, \c nbins.
 *
 * Bin index convention:
 * - bin 0 is the "underflow bin", it corresponds to the range (-infinity, xmin)
 * - bin \c i with \c i &gt;=1 and \c i &lt;= nbins holds weight for the range [xmin + (i-1)*binwidth, xmin + i*binwidth), where binwidth = (xmax-xmin)/nbins.
 * - bin nbins+1 is the "overflow bin", it corresponds to the range [xmax, infinity).
 */
class Histogram {
private:
    double* histodata;
    //number of bins for the histogram, *excluding* underflow and overflow:
    size_t nbins;
    double xmin, xmax;
    
    void initFromHisto(const Histogram & h);
    
    // throws the actual exception, with the message.
    void fail_check_compatibility(const Histogram & h) const;
    
public:
    /** \brief Create an empty Histogram with \c bins bins with range (\c xmin, \c xmax )
     */
    explicit Histogram(size_t bins=0, double xmin=0, double xmax=1);

    /// Copy constructor. Copies the contents of \c rhs into this.
    Histogram(const Histogram& rhs);

    /// Copy assignment. Copies the content of \c rhs into this.
    void operator=(const Histogram& rhs);

    ///Destructor. De-allocates internal memory of the histogram data
    ~Histogram();

    /** \brief Set all entries to zero
     *
     * If \c nbins!=0, also change the number of bins and range to the parameters given. Otherwise,
     * the binning and range remain unchanged and bin contents is just reset to zero.
     *
     * If nbins &gt; 0 is specified, you also <em>must</em> supply valid values for \c xmin and \c xmax.
     */
    void reset(size_t nbins=0, double xmin=0, double xmax=0);

    /** \brief Set all bin entries to 1.0
     *
     * This is mainly useful for resetting Histograms subsequently used for multiplication.
     */
    void reset_to_1();

    /** \brief Returns the bincontents of bin i
     *
     * See bin index convention in the class documentation. This function does no range checking.
     */
    double get(size_t i) const{
        return histodata[i];
    }

    /** \brief Set bin entry i to weight.
     *
     * See bin index convention in the class documentation. This function does no range checking.
     */
    void set(size_t i, double weight){
       histodata[i] = weight;
    }

    /** \brief Return a pointer to the raw Histogram data
     *
     * The bin index convention applies to the result the same as for the get() and set() routines,
     * i.e., get(i) is the same as getData()[i].
     */
    const double* getData() const{
        return histodata;
    }
    double* getData() {
        return histodata;
    }

    /// Get the number of bins of this Histogram
    size_t get_nbins() const{
       return nbins;
    }

    /// Get the minimum x value for this Histogram
    double get_xmin() const{
       return xmin;
    }

    /// Get the maximum x value of this Histogram
    double get_xmax() const{
       return xmax;
    }

    /** \brief Add weight to the bin corresponding to xvalue
     *
     * In case of xvalue &lt; xmin, \c weigt is added to the underflow bin 0.
     * If xvalue &gt; xmax, \c weight is added to bin nbins+1.
     */
    void fill(double xvalue, double weight);

    /// Get the sum of bin contents of all bins of the Histogram.
    double get_sum_of_bincontents() const{
        double result = 0.0;
        for(size_t i=0; i<=nbins+1; ++i){
            result += histodata[i];
        }
        return result;
    }

    /** \brief Returns the x value for the bin center of bin ibin
     *
     *  The returned value is only valid for ibin>=1 && ibin&lt;=nbins
     */
    double get_bincenter(int ibin) const{
        return xmin + (ibin-0.5) * (xmax - xmin) / nbins;
    }

    /** \brief Multiply the Histogram \c other with this, bin-by-bin
     *
     * Throws an \c InvalidArgumentException if \c other is not compatible with this Histogram.
     * \c &other must be different from \c this.
     */
    void operator*=(const Histogram & other);

    /// Multiply the entry of each bin by \c a.
    void operator*=(const double a){
       utils::mul_fast(histodata, a, nbins+2);
    }

    /** \brief check compatibility of \c this to the \c other Histogram.
     *
     * A Histogram is considered compatible if it has the exact same number of bins and range. In case if incompatibility, an
     * \c InvalidArgumentException is thrown.
     */
    //note: this is split into checking and reporting to increase inlining probability
    void check_compatibility(const Histogram & h) const{
        if (nbins != h.nbins || xmin!=h.xmin || xmax != h.xmax){
           fail_check_compatibility(h);
        }
    }
    
    /** \brief Calculate this = this * (nominator/denominator)^exponent, bin-by-bin.
    *
    * An \c InvalidArgumentException is thrown if either \c nominator or \c denominator are not compatible with this.
    */
    void multiply_with_ratio_exponented(const Histogram & nominator, const Histogram & denominator, double exponent);
    
    /** \brief Add another Histogram to this
     *
     * Throws an \c InvalidArgumentException if \c other is not compatible to this.
     */
    void operator+=(const Histogram & other){
        utils::add_fast(histodata, other.histodata, nbins+2);
    }

    /** \brief Calculate this = this + coeff * other.
     *
     * Throws an \c InvalidArgumentException if \c other is not compatible with this.
     */
    void add_with_coeff(double coeff, const Histogram & other){
       utils::add_fast_with_coeff(histodata, other.histodata, coeff, nbins+2);
    }
    
};

}

#endif

