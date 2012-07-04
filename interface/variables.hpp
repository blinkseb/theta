#ifndef VARIABLES_HPP
#define VARIABLES_HPP

#include <set>
#include <map>
#include <string>
#include <vector>
#include <ostream>

#include "interface/decls.hpp"
#include "interface/exception.hpp"

#include <cmath>

#include <boost/utility.hpp>

namespace theta {
    
    /** \brief To refer to a certain parameter or observable,\c ParId and \c ObsId instances are used throughout %theta.
     *
     *  ParId and ObsId are used internally everyehere where a user would write the parameter / observable name
     *  to refer to a certain parameter / observable.
     *
     *  The \c VarId class defines equality and less-than order relation, so they can be used as key
     *  in a map.
     *
     *  Any associated informations such as default value and ranges are stored
     *  in a VarIdManager instance. Concrete parameter values are defined via an instance
     *  of the ParValues class.
     */
    template<typename tag>
    class VarId{
    friend class VarIdManager;
    friend class ParValues;
    template<typename T>  friend class DataT;
    friend std::ostream & operator<<(std::ostream & out, const VarId & vid){
        return out << vid.id;
    }
    public:
        //@{
        /** \brief Implements the order and equality semantics.
         */
        bool operator==(const VarId & rhs) const{
            return id==rhs.id;
        }
        bool operator!=(const VarId & rhs) const{
            return id!=rhs.id;
        }        
        bool operator<(const VarId & rhs) const{
            return id<rhs.id;
        }
        //@}
        
    private:
        size_t id;
        explicit VarId(size_t i): id(i){}
    };
    
    // An alternative to tag structs would be using common inheritance from
    // a (instantiable) VarId base class.
    //
    // However, using typedefs with tags instead will create truely new C++ types
    // and will make comparion of ObsId and VarId objects impossible, as it should
    // be (CCS, Item 14).
    /// \brief Empty tag struct to create ParId from the VarId template.
    struct ParIdTag{};
    /// \brief Empty tag struct to create ObsId from the VarId template.
    struct ObsIdTag{};
    
    /** Identification object for model parameters.
     * 
     * ParId instances are usually managed via a VarIdManager instance, where
     * additional information about each parameter is stored.
     * 
     * \sa VarId ObsId
     */    
    typedef VarId<ParIdTag> ParId;
    
    /** \brief Identification object for model observables.
     * 
     * ObsId objects are typically managed via a VarIdManager instance, where 
     * additional information about each observable is stored.
     * 
     * \sa VarId ParId
     */    
    typedef VarId<ObsIdTag> ObsId;
    
    /** A container for \c ParId or \c ObsId values. Set the template parameter id_type
     * to either \c ParId or \c ObsId
     *
     * The interface allows an STL-like iteration over the VarIds contained in this object.
     * The iteration will visit the elements in their natural order.
     * 
     * \param id_type The type to build the container for. Only ObsId and ParId are used here.
     */
    template<class id_type>
    class VarIds {
    public:
        /// \brief Sort of an iterator. Not necessarily one of the standard flavor, but can be used to go through all ids.
        typedef typename std::set<id_type>::const_iterator const_iterator;

        /// \brief Get the an iterator pointing to the first element. 
        const_iterator begin()const {
            return vars.begin();
        }
        
        /// \brief Get the an iterator pointing past the last element.
        const_iterator end() const {
            return vars.end();
        }

        /** \brief Insert a new id
         */
        void insert(const id_type & id) {
            vars.insert(id);
        }
        
        /** \brief Insert new ids.
         *
         * Insert [first, last) in this container.
         */
        void insert_all(const VarIds & vids) {
            vars.insert(vids.vars.begin(), vids.vars.end());
        }

        /** \brief Test whether an id is contained.
         * 
         * \param id The id object to test. 
         * \return Whether \c id is contained.
         */ 
        bool contains(const id_type & id)const {
            return vars.find(id) != vars.end();
        }

        /** \brief Test equality with other VarIds object.
         *
         * Two VarIds are the same if and only if the set of contained VarId s is the same.
         */
        bool operator==(const VarIds<id_type> & rhs) const{
            return vars == rhs.vars;
        }

        /// The number of contained ids
        size_t size() const {
            return vars.size();
        }
    private:
        std::set<id_type> vars;
    };
    
    /// \brief Template instantiation for a set of observables
    typedef VarIds<ObsId> ObsIds;
    /// \brief Template instantiation for a set of parameters
    typedef VarIds<ParId> ParIds;

    class ParValues;
    
    /** \brief Manager class for parameter and observable information
     *
     * This class provides methods to save the information given in the "parameters"
     * and "observables" setting groups.
     *
     * This class
     * <ul>
     * <li>keeps track of the association between parameter / observable names and ParId / ObsId objects</li>
     * <li>saves the configured range and binning (for observables) and the range / default value (for parameters)</li>
     * </ul>
     *
     * Note that there does not exist any global "current value" of a variable.
     */
    class VarIdManager {
        friend class ParValues;
    public:
        //@{
        /** \brief Creates a new parameter or observable ids (ParId, ObsId) and associates it with the given name.
         *
         * Observable names on themselves and parameter names on themselves must be unique;
         * using the same name for an observable and a parameter is allowed.
         *
         * \c type is the parameter type. The meaning is user-defined; the value of "type"
         * is just saved here and returned by get_type.
         *
         * If the name is already used for another parameter / observable, an InvalidArgumentException is thrown.
         * In case of nbins==0 or xmax < xmin, an InvalidArgumentException will be thrown.
         */
        ParId create_par_id(const std::string & name, const std::string & type = "par");
        ObsId create_obs_id(const std::string & name, size_t nbins, double xmin, double xmax);
        //@}
                
        //@{
        /** \brief Return the name of the given ParId or ObsId.
         *
         * If the id is not known, a std::invalid_argument is thrown.
         */
        std::string get_name(const ParId & id) const;
        std::string get_name(const ObsId & id) const;
        //@}
        
        /// Return the type of this parameter
        std::string get_type(const ParId & id) const;
        
        //@{
        /** \brief Return the number of bins and range for an observable identified by the ObsId id.
         */
        size_t get_nbins(const ObsId & id) const;
        const std::pair<double, double> & get_range(const ObsId & id) const;
        //@}
        
        //@{
        /** \brief Return the ParId / ObsId with the given name
         *
         * If the name is not known, a std::invalid_argument is thrown.
         */
        ParId get_par_id(const std::string & name) const;
        ObsId get_obs_id(const std::string & name) const;
        //@}
        
        //@{
        /** Returns all currently defined ParId or ObsId identifiers as ParIds or ObsIds
         */
        ParIds get_all_parameters() const;
        ObsIds get_all_observables() const;
        //@}
        
        /** \brief Create an empty VarIdManager with no registered variables.
         */
        VarIdManager(): next_pid_id(0), next_oid_id(0) {
        }

    private:
        //ParIds:
        std::map<size_t, std::string> pid_to_name;
        std::map<size_t, std::string> pid_to_type;
        std::map<std::string, size_t> name_to_pid;
        size_t next_pid_id;
        //ObsIds:
        std::map<size_t, std::string> oid_to_name;
        std::map<std::string, size_t> name_to_oid;
        std::map<size_t, std::pair<double, double> > oid_to_range;
        std::map<size_t, size_t> oid_to_nbins;
        size_t next_oid_id;
    };

    /** \brief A mapping-like class storing parameter values.
     * 
     * Conceptually, represents a mapping from ParId instances to double values.
     * It can only hold non-NAN values (a NAN-value is treated as non-existent). 
     *
     * Used to pass <b>actual</b> values of parameters to functions, as opposed to
     * passing a set of parameters, where parameter <i>identity</i> is sufficient.
     * For the latter, ParIds and ObsIds objects are used.
     */
    class ParValues {
    private:
        void fail_get(const ParId & pid) const;
    public:
        /** \brief Default constructor which creates an empty container.
        */
        ParValues():values(10, NAN){}
        
        /** \brief Constructor optimized for parameter information from \c vm.
         *
         * This is semantically equivalent to the default constructor. Using this
         * constructor makes possible some optimizations based on the total number of
         * parameters.
         */
        explicit ParValues(const VarIdManager & vm): values(vm.next_pid_id, NAN){}
        
        /** \brief Constructor initializing the values according to an array of doubles
         *
         * The resulting ParValues instance is initialized with the given data
         * by iterating over par_ids and using the value at that position from data.
         * This is the convention to convert array data to ParValues used in theta.
         *
         * Assumed that data contains (at least) par_ids.size() values. Otherwise,
         * behaviour is undefined.
         */
        ParValues(const double * data, const ParIds & par_ids){
           //to reallocate only once, find the maximum id:
           size_t s = 1;
           for(ParIds::const_iterator it=par_ids.begin(); it!=par_ids.end(); ++it){
               s = std::max<size_t>(s, it->id + 1);
           }
           values.resize(s, NAN);
           size_t i = 0;
           for(ParIds::const_iterator it=par_ids.begin(); it!=par_ids.end(); ++it, ++i){
               values[it->id] = data[i];
           }
        }
        
        /** \brief Constructor initializing the values according to an existing ParValues instance, only using the given parameters
         *
         * Set the values according to rhs, but only those parameters given in \c pids.
         */
        ParValues(const ParValues & rhs, const ParIds & pids): values(rhs.values.size(), NAN){
            for(ParIds::const_iterator it = pids.begin(); it!=pids.end(); ++it){
                values[it->id] = rhs.values[it->id];
            }
        }
        
        /** \brief Set a value.
         *
         * Set the value of the ParId \c pid to \c val. Setting a parameter 
         * to NAN means to delete it from the list (i.e., get() will throw a
         * invalid_argument for that parameter). This makes it possible to "clear" a parameter
         * from a VarValues instance after setting it.
         *
         * Returns a reference to this \c ParValues object to allow 
         * chaining calls like \c values.set(a, 0.0).set(b, 1.0).set(c, 2.0) ...
         *
         * \param pid Identified the parameter to assign a new value to.
         * \param val The new value for the parameter.
         */
        ParValues & set(const ParId & pid, double val){
            const size_t id = pid.id;
            if(id >= values.size()){
                values.resize(id+1, NAN);
            }
            values[id] = val;
            return *this;
         }
         
         /** \brief Set all values contained in rhs.
          *
          * This is equivalent to calling set(pid, val) for each pair (pid, val) contained in rhs.
          */
         void set(const ParValues & rhs){
             if(rhs.values.size() > values.size()){
                 values.resize(rhs.values.size(), NAN);
             }
             for(size_t i=0; i<rhs.values.size(); ++i){
                 if(std::isnan(rhs.values[i]))continue;
                 values[i] = rhs.values[i];
             }
         }

        /** \brief Retrieve the current value of a parameter.
         *
         *  Throws a invalid_argument if no value was set for \c pid in this \c ParValues.
         *
         *  \param pid The parameter for which the value should be returned.
         *  \return The current value for the parameter \c pid.
         */
        double get(const ParId & pid) const{
            double result = 0.0;
            const size_t id = pid.id;
            if(id >= values.size() || std::isnan(result = values[id])){
                //do failure outside this function to keep this function small to increase inlining probability
                fail_get(pid);
            }
            return result;
        }
        
        /** \brief Retrieve the current value of a parameter, returning a default if the parameter was not set.
         *
         *  \c pid is the parameter for which the value should be returned.
         *
         *  Returns the current value for the parameter \c pid, or the \c default_value if no value for \c pid is available.
         */
        double get(const ParId & pid, double default_value) const{
            const size_t id = pid.id;
            if(id < values.size() && !std::isnan(values[id])){
                return values[id];
            }
            else{
                return default_value;
            }
        }

        /// fast replacement for get which does not perform boundary checking
        double get_unchecked(const ParId & pid) const{
            return values[pid.id];
        }

        /** \brief Returns whether there is a value for \c pid.
         */
        bool contains(const ParId & pid) const{
            const size_t id = pid.id;
            return id < values.size() && !std::isnan(values[id]);
        }
        
        /** \brief Test whether all given parameters have a value set
         */
        bool contains_all(const ParIds & pids) const{
            for(ParIds::const_iterator it=pids.begin(); it!=pids.end(); ++it){
                const size_t id = it->id;
                if(id > values.size() || std::isnan(values[id])) return false;
            }
            return true;
        }
        
        /** \brief Clear all parameter values
         */
        void clear(){
            std::fill(values.begin(), values.end(), NAN);
        }

    private:
        //values are stored using the ParId.id as index
        std::vector<double> values;
    };

}

#endif 
