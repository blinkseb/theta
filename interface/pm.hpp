#ifndef PM_HPP
#define PM_HPP

#include <boost/shared_ptr.hpp>
#include <map>
#include <string>
#include <typeinfo>

#include <sstream>

#include "interface/exception.hpp"

/** \brief container to save shared pointers of arbitrary type indexed by an arbitratry string and type
 *
 * Elements in the container are indexed by both a supplied instance name and the type
 * passed to the get and set function templates, i.e., both the name and the type have to be the same.
 *
 * Setting an instance_name to a NULL pointer is in effect deleting the element: the associated internal shared_ptr
 * instance is released and a subsequent call to \c get will throw an exception. There is no way to
 * distinguish whether a NULL pointer has been set for an instance name or nothing has been set at all.
 */
class PropertyMap{
public:

   /** \brief Retrieve a value previously set
    *
    * throws a std::exception if no value is stored with this combination of instance name and type T.
    */
   template<typename T>
   boost::shared_ptr<T> get(const std::string & instance_name = "default") const;
   
   /** \brief Store a value
    *
    * Stores (and possibly overwrites) a value using the instance_name and gievn type T.
    */
   template<typename T>
   void set(const std::string & instance_name, const boost::shared_ptr<T> & value);
   
   template<typename T>
   bool exists(const std::string & instance_name) const;
   
   virtual ~PropertyMap(){}
private:
   struct nametype{
       std::string name;
       const std::type_info & type;
       nametype(const std::string & name_, const std::type_info & type_):name(name_), type(type_){}
       bool operator<(const nametype & rhs) const{
           return (name < rhs.name) || (name == rhs.name && type.before(rhs.type));
       }
   };
   std::map<nametype, boost::shared_ptr<void> > instances;
};


template<typename T>
boost::shared_ptr<T> PropertyMap::get(const std::string & instance_name) const {
    nametype nt(instance_name, typeid(T));
    std::map<nametype, boost::shared_ptr<void> >::const_iterator it=instances.find(nt);
    if(it==instances.end()){
        std::stringstream ss;
        ss << "PropertyMap: instance with name '" << instance_name << "' with type " << typeid(T).name() << " not found";
        throw theta::NotFoundException(ss.str());
    }
    return boost::static_pointer_cast<T>(it->second);
}

template<typename T>
bool PropertyMap::exists(const std::string & instance_name) const {
    nametype nt(instance_name, typeid(T));
    return instances.find(nt) != instances.end();
}

template<typename T>
void PropertyMap::set(const std::string & instance_name, const boost::shared_ptr<T> & value){
    nametype nt(instance_name, typeid(T));
    if(value.get()) instances[nt] = value;
    else instances.erase(nt);
}

#endif

