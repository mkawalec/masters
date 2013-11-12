#ifndef turb_Base_cpp
#define turb_Base_cpp

#include "exceptions.hpp"

#include <string>
#include <list>

namespace turb {

    template <typename T>
    std::list<T*> Base<T>::available;

    template <typename T>
    std::string Base<T>::class_name;

    template <typename T>
    void Base<T>::unregister(T *instance)
    {
        Base::available.remove(instance);
    }

    template <typename T>
    T* Base<T>::choose(std::string name)
    {
        for(typename std::list<T*>::iterator it = T::available.begin();
            it != T::available.end(); ++it) {
            if ((*it)->name == name) {
                return *it;
            }
        }

        throw InstanceNotFound(T::class_name + " " + name + " is not currently available");
    }

    template <typename T>
    std::string Base<T>::list_available()
    {
        std::string output_methods_desc;
        
        for (typename std::list<T*>::iterator it = T::available.begin();
             it != T::available.end(); ++it) {
            output_methods_desc += (*it)->name + ": \t" +
                (*it)->description + "\n\n";
        }

        return output_methods_desc;
    }
}


#endif 
