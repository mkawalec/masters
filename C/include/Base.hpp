#ifndef turb_Base_h
#define turb_Base_h

#include <list>
#include <string>

namespace turb {

    template <typename T>
    void unused(T &&) { }

    template <typename T>
    class Base {
    protected:
        void unregister(T *instance);

    public:
        static std::list<T*> available;
        static std::string list_available();
        static T* choose(std::string name);

        std::string name;
        std::string description;
        std::string class_name;
    };
}

#include "Base.cpp"

#endif
