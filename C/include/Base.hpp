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
        static T* choose(std::string name);

        static std::string list_available();
        virtual std::string additional_info() { return "";}

        std::string name;
        std::string description;
        static std::string class_name;

        virtual Base<T>* clone() const = 0;
    };
}

#include "Base.cpp"

#endif
