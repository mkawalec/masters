#ifndef turb_Serializer_h
#define turb_Serializer_h

#include "Integrator.hpp"

#include <memory>
#include <fstream>
#include <list>
#include <string>

namespace turb {

    class Serializer {
    protected:
        void unregister(Serializer *instance);

    public:
        static std::list<Serializer*> available;
        std::string name;
        std::string description;

        static Serializer* choose_serializer(std::string name);

        virtual void serialize(Integrator *integrator, 
                std::ofstream *output, double time) = 0;
    };
}

#endif
