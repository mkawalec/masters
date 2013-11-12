#ifndef turb_ComponentSerializer_h
#define turb_ComponentSerializer_h

#include "Serializer.hpp"
#include "Integrator.hpp"

#include <fstream>

namespace turb {

    class ComponentSerializer : public Serializer {
    public:
        ComponentSerializer();
        ~ComponentSerializer() { unregister(this); }

        void serialize(Integrator *instance, 
                std::ofstream *output, double time);

        Serializer* clone() const { return new ComponentSerializer(*this); }
    };
}

#endif     
