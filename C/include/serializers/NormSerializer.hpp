#ifndef turb_NormSerializer_h
#define turb_NormSerializer_h

#include "Serializer.hpp"
#include "Integrator.hpp"
#include "helpers.hpp"

#include <fstream>

namespace turb {

    class NormSerializer : public Serializer {
    public:
        NormSerializer();
        ~NormSerializer() { unregister(this); }

        void serialize(Integrator *instance, 
                std::ofstream *output, double time);
    };
}

#endif     
