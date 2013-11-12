#ifndef turb_Serializer_h
#define turb_Serializer_h

#include "Integrator.hpp"
#include "helpers.hpp"

#include <fstream>

namespace turb {

    class Serializer : public Base<Serializer> {
    public:
        virtual void serialize(Integrator *integrator, 
                std::ofstream *output, double time) = 0;
    };
}

#endif
