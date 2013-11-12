#ifndef turb_Serializer_h
#define turb_Serializer_h

#include "Integrator.hpp"
#include "Base.hpp"

#include <fstream>

namespace turb {

    class Serializer : public Base<Serializer> {
    public:
        virtual void serialize(Integrator *integrator, 
                std::ofstream *output, double time) = 0;

        virtual Serializer* clone() const = 0;
    };
}

#endif
