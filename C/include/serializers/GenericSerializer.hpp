#ifndef turb_GenericSerializer_h
#define turb_GenericSerializer_h

#include "Serializer.hpp"

#include <string>
#include <fstream>

namespace turb {

    class GenericSerializer : public Serializer {
    public:
        GenericSerializer();
        virtual ~GenericSerializer() { unregister(); }

        void serialize(Integrator *nothing,
                std::ofstream *output, void *output_data);
        Serializer* clone() const { return new GenericSerializer(*this); }
    };
}

#endif
