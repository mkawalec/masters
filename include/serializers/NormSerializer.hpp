#ifndef turb_NormSerializer_h
#define turb_NormSerializer_h

#include "Serializer.hpp"
#include "Integrator.hpp"

#include <fstream>

namespace turb {

    /*! \brief Outputs L2(u) to the output stream
     */
    class NormSerializer : public Serializer {
    public:
        NormSerializer();
        virtual ~NormSerializer() { unregister(); }

        void serialize(Integrator *instance, 
                std::ofstream *output, void *time);
        Serializer* clone() const { return new NormSerializer(*this); }
    };
}

#endif     
