#ifndef turb_ComponentSerializer_h
#define turb_ComponentSerializer_h

#include "Serializer.hpp"
#include "Integrator.hpp"

#include <fstream>

namespace turb {

    /*! \brief Dumps Fourier components of u
     *      to the specified output file.
     */
    class ComponentSerializer : public Serializer {
    public:
        ComponentSerializer();
        virtual ~ComponentSerializer() { unregister(); }

        void serialize(Integrator *instance, 
                std::ofstream *output, void *time);

        Serializer* clone() const { return new ComponentSerializer(*this); }
    };
}

#endif     
