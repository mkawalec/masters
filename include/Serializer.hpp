#ifndef turb_Serializer_h
#define turb_Serializer_h

#include "Integrator.hpp"
#include "Base.hpp"

#include <fstream>

namespace turb {

    /*! \brief Serializes the program data to some
     *      output stream
     *
     *  Every Serializer needs to be stateless,
     *  keeping any kind of state information inside
     *  a Serializer class will most probably result
     *  in nasty and hard to trace bugs
     */
    class Serializer : public Base<Serializer> {
    public:
        /*! \brief Function called to persist current
         *      integrator state.
         *
         *  \param integrator integrator instance from
         *      which data will be taken
         *  \param output output stream
         *  \param time current time, its role depends
         *      on output method used
         */
        virtual void serialize(Integrator *integrator, 
                std::ofstream *output, void *time) = 0;

        /*! \brief Base::clone overload to return
         *      a correct pointer type
         */
        virtual Serializer* clone() const = 0;
    };
}

#endif
