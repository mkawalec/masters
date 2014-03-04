#ifndef turb_SimpleComputer_h
#define turb_SimpleComputer_h

#include "Serializer.hpp"
#include "Computer.hpp"

namespace turb {

    /*! \brief Runs the integrator for one full cycle
     *      output data that is Serializer-dependent
     */
    class SimpleComputer : public Computer {
    protected:
        void compute();

    public:
        SimpleComputer();
        virtual ~SimpleComputer() { unregister(); }

        Computer* clone() const {
            auto tmp = new SimpleComputer(*this);
            tmp->is_clone = true;
            return tmp;
        }
    };
}



#endif
