#ifndef turb_SimpleComputer_h
#define turb_SimpleComputer_h

#include "Serializer.hpp"
#include "Computer.hpp"
#include "helpers.hpp"

#include <thread>

namespace turb {
    
    class SimpleComputer : public Computer {
    protected:
        void compute();

    public:
        SimpleComputer();
        ~SimpleComputer() { unregister(this); }

        std::thread run();
    };
}
        


#endif
