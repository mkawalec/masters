#ifndef turb_SimpleComputer_h
#define turb_SimpleComputer_h

#include "Computer.hpp"
#include "helpers.hpp"
#include "Serializer.hpp"

#include <thread>

namespace turb {
    
    class SimpleComputer : public Computer {
    private:
        size_t iteration;

    public:
        SimpleComputer();
        SimpleComputer(Serializer *serializer);
        ~SimpleComputer() { unregister(this); }

        std::thread run();
    };
}
        


#endif
