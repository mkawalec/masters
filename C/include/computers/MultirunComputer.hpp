#ifndef turb_MultirunComputer_h
#define turb_MultirunComputer_h

#include "Computer.hpp"

#include <fstream>

namespace turb {

    class MultirunComputer : public Computer {
    protected:
        virtual void run_single(std::ofstream *output) = 0;
        void compute();

    };
}

#endif
