#ifndef turb_MultirunComputer_h
#define turb_MultirunComputer_h

#include "Computer.hpp"

#include <fstream>

namespace turb {

    template <typename T>
    class MultirunComputer : public Computer {
    protected:
        void compute();

    public:
        virtual void compute_single(std::ofstream *output) = 0;

    };
}

#include "computers/MultirunComputer.cpp"

#endif
