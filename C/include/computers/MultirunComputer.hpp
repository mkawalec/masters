#ifndef turb_MultirunComputer_h
#define turb_MultirunComputer_h

#include "Computer.hpp"

#include <fstream>
#include <vector>

namespace turb {

    template <typename T>
    class MultirunComputer : public Computer {
    private:
        void fit_it(std::vector<double> *decay_times);
        double fit_part = 0.3;

    protected:
        void compute();

    public:
        virtual double compute_single(std::ofstream *output) = 0;

    };
}

#include "computers/MultirunComputer.cpp"

#endif
