#ifndef turb_DecayPathComputer_h
#define turb_DecayPathComputer_h

#include "computers/MultirunComputer.hpp"

#include <fstream>


namespace turb {

    class DecayPathComputer : public MultirunComputer<DecayPathComputer> {
    private:
        double decay_threshold, fast_threshold;

    public:
        DecayPathComputer();
        virtual ~DecayPathComputer() { unregister(this); }

        void compute_single(std::ofstream *output);
        Computer* clone() const { return new DecayPathComputer(*this); }
    };
}

#endif
