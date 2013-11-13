#ifndef turb_DecayMultirunComputer_h
#define turb_DecayMultirunComputer_h

#include "computers/MultirunComputer.hpp"
#include "Computer.hpp"

#include <fstream>

namespace turb {

    class DecayMultirunComputer : public MultirunComputer<DecayMultirunComputer> {
    private:
        double decay_threshold;

    public:
        DecayMultirunComputer();
        virtual ~DecayMultirunComputer() { unregister(this); }
        void compute_single(std::ofstream *output);

        Computer* clone() const { return new DecayMultirunComputer(*this); }
    };
}

#endif
