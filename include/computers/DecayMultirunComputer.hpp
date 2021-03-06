#ifndef turb_DecayMultirunComputer_h
#define turb_DecayMultirunComputer_h

#include "computers/MultirunComputer.hpp"
#include "Computer.hpp"

#include <fstream>

namespace turb {

    /*! \brief Runs the integrator fully many times
     *      and output the decay times for each run
     */
    class DecayMultirunComputer : public MultirunComputer<DecayMultirunComputer> {
    private:
        double decay_threshold;

    public:
        DecayMultirunComputer();
        virtual ~DecayMultirunComputer() { unregister(); }
        double compute_single(std::ofstream *output, MultirunComputer *base);

        Computer* clone() const {
            auto tmp = new DecayMultirunComputer(*this);
            tmp->is_clone = true;
            return tmp;
        }
    };
}

#endif
