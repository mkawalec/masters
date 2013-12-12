#ifndef turb_DecayPathComputer_h
#define turb_DecayPathComputer_h

#include "computers/MultirunComputer.hpp"

#include <fstream>


namespace turb {

    /*! \brief Runs the integrator many times and 
     *      output decay paths for the runs that decayed
     *      before a certain time.
     */
    class DecayPathComputer : public MultirunComputer<DecayPathComputer> {
    private:
        double decay_threshold, fast_threshold;

    public:
        DecayPathComputer();
        virtual ~DecayPathComputer() { unregister(); }

        double compute_single(std::ofstream *output);
        Computer* clone() const { return new DecayPathComputer(*this); }

        void set_options();
    };
}

#endif
