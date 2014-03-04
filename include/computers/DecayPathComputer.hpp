#ifndef turb_DecayPathComputer_h
#define turb_DecayPathComputer_h

#include "computers/MultirunComputer.hpp"
#include "Searcher.hpp"

#include <fstream>


namespace turb {

    /*! \brief Runs the integrator many times and
     *      output decay paths for the runs that decayed
     *      before a certain time.
     */
    class DecayPathComputer : public MultirunComputer<DecayPathComputer> {
    private:
        double decay_threshold, fast_threshold;
        Searcher *searcher = NULL;
        bool search;

    public:
        DecayPathComputer();
        virtual ~DecayPathComputer() { unregister(); delete searcher;}

        double compute_single(std::ofstream *output, MultirunComputer *base);

        Computer* clone() const {
            auto tmp = new DecayPathComputer(*this);
            tmp->is_clone = true;
            return tmp;
        }

        void set_options();
    };
}

#endif
