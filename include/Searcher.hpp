#ifndef turb_Searcher_h
#define turb_Searcher_h

#include "Integrator.hpp"
#include <vector>

namespace turb {

    class Searcher {
    private:
        Integrator *integrator;

    public:
        Searcher(Integrator *integrator) : integrator(integrator) { };

        std::vector<double> run();
    };
}

#endif
