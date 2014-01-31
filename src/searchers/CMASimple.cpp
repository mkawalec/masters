#include "searchers/CMASimple.hpp"
#include "Searcher.hpp"
#include "Integrator.hpp"

#include <armadillo>
#include <vector>


namespace turb {

    CMASimple::CMASimple()
    {
        name = "CMA-simple";
        class_name = "Searcher";
        description = "Searches using the Covariant Matrix "
            "Adaptation";

        compatible_integrators.push_back("paper");

        Searcher::available.push_back(this);
    }

    CMASimple::~CMASimple()
    {

    }

    std::vector<double> CMASimple::run()
    {

    }


    void CMASimple::allocate(Integrator *integrator)
    {

    }

}
