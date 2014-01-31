#ifndef turb_CMA_Simple_hpp
#define turb_CMA_Simple_hpp

#include "searchers/SimpleSearcher.hpp"


namespace turb {

    class CMASimple : public SimpleSearcher {

    public:
        CMASimple();
        ~CMASimple();

        std::vector<double> run();
        void allocate(Integrator *integrator);

        Searcher* clone() const { return new CMASimple(*this);}
    };
}

#endif

