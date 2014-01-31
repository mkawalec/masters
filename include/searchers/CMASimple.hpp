#ifndef turb_CMA_Simple_hpp
#define turb_CMA_Simple_hpp

#include "searchers/SimpleSearcher.hpp"

#include <armadillo>
using namespace arma;


namespace turb {

    class CMASimple : public SimpleSearcher {
    private:
        int N, lambda, mu;
        double sigma = 0.3,
               stop_fitness = 1e-10,
               stop_iters = 100,
               mueff,                   // Variance-effectiveness
               cc,                      // Time constant for cumulation for C
               cs,                      // t-const for cumulation for sigma control
               c1,                      // learning rate for rank-one update of C
               cmu,                     // for rank-mu update
               damps,                   // damping for sigma
               eigenval,
               chiN;

        mat *weights,                   // array for weighted recombination
            *xmean,                     // initial point for variables
            *pc,
            *ps,
            *B,
            *D,
            *C,
            *invsqrtC;

        void set_params();

    public:
        CMASimple();
        ~CMASimple();

        std::vector<double> run();
        void allocate(Integrator *integrator);

        Searcher* clone() const { return new CMASimple(*this);}
    };
}

#endif

