#ifndef turb_CMA_Simple_hpp
#define turb_CMA_Simple_hpp

#include "searchers/SimpleSearcher.hpp"

#include <armadillo>
using namespace arma;


namespace turb {

    class CMASimple : public SimpleSearcher {
    private:
        int N, lambda, mu;
        double sigma = 0.5,
               stop_fitness = 1,
               stop_iters,
               mueff,                   // Variance-effectiveness
               cc,                      // Time constant for cumulation for C
               cs,                      // t-const for cumulation for sigma control
               c1,                      // learning rate for rank-one update of C
               cmu,                     // for rank-mu update
               damps,                   // damping for sigma
               eigenval,
               chiN;

        vec **values,
            *weights,
            *D,
            *xmean,
            *xold,
            *pc,
            *ps;

        mat *B,
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

