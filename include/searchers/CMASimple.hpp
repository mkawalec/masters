#ifndef turb_CMA_Simple_hpp
#define turb_CMA_Simple_hpp

#include "searchers/SimpleSearcher.hpp"

#include <armadillo>
#include <iostream>
using namespace arma;


namespace turb {

    class CMASimple : public SimpleSearcher {
    private:
        int N, lambda, mu;
        double sigma = 0.3,
               stop_fitness = 3e-2,
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

        Searcher* clone() const
        {
            std::cout << "cloned " << check_filename << std::endl;
            Searcher *cloned = new CMASimple(*this);
            cloned->check_filename = check_filename;
            return cloned;
        }
    };
}

#endif

