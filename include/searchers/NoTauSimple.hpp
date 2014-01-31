#ifndef turb_NoTauSimple_h
#define turb_NoTauSimple_h

#include "Searcher.hpp"
#include "Integrator.hpp"
#include "Jacobian.hpp"
#include "searchers/SimpleSearcher.hpp"

namespace turb {

    class NoTauSimple : public SimpleSearcher {
    private:
        typedef double jacobian_type;

        Jacobian<jacobian_type> *jacobian;

        size_t iterations = 20;
        double threshold = 2e-1;
        double overflow = 1e5;
        double h = 1e-6;

        void get_jacobian();

        // Gaussian elimination, must operate on
        // non-singular matrices, will break badly if
        // a singular matrix is provided
        void gauss(double *f, double *result);

    public:
        NoTauSimple();
        virtual ~NoTauSimple();

        std::vector<double> run();
        void allocate(Integrator *integrator);

        Searcher* clone() const { return new NoTauSimple(*this);}
    };
}

#endif
