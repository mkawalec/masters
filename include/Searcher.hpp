#ifndef turb_Searcher_h
#define turb_Searcher_h

#include "Integrator.hpp"
#include <vector>

namespace turb {

    class Searcher {
    private:
        Integrator *integrator;
        double *f, *du, *jacobian, 
               *f_val1, *f_val2, *dx;
        fftw_complex *d_cu;

        size_t iterations = 50;
        double threshold = 0.0001;
        double overflow = 1e20;
        double h = 0.0001;

        fftw_plan du_c, du_r;

        void get_jacobian();
        void F(double *input, double *result);

        // TODO: Actually implement the Gauss' elimination
        void gauss(double *f, double *result);

    public:
        Searcher(Integrator *integrator);
        ~Searcher();

        std::vector<double> run();
    };
}

#endif
