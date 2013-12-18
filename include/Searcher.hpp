#ifndef turb_Searcher_h
#define turb_Searcher_h

#include "Integrator.hpp"
#include "Jacobian.hpp"
#include <vector>

namespace turb {

    class Searcher {
    private:
        Integrator *integrator;
        double *f, *du, *f_val1, *f_val2, *dx;
        Jacobian<long double> *jacobian;
        fftw_complex *d_cu;

        size_t iterations = 50;
        double threshold = 0.0001;
        double overflow = 1e20;
        double h = 0.0001;

        fftw_plan du_c, du_r;

        void get_jacobian();
        void F(double *input, double *result);

        // Gaussian elimination, must operate on 
        // non-singular matrices, will break badly if
        // a singular matrix is provided
        void gauss(double *f, double *result);

        // An in-place quicksort
        void sort_jacobian(int start, int end);
        int get_prefix(size_t line_index, double *where);

        public:
            Searcher(Integrator *integrator);
            ~Searcher();

            std::vector<double> run();
        };

}

#endif
