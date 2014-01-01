#ifndef turb_Searcher_h
#define turb_Searcher_h

#include "Integrator.hpp"
#include "Jacobian.hpp"
#include <vector>

namespace turb {

    class Searcher {
    private:
        typedef __float128 jacobian_type;

        Integrator *integrator;
        double *f, *du, *f_val1, *f_val2, *dx,
               *d2_v, *d2_u, *d4_u;
        Jacobian<jacobian_type> *jacobian;
        fftw_complex *d_cu, *d2_cv, *d2_cu,
                     *d4_cu;

        size_t iterations = 20;
        double threshold = 1e-4;
        double overflow = 1e10;
        double h = 0.0001;

        fftw_plan du_c, du_r, d2v_c, d2v_r,
                  d2u_c, d2u_r, d4u_c, d4u_r;

        void get_jacobian();
        void F(double *input, double *result);

        // Gaussian elimination, must operate on
        // non-singular matrices, will break badly if
        // a singular matrix is provided
        void gauss(double *f, double *result);

        // An in-place quicksort
        void sort_jacobian(int start, int end);
        int get_prefix(size_t line_index, double *where);

    protected:
        virtual void compute_F(double *input, double *result);

    public:
        Searcher(Integrator *integrator);
        virtual ~Searcher();

        std::vector<double> run();
        void init();
    };

}

#endif
