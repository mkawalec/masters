#ifndef turb_Searcher_h
#define turb_Searcher_h

#include "Integrator.hpp"
#include "Jacobian.hpp"
#include <vector>

namespace turb {

    class Integrator;

    class Searcher {
    private:
        typedef double jacobian_type;

        Integrator *integrator;
        double *f, *du, *f_val1, *dx,
               *d2_v, *d2_u, *d4_u,
               *dv;
        Jacobian<jacobian_type> *jacobian;
        fftw_complex *d_cu, *d2_cv, *d2_cu,
                     *d4_cu, *d_cv;

        size_t iterations = 20;
        double threshold = 2e-1;
        double overflow = 1e5;
        double h = 1e-6;

        fftw_plan du_c, du_r, d2v_c, d2v_r,
                  d2u_r, d4u_r, dv_c, dv_r;

        void get_jacobian();
        void F(double *input, double *result);

        // Gaussian elimination, must operate on
        // non-singular matrices, will break badly if
        // a singular matrix is provided
        void gauss(double *f, double *result);

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
