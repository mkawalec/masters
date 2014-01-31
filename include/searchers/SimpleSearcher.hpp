#ifndef turb_SimpleSearcher_h
#define turb_SimpleSearcher_h

#include <fftw3.h>

#include "Searcher.hpp"
#include "Integrator.hpp"

namespace turb {

    class SimpleSearcher : public Searcher {
    protected:
        double *f, *du, *f_val1, *dx,
               *d2_v, *d2_u, *d4_u,
               *dv;
        fftw_complex *d_cu, *d2_cv, *d2_cu,
                     *d4_cu, *d_cv;

        fftw_plan du_c, du_r, d2v_c, d2v_r,
                  d2u_r, d4u_r, dv_c, dv_r;

        Integrator *integrator;

        void F(double *__restrict__ input, double *__restrict__ result);
        void compute_F(double *__restrict__ input, double *__restrict__ result);

    public:
        void init();
        virtual void allocate(Integrator *integrator);
        virtual ~SimpleSearcher();
    };
}

#endif
