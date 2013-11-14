#ifndef turb_integrator_h
#define turb_integrator_h

#include <fftw3.h>
#include <fstream>
#include "helpers.hpp"

namespace turb {

    /*! \brief The class doing actual computation through
     *      integrating a non-linear PDE
     *
     *  It is NOT thread-safe, don't try to use the same
     *  Integrator instance in multiple threads.
     */
    class Integrator {
        /*! \brief The complex arrays holding the intermediate
         *      and final results of computations
         */
        protected:
        fftw_plan e_u, e_v, i_u, i_v,
                  f_u, f_v, f_du, b_u, b_v;


        void initialize_operators();

        /// Computes one timestep of the nonlinear transform
        void compute_nonlinear();
        void compute_linear();

        virtual void override_initialize();

        virtual void initialize_function(double x, double *results);
        /*! \brief The overloadable nonlinear transform applied 
         *      in the real (non-Fourier) space
         */
        virtual void nonlinear_transform(size_t i, double *results);

        // The much needed temprorary array of two doubles
        double *temp_array;

        public:
        Integrator(size_t dim_power, double timestep, double domain=2*M_PI);
        virtual ~Integrator();
        void initialize();

        virtual void apply_step();
        void forward_transform();

        double *u, *v, *du, *Lu, *Lv;
        double dt, domain_size, e, a, b, D, R;
        fftw_complex *c_u, *c_v, *dc_u;
        size_t size_real, size_complex;
    };
}

#endif
