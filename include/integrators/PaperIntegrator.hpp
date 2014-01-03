#ifndef turb_PaperIntegrator_h
#define turb_PaperIntegrator_h

#include "helpers.hpp"
#include "Integrator.hpp"

#include <fftw3.h>
#include <fstream>

namespace turb {

    /*! \brief The class doing actual computation through
     *      integrating a non-linear PDE
     *
     *  It is NOT thread-safe, don't try to use the same
     *  Integrator instance in multiple threads.
     */
    class PaperIntegrator : Integrator {
        /*! \brief The complex arrays holding the intermediate
         *      and final results of computations
         */
        protected:
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
        PaperIntegrator(size_t dim_power, double timestep, double domain=2*M_PI);
        virtual ~PaperIntegrator();
        void initialize();

        void forward_transform();

        double *u, *v, *du, *Lu, *Lv;
        double e=-0.1, a=0.125, b=-0.004, D=40, R=1.04;

        fftw_complex *c_u, *c_v, *dc_u;
    };
}

#endif
