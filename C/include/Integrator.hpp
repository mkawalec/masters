#ifndef turb_integrator_h
#define turb_integrator_h

#include <fftw3.h>
#include <fstream>
#include "helpers.hpp"

namespace turb {

    class Integrator {
        /** The complex arrays holding the intermediate
         *  and final results of computations
         */
        protected:
        fftw_plan e_u, e_v, i_u, i_v,
                  f_u, f_v, f_du, b_u, b_v;

        size_t size_real, size_complex;

        void initialize_operators();

        /// Computes one timestep of the nonlinear transform
        void compute_nonlinear();
        void compute_linear();

        virtual void override_initialize() {};

        virtual void initialize_function(double x, double *results);
        /** The overloadable nonlinear transform applied 
         *  in the real (non-Fourier) space
         */
        virtual void nonlinear_transform(size_t i, double *results);

        public:
        Integrator(size_t dim_power, double timestep, double domain=2*M_PI);
        ~Integrator();
        void initialize();

        virtual void apply_step();
        virtual void serialize(std::ofstream *output, double current_time);

        double *u, *v, *du, *Lu, *Lv;
        double dt, domain_size;
        fftw_complex *c_u, *c_v, *dc_u;
    };
}

#endif
