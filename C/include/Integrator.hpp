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
        Integrator(size_t dim_power, double timestep);
        ~Integrator();
        void initialize();

        virtual void apply_step();
        virtual void serialize(std::ofstream *output, double current_time);

        double *u, *v, *du, *Lu, *Lv;
        double dt;
        fftw_complex *c_u, *c_v, *dc_u;
    };

    /** Enables checking of the decay rate on 0-th 
     *  component of u (which should be 2 * (1 + e))
     */
    class TestDecayIntegrator : public Integrator {
        protected:
        void override_initialize();

        public:
        TestDecayIntegrator(size_t dim_power, double timestep) : Integrator(dim_power, timestep) {};

        void serialize(std::ofstream *output, double current_time);
        void apply_step();
    };

    class TestMultIntegrator : public Integrator {
        protected:
        void nonlinear_transform(size_t i, double *results);
        void initialize_function(double x, double *results);

        public:
        TestMultIntegrator(size_t dim_power, double timestep) : Integrator(dim_power, timestep) {};

        void serialize(std::ofstream *output, double current_time);
        void apply_step();
    };

    class TestStabilityIntegrator : public Integrator {
        protected:
        void nonlinear_transform(size_t i, double *results);

        public:
        TestStabilityIntegrator(size_t dim_power, double timestep) : Integrator(dim_power, timestep) {};

        void apply_step();
    };
        
}

#endif
