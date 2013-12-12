#ifndef turb_test_integrators_h
#define turb_test_integrators_h

#include "Integrator.hpp"
#include <fstream>

namespace turb {

    /*! \brief Enables checking of the decay rate on 0-th 
     *      component of u (which should be 2 * (1 + e))
     */
    class TestDecayIntegrator : public Integrator {
        protected:
        void override_initialize();
        void initialize_function(double x, double *results);

        public:
        TestDecayIntegrator(size_t dim_power, double timestep) : Integrator(dim_power, timestep) {};

        void serialize(std::ofstream *output, double current_time);
        void apply_step();
    };

    /*! \brief Checks if function multiplication behaves
     *      correctly under Fourier transform
     */
    class TestMultIntegrator : public Integrator {
        protected:
        void nonlinear_transform(size_t i, double *results);
        void initialize_function(double x, double *results);

        public:
        TestMultIntegrator(size_t dim_power, double timestep) : Integrator(dim_power, timestep) {};

        void serialize(std::ofstream *output);
        void apply_step();
        void override_initialize() {};
    };

    /*! \brief Checks if functions stay stable under Fourier
     *      transform.
     */
    class TestStabilityIntegrator : public Integrator {
        protected:
        void nonlinear_transform(size_t i, double *results);
        void initialize_function(double x, double *results);

        public:
        TestStabilityIntegrator(size_t dim_power, double timestep) : Integrator(dim_power, timestep) {};

        void apply_step();
    };

    /*! \brief Checks if the non-linear part of the integrator
     *      behaves correctly.
     */
    class TestNonLinear : public Integrator {
        protected:
        void initialize_function(double x, double *results);

        public:
        TestNonLinear(size_t dim_power, double timestep, double domain) : Integrator(dim_power, timestep, domain) {};
        void serialize(std::ofstream *output, double current_time);
        void apply_step();
        void override_initialize() {};
    };

}


#endif
