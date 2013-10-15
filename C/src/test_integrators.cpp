#include "test_integrators.hpp"
#include <fstream>
#include <fftw3.h>


namespace turb {

        /* *****      TestDecayIntegrator      ***** */
        void TestDecayIntegrator::serialize(std::ofstream *output, double current_time)
        {
            *output << current_time << " " << *(c_u[0]) << std::endl;
        }

        void TestDecayIntegrator::apply_step() 
        {
            compute_linear();
        }

        void TestDecayIntegrator::override_initialize()
        {
            fftw_execute(e_u); fftw_execute(e_v);
            normalize(u, size_real); normalize(v, size_real);

            for (size_t i = 0; i < size_complex; ++i) {
                double *tmp_u = c_u[i];
                *tmp_u = i + 1;
                *(tmp_u + 1) = 0;
            }
        }

        /* *****        TestMultIntegrator      ***** */

        void TestMultIntegrator::serialize(std::ofstream *output, double current_time)
        {
            fftw_execute(e_u); fftw_execute(e_v);
            normalize(u, size_real); normalize(v, size_real);
            for (size_t i = 0; i < size_real; ++i) 
                *output << i << " " << u[i] << " " << v[i] << std::endl;

            *output << std::endl;
        }

        void TestMultIntegrator::apply_step() 
        {
            compute_nonlinear();
        }

        void TestMultIntegrator::initialize_function(double x, double *results)
        {
            results[0] = cos(x);
            results[1] = cos(2 * x);
        }


        void TestMultIntegrator::nonlinear_transform(size_t i, double *results)
        {
            results[0] = u[i] * v[i];
            results[1] = v[i];
        }

        /* *****        TestStabilityIntegrator     ***** */
        void TestStabilityIntegrator::apply_step()
        {
            compute_nonlinear();
        }

        void TestStabilityIntegrator::nonlinear_transform(size_t i, double *results)
        {
            results[0] = u[i];
            results[1] = v[i];
        }
}
