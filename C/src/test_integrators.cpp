#include "test_integrators.hpp"
#include <fstream>
#include <iostream>
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

        void TestDecayIntegrator::initialize_function(double x, double *result)
        {
            result[0] = 2.0 * cos(x) + 0.03 * cos(11 * x / 12.0);
            result[1] = 0.0;
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

        void TestStabilityIntegrator::initialize_function(double x, double *result)
        {
            result[0] = 2.0 * cos(x) + 0.03 * cos(11 * x / 12.0);
            result[1] = 0.0;
        }

        /* *****        TestNonLiner            ***** */
        void TestNonLinear::serialize(std::ofstream *output, double current_time)
        {
            double sum_u=0, sum_v=0;
            for (size_t i = 0; i < 15; ++i) {
                double *temp_u = c_u[i];
                double *temp_v = c_v[i];
                std::cout << i << " " << *temp_u << " " << *(temp_u + 1) << "\t\t" 
                    << *temp_v << " " << *(temp_v + 1) << std::endl;

                sum_u += sqrt(pow(*temp_u, 2) + pow(*(temp_u + 1), 2));
                sum_v += sqrt(pow(*temp_v, 2) + pow(*(temp_v + 1), 2));
            }
            std::cout << "\n\n";

            *output << current_time << " " << sum_u << " " << sum_v << std::endl;
        }

        void TestNonLinear::apply_step()
        {
            std::cout << "Applting a step" << std::endl;
            compute_nonlinear();
        }

        void TestNonLinear::initialize_function(double x, double *result)
        {
            result[0] = 2.0 * cos(x) * cos(11 * x / 12.0);
            result[1] = sin(2 * x);
        }

}
