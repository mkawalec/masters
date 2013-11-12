#include "Integrator.hpp"
#include "helpers.hpp"

#include <fftw3.h>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

namespace turb {

        Integrator::Integrator(size_t dim_power, double timestep, double domain) : 
            dt(timestep), domain_size(domain)
        {
            // Making the code compatible with the
            // ancient FFTW version on cplab computers
            FILE *f = fopen(".wisdom", "r");
            if (f) {
                fftw_import_wisdom_from_file(f);
                fclose(f);
            }

            size_real = pow(2, dim_power);
            size_complex = size_real / 2 + 1;

            temp_array = (double*) fftw_malloc(2 * sizeof(double));

            u = (double*) fftw_malloc(size_real * sizeof(double));
            v = (double*) fftw_malloc(size_real * sizeof(double));
            du = (double*) fftw_malloc(size_real * sizeof(double));
            c_u = (fftw_complex*) fftw_malloc(size_complex * sizeof(fftw_complex));
            c_v = (fftw_complex*) fftw_malloc(size_complex * sizeof(fftw_complex));
            dc_u = (fftw_complex*) fftw_malloc(size_complex * sizeof(fftw_complex));

            // Linear operators acting on u and v
            Lu = (double*) fftw_malloc(size_complex *
                    sizeof(double) * 2.0 / 3);
            Lv = (double*) fftw_malloc(size_complex * 
                    sizeof(double) * 2.0 / 3);

            // Plans to be applied initially
            i_u = fftw_plan_dft_r2c_1d(size_real, u, c_u,
                    FFTW_PATIENT);
            i_v = fftw_plan_dft_r2c_1d(size_real, v, c_v,
                    FFTW_PATIENT);

            // Final transformations
            e_u = fftw_plan_dft_c2r_1d(size_real, c_u, u,
                    FFTW_PATIENT | FFTW_PRESERVE_INPUT);
            e_v = fftw_plan_dft_c2r_1d(size_real, c_v, v,
                    FFTW_PATIENT | FFTW_PRESERVE_INPUT);

            // Forward plans
            f_u = fftw_plan_dft_c2r_1d(size_real, c_u, u,
                    FFTW_PATIENT);
            f_v = fftw_plan_dft_c2r_1d(size_real, c_v, v,
                    FFTW_PATIENT);
            f_du = fftw_plan_dft_c2r_1d(size_real, dc_u, du,
                    FFTW_PATIENT);

            // Backward plans
            b_u = fftw_plan_dft_r2c_1d(size_real, u, c_u,
                    FFTW_PATIENT);
            b_v = fftw_plan_dft_r2c_1d(size_real, v, c_v,
                    FFTW_PATIENT);

            f = fopen(".wisdom", "w");
            fftw_export_wisdom_to_file(f);
            fclose(f);
        }

        Integrator::~Integrator()
        {
            fftw_free(u); fftw_free(v); fftw_free(du);
            fftw_free(c_u); fftw_free(c_v); fftw_free(dc_u);

            fftw_free(Lu); fftw_free(Lv);
            fftw_free(temp_array);
        }

        /** Computes the linear operators acting
         *  on u and v - precomputing saves a lot of
         *  cycles later on.
         */
        void Integrator::initialize_operators() 
        {
            for (size_t i = 0; i < size_complex * 2.0 / 3; ++i) {
                double k = (double) i * 2 * M_PI / domain_size;
                double Lx, Ly;

                Lx = - pow(k, 4) + 2*pow(k, 2) - (1 - e);
                Ly = - D * pow(k, 2) - 1;

                Lu[i] = (1 + 0.5 * dt * Lx) / (1 - 0.5 * dt * Lx);
                Lv[i] = (1 + 0.5 * dt * Ly) / (1 - 0.5 * dt * Ly);
            }
        }

        void Integrator::initialize_function(double x, double *result)
        {
            result[0] = 2.0 * cos(x) + 0.03 * cos(11 * x / 12.0);
            result[1] = 0.0;
        }

        void Integrator::initialize() 
        {
            /** The outer wrapper for initializing, calls the
             *  initialize_function to get the actual values it sets
             */
            for (size_t i = 0; i < size_real; ++i) {
                initialize_function((double) i / size_real * domain_size, temp_array);
                u[i] = temp_array[0];
                v[i] = temp_array[1];
            }

            fftw_execute(i_u);
            fftw_execute(i_v);

            /** Make padidng the padding and scale the 
             *  transformed arrays by sqrt(N)
             */
            override_initialize();

            //double scale_factor = sqrt(2 * M_PI) / size_real;
            double scale_factor = 1 / sqrt(size_real);
            for (size_t i = 0; i < size_complex; ++i) {
                double *tmp_u = c_u[i];
                double *tmp_v = c_v[i];

                if (i >= size_complex * 2.0 / 3) {
                    *tmp_u = 0.0;
                    *(tmp_u + 1) = 0.0;
                    *tmp_v = 0.0;
                    *(tmp_v + 1) = 0.0;
                } else {
                    *tmp_u *= scale_factor;
                    *(tmp_u + 1) *= scale_factor;
                    *tmp_v *= scale_factor;
                    *(tmp_v + 1) *= scale_factor;
                }
            }

            initialize_operators();
        }

        /** \brief Used to override the default
         *      initialization
         */
        void Integrator::override_initialize()
        {
            boost::random::mt19937 rng;
            timeval tv;
            gettimeofday(&tv, NULL);
            rng.seed(1000000 * tv.tv_sec + tv.tv_usec);

            boost::random::uniform_real_distribution<> random_data(0, 
                    3 * pow(size_real, 0.85));
            for (size_t i = 0; i < size_complex; ++i) {
                double *tmp_u = c_u[i];

                *tmp_u          = random_data(rng);
                *(tmp_u + 1)    = random_data(rng);
            }
        }

        /** Computes the linear transformation, which is
         *  essencially just a multiplication by a precomputed
         *  vector. Nothing to see here, move along...
         */
        void Integrator::compute_linear()
        {
            for (size_t i = 0; i < size_complex * 2.0 / 3; ++i) {
                double *tmp_u = c_u[i];
                double *tmp_v = c_v[i];

                *tmp_u *= *(Lu + i);
                *(tmp_u + 1) *= *(Lu + i);
                *tmp_v *= *(Lv + i);
                *(tmp_v + 1) *= *(Lv + i);
            }
        }

        void Integrator::nonlinear_transform(size_t i, double *result)
        {
            double temp_u = (1 - e) *
                (a * v[i] + b * pow(v[i], 2)) * u[i] +
                u[i] * du[i];

            result[0] = u[i] + dt * temp_u;
            result[1] = v[i] + dt * R * pow(u[i], 2);
        }

        void Integrator::compute_nonlinear()
        {
            // Compute the derivative
            double inv_domain_size = 1 / domain_size;
            for (size_t i = 0; i < size_complex; ++i) {
                double *temp_dc_u = dc_u[i];
                double *temp_c_u = c_u[i];
                double k = (double) i * 2 * M_PI * inv_domain_size;

                *temp_dc_u = - k * *(temp_c_u + 1);
                *(temp_dc_u + 1) = k * *temp_c_u;
            }

            fftw_execute(f_u); fftw_execute(f_v); fftw_execute(f_du);

            // Normalize the transform
            double scale_factor = 1 / sqrt(size_real);
            for (size_t i = 0; i < size_real; ++i) {
                *(u + i) *= scale_factor;
                *(v + i) *= scale_factor;
                *(du + i) *= scale_factor;
            }

            // Multiply the relevant parts according to the nonlinear equation
            for (size_t i = 0; i < size_real; ++i) {
                nonlinear_transform(i, temp_array);
                u[i] = temp_array[0];
                v[i] = temp_array[1];
            }

            fftw_execute(b_u);
            fftw_execute(b_v);

            // Normalize and nullify the padding
            for (size_t i = 0; i < size_complex; ++i) {
                double *u = c_u[i];
                double *v = c_v[i];

                if (i >= size_complex * 2.0 / 3) {
                    *u = 0.0;
                    *(u + 1) = 0.0;
                    *v = 0.0;
                    *(v + 1) = 0.0;
                } else { 
                    *u *= scale_factor;
                    *(u + 1) *= scale_factor;
                    *v *= scale_factor;
                    *(v + 1) *= scale_factor;
                }
            }
        }

        void Integrator::apply_step() 
        {
            compute_nonlinear();
            compute_linear();
        }

        void Integrator::forward_transform()
        {
            fftw_execute(e_u); fftw_execute(e_v);
            normalize(u, size_real); normalize(v, size_real);
        }
}

