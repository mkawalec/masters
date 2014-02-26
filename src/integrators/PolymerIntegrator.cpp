#include "integrators/PolymerIntegrator.hpp"
#include "helpers.hpp"

#include <mpi.h>
#include <fftw3.h>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include <boost/program_options.hpp>
#include <boost/date_time.hpp>
namespace po = boost::program_options;
namespace pt = boost::posix_time;

namespace turb {

    PolymerIntegrator::PolymerIntegrator()
    {
        // Set the info for runtime class discovery
        name = "polymer";
        class_name = "integrator";
        description = " integrator trying to model polymers.";
        set_options();
        Integrator::available.push_back(this);

        search = false;
    }

    void PolymerIntegrator::set_options()
    {
        options.reset(new po::options_description(name + " options"));
        options->add_options()
            ("e,e", po::value<double>(&e)->default_value(-0.1),
             "an integration parameter")
            ("a,a", po::value<double>(&a)->default_value(0.125),
             "an integration parameter")
            ("b,b", po::value<double>(&b)->default_value(-0.004),
             "an integration parameter")
            ("D,D", po::value<double>(&D)->default_value(40),
             "an integration parameter")
            ("R,R", po::value<double>(&R)->default_value(1.04),
             "integration parameter corresponding to a Reynolds number")
            ("lambda", po::value<double>(&lambda)->default_value(1),
             "elasticity constant")
            ;
    }

    void PolymerIntegrator::clear(size_t dim_power, double dt, double domain_size)
    {
        fftw_free(u); fftw_free(v); fftw_free(du);
        fftw_free(tau); fftw_free(dtau);

        fftw_free(c_u); fftw_free(c_v); fftw_free(dc_u);
        fftw_free(c_tau); fftw_free(c_dtau);

        fftw_free(Lv);
        fftw_free(temp_array);

        delete L_u_tau;

        initialize(dim_power, dt, domain_size);
    }


    PolymerIntegrator::~PolymerIntegrator()
    {
        fftw_free(u); fftw_free(v); fftw_free(du);
        fftw_free(tau); fftw_free(dtau);

        fftw_free(c_u); fftw_free(c_v); fftw_free(dc_u);
        fftw_free(c_tau); fftw_free(c_dtau);

        fftw_free(Lv);
        delete L_u_tau;
        fftw_free(temp_array);
    }

    void PolymerIntegrator::allocate(size_t dim_power, double timestep, double domain)
    {
        dt = timestep;
        domain_size = domain;
        dim_power = dim_power;

        fftw_import_wisdom_from_filename(get_wisdom_filename().c_str());

        size_real = pow(2, dim_power);
        size_complex = size_real / 2 + 1;

        temp_array = (double*) fftw_malloc(3 * sizeof(double));

        u = (double*) fftw_malloc(size_real * sizeof(double));
        v = (double*) fftw_malloc(size_real * sizeof(double));
        tau = (double*) fftw_malloc(size_real * sizeof(double));

        du = (double*) fftw_malloc(size_real * sizeof(double));
        dtau = (double*) fftw_malloc(size_real * sizeof(double));

        c_u = (fftw_complex*) fftw_malloc(size_complex * sizeof(fftw_complex));
        c_v = (fftw_complex*) fftw_malloc(size_complex * sizeof(fftw_complex));
        dc_u = (fftw_complex*) fftw_malloc(size_complex * sizeof(fftw_complex));
        c_tau = (fftw_complex*) fftw_malloc(size_complex * sizeof(fftw_complex));
        c_dtau = (fftw_complex*) fftw_malloc(size_complex * sizeof(fftw_complex));

        // Linear operators acting on u and v
        Lv = (double*) fftw_malloc(size_complex *
                sizeof(double) * 2.0 / 3);
        L_u_tau = new matrix<std::complex<double> >(2 * size_complex, 2 * size_complex);

        // Plans to be applied initially
        i_u = fftw_plan_dft_r2c_1d(size_real, u, c_u,
                FFTW_MEASURE);
        i_v = fftw_plan_dft_r2c_1d(size_real, v, c_v,
                FFTW_MEASURE);
        i_tau = fftw_plan_dft_r2c_1d(size_real, tau, c_tau,
                FFTW_MEASURE);

        // Final transformations
        e_u = fftw_plan_dft_c2r_1d(size_real, c_u, u,
                FFTW_MEASURE | FFTW_PRESERVE_INPUT);
        e_v = fftw_plan_dft_c2r_1d(size_real, c_v, v,
                FFTW_MEASURE | FFTW_PRESERVE_INPUT);
        e_tau = fftw_plan_dft_c2r_1d(size_real, c_tau, tau,
                FFTW_MEASURE | FFTW_PRESERVE_INPUT);

        // Forward plans
        f_u = fftw_plan_dft_c2r_1d(size_real, c_u, u,
                FFTW_MEASURE);
        f_v = fftw_plan_dft_c2r_1d(size_real, c_v, v,
                FFTW_MEASURE);
        f_tau = fftw_plan_dft_c2r_1d(size_real, c_tau, tau,
                FFTW_MEASURE);
        f_du = fftw_plan_dft_c2r_1d(size_real, dc_u, du,
                FFTW_MEASURE);
        f_dtau = fftw_plan_dft_c2r_1d(size_real, c_dtau, dtau,
                FFTW_MEASURE);

        // Backward plans
        b_u = fftw_plan_dft_r2c_1d(size_real, u, c_u,
                FFTW_MEASURE);
        b_v = fftw_plan_dft_r2c_1d(size_real, v, c_v,
                FFTW_MEASURE);
        b_tau = fftw_plan_dft_r2c_1d(size_real, tau, c_tau,
                FFTW_MEASURE);

        fftw_export_wisdom_to_filename(get_wisdom_filename().c_str());
    }

    /** Computes the linear operators acting
     *  on u and v - precomputing saves a lot of
     *  cycles later on.
     */
    void PolymerIntegrator::initialize_operators()
    {
        for (size_t i = 0; i < size_complex * 2.0 / 3; ++i) {
            double k = (double) i * 2 * M_PI / domain_size;
            double Ly;

            Ly = - D * pow(k, 2) - 1;

            Lv[i] = (1 + 0.5 * dt * Ly) / (1 - 0.5 * dt * Ly);

            // Set the L_u_tau operator
        }
    }

    void PolymerIntegrator::initialize_function(double x, double *result)
    {
        result[0] = 2.0 * cos(x) + 0.03 * cos(11 * x / 12.0);
        result[1] = 0.0;
        result[2] = 0.0;
    }

    void PolymerIntegrator::initialize(size_t dim_power, double timestep, double domain)
    {
        allocate(dim_power, timestep, domain);

        /** The outer wrapper for initializing, calls the
         *  initialize_function to get the actual values it sets
         */
        for (size_t i = 0; i < size_real; ++i) {
            initialize_function((double) i / size_real * domain_size, temp_array);
            u[i] = temp_array[0];
            v[i] = temp_array[1];
            tau[i] = temp_array[2];
        }

        fftw_execute(i_u);
        fftw_execute(i_v);
        fftw_execute(i_tau);

        /** Make padidng the padding and scale the
         *  transformed arrays by sqrt(N)
         */
        override_initialize();

        double scale_factor = 1 / sqrt(size_real);
        for (size_t i = 0; i < size_complex; ++i) {
            double *tmp_u = c_u[i];
            double *tmp_v = c_v[i];
            double *tmp_tau = c_tau[i];

            if (i >= size_complex * 2.0 / 3) {
                *tmp_u = 0.0;
                *(tmp_u + 1) = 0.0;
                *tmp_v = 0.0;
                *(tmp_v + 1) = 0.0;
                *tmp_tau = 0.0;
                *(tmp_tau + 1) = 0.0;
            } else {
                *tmp_u *= scale_factor;
                *(tmp_u + 1) *= scale_factor;
                *tmp_v *= scale_factor;
                *(tmp_v + 1) *= scale_factor;
                *tmp_tau *= scale_factor;
                *(tmp_tau + 1) *= scale_factor;
            }
        }

        initialize_operators();
    }

    /** \brief Used to override the default
     *      initialization
     */
    void PolymerIntegrator::override_initialize()
    {
        pt::ptime t = pt::microsec_clock::local_time();
        pt::ptime tick = pt::second_clock::local_time();
        pt::time_duration diff = t - tick;

        std::mt19937 generator(diff.total_nanoseconds());
        double divisor = 3 * pow(size_real, 0.85) / generator.max();

        for (size_t i = 0; i < size_complex; ++i) {
            double *tmp_u = c_u[i];

            *tmp_u          = generator() * divisor;
            *(tmp_u + 1)    = generator() * divisor;
        }
    }

    /** Computes the linear transformation, which is
     *  essentially just a multiplication by a precomputed
     *  vector. Nothing to see here, move along...
     */
    void PolymerIntegrator::compute_linear()
    {
        for (size_t i = 0; i < size_complex * 2.0 / 3; ++i) {
            double *tmp_u = c_u[i];
            double *tmp_v = c_v[i];
            double *tmp_tau = c_tau[i];

            *tmp_v *= *(Lv + i);
            *(tmp_v + 1) *= *(Lv + i);
        }
    }

    void PolymerIntegrator::nonlinear_transform(size_t i, double *result)
    {
        double temp_u = (1 - e) *
            (a * v[i] + b * pow(v[i], 2)) * u[i] +
            u[i] * du[i] - dtau[i];

        double temp_tau =
            1 / lambda * du[i] - u[i] * dtau[i] + 2 * tau[i] * du[i];

        result[0] = u[i] + dt * temp_u;
        result[1] = v[i] + dt * R * pow(u[i], 2);
        result[2] = tau[i] + dt * temp_tau;
    }

    void PolymerIntegrator::compute_nonlinear()
    {
        // Compute the derivative
        double inv_domain_size = 1 / domain_size;
        for (size_t i = 0; i < size_complex; ++i) {
            double *temp_dc_u = dc_u[i],
                   *temp_c_u = c_u[i],
                   *temp_dc_tau = c_dtau[i],
                   *temp_c_tau = c_tau[i];

            double k = (double) i * 2 * M_PI * inv_domain_size;

            *temp_dc_u = - k * *(temp_c_u + 1);
            *(temp_dc_u + 1) = k * *temp_c_u;

            *temp_dc_tau = -k * *(temp_c_tau + 1);
            *(temp_dc_tau + 1) = k * *temp_c_tau;
        }

        fftw_execute(f_u); fftw_execute(f_v); fftw_execute(f_du);
        fftw_execute(f_tau); fftw_execute(f_dtau);

        // Normalize the transform
        double scale_factor = 1 / sqrt(size_real);
        for (size_t i = 0; i < size_real; ++i) {
            *(u + i) *= scale_factor;
            *(v + i) *= scale_factor;
            *(du + i) *= scale_factor;
            *(tau + i) *= scale_factor;
            *(dtau + i) *= scale_factor;
        }

        // Multiply the relevant parts according to the nonlinear equation
        for (size_t i = 0; i < size_real; ++i) {
            nonlinear_transform(i, temp_array);
            u[i] = temp_array[0];
            v[i] = temp_array[1];
            tau[i] = temp_array[2];
        }

        fftw_execute(b_u);
        fftw_execute(b_v);
        fftw_execute(b_tau);


        // Normalize and nullify the padding
        for (size_t i = 0; i < size_complex; ++i) {
            double *u = c_u[i];
            double *v = c_v[i];
            double *tau = c_tau[i];

            if (i >= size_complex * 2.0 / 3) {
                *u = 0.0;
                *(u + 1) = 0.0;
                *v = 0.0;
                *(v + 1) = 0.0;
                *tau = 0.0;
                *(tau + 1) = 0.0;
            } else {
                *u *= scale_factor;
                *(u + 1) *= scale_factor;
                *v *= scale_factor;
                *(v + 1) *= scale_factor;
                *tau *= scale_factor;
                *(tau + 1) *= scale_factor;
            }
        }
    }

    void PolymerIntegrator::forward_transform()
    {
        fftw_execute(e_u);
        fftw_execute(e_v);
        fftw_execute(e_tau);

        normalize(u, size_real);
        normalize(v, size_real);
        normalize(tau, size_real);
    }

    PolymerIntegrator *polymer_integrator_instance = new PolymerIntegrator();
}

