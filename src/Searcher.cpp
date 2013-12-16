#include "Searcher.hpp"
#include "Integrator.hpp"
#include "helpers.hpp"
#include "exceptions.hpp"

#include <vector>
#include <fftw3.h>
#include <cmath>
#include <iostream>
#include <cstdio>

namespace turb {

    Searcher::Searcher(Integrator *integrator) :
        integrator(integrator)
    {
        f = (double*) 
            fftw_malloc(sizeof(double) * 2 * integrator->size_real);
        du = (double*)
            fftw_malloc(sizeof(fftw_complex) * integrator->size_real);
        d_cu = (fftw_complex*)
            fftw_malloc(sizeof(fftw_complex) * integrator->size_complex);
        jacobian = (double*)
            fftw_malloc(sizeof(double) * pow(2 * integrator->size_real, 2));

        f_val1 = (double*)
            fftw_malloc(sizeof(double) * 2 * integrator->size_real);
        f_val2 = (double*)
            fftw_malloc(sizeof(double) * 2 * integrator->size_real);

        dx = (double*)
            fftw_malloc(sizeof(double) * 2 * integrator->size_real);

        for (size_t i = 0; i < 2 * integrator->size_real; ++i) {
            if (i < integrator->size_real)
                f[i] = integrator->u[i];
            else 
                f[i] = integrator->v[i-integrator->size_real];
        }

        fftw_import_wisdom_from_filename(".wisdom");
        du_c = fftw_plan_dft_r2c_1d(integrator->size_real, f, d_cu,
                FFTW_PATIENT);
        du_r = fftw_plan_dft_c2r_1d(integrator->size_real, d_cu, du,
                FFTW_PATIENT);

        fftw_export_wisdom_to_filename(".wisdom");
    }

    void Searcher::F(double *input, double *result)
    {
        double inv_domain_size = 1 / integrator->domain_size;
        double inv_size = 1 / (double) integrator->size_real;
        double tmp[2];

        // Execute the plan transforming u to complex form
        fftw_execute(du_c);

        // Computing the derivative
        for (size_t i = 0; i < integrator->size_complex; ++i) {
            double k = (double) i * 2 * M_PI * inv_domain_size;
            double *inp = d_cu[i];
            tmp[0] = -k * *(inp + 1);
            tmp[1] = k * *inp;

            *inp = tmp[0] * inv_size;
            *(inp + 1) = tmp[1] * inv_size;
        }

        // Transform the derivative back into the real form
        fftw_execute(du_r);

        size_t size = integrator->size_real;
        for (size_t i = 0; i < integrator->size_real; ++i) {
            double k = (double) i * 2 * M_PI * inv_domain_size;
            double u = input[i];
            double v = input[i + size];

            double u_mult = -pow(k, 4) + 2 * pow(k, 2) - 
                1 + integrator->e + (1 - integrator->e) * 
                (integrator->a * v + integrator->b * pow(v, 2)) +
                du[i];
            double v_mult = -integrator->D * pow(k, 2) - 1;

            result[i] = u * u_mult;
            result[i + size] = v * v_mult + integrator->R * pow(u, 2);
        }
    }

    void Searcher::jacobian_mult(double *what, double *result)
    {
        size_t size = 2 * integrator->size_real;

        for (size_t j = 0; j < size; ++j) {
            result[j] = 0;

            for (size_t k = 0; k < size; ++k) 
                result[j] += jacobian[j + size * k] * what[k];
        }
    }

    std::vector<double> Searcher::run()
    {
        //std::cout << "in!"<< std::endl;
        size_t size = 2 * integrator->size_real;
        double *tmp_f = (double*) fftw_malloc(sizeof(double) * size);
        double norm;

        for (size_t i = 0; i < iterations; ++i) {
            get_jacobian();

            // Get the gradient
            jacobian_mult(f, dx);
            //std::cout << " " << l2_norm(f, size) << " " << l2_norm(dx, size) << std::endl;

            // Compute the value of F at current position
            F(f, f_val1);
            //std::cout << l2_norm(f_val1, size) << " " << l2_norm(f, size) << std::endl;
            jacobian_mult(dx, f_val2);

            for (size_t j = 0; j < size; ++j)
                f_val1[j] += f_val2[j];
            if ((norm = l2_norm(f_val1, size)) < threshold)
                return std::vector<double>(f_val1, f_val1 + size);
            else if (norm > overflow)
                throw NoResult();

            // Set a new value of f
            double *tmp = f;
            f = dx;
            dx = tmp;
        }

        std::cout << "nothing!" << std::endl;
        fftw_free(tmp_f);
        throw NoResult();
    }

    void Searcher::get_jacobian()
    {
        size_t size = 2 * integrator->size_real;
        double *tmp_f = (double*) fftw_malloc(sizeof(double) * size);
        double inv_h = 1 / h;

        for (size_t i = 0; i < size; ++i) {
            // Initialize the tmp array
            for (size_t j = 0; j < size; ++j)
                tmp_f[j] = f[j];
            tmp_f[i] += h;

            // Calculate the values of F
            F(tmp_f, f_val1);
            F(f, f_val2);
            //for (size_t j = 0; j < size; j++) 
            //    std::cout << f_val1[j] << " ";

            // Subtract and invert
            for (size_t j = 0; j < size; ++j) {
                jacobian[size * i + j] = 
                    (f_val1[j] - f_val2[j]) * inv_h;
                //std::cout << jacobian[size * i + j] << " ";
            }
            //std::cout << std::endl;

        }
        //system("sleep 30s");
        fftw_free(tmp_f);
    }



    Searcher::~Searcher()
    {
        fftw_free(f);
        fftw_free(du);
        fftw_free(d_cu);
        fftw_free(jacobian);
        fftw_free(f_val1);
        fftw_free(f_val2);
        fftw_free(dx);
    }
}
