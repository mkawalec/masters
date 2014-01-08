#include "Searcher.hpp"
#include "Integrator.hpp"
#include "helpers.hpp"
#include "exceptions.hpp"

#include <vector>
#include <fftw3.h>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <cstdlib>

namespace turb {

    Searcher::Searcher(Integrator *integrator) :
        integrator(integrator)
    {
        srand(1234);
        f = (double*)
            fftw_malloc(sizeof(double) * 2 * integrator->size_real);
        du = (double*)
            fftw_malloc(sizeof(fftw_complex) * integrator->size_real);
        dv = (double*)
            fftw_malloc(sizeof(fftw_complex) * integrator->size_real);
        d_cu = (fftw_complex*)
            fftw_malloc(sizeof(fftw_complex) * integrator->size_complex);
        d_cv = (fftw_complex*)
            fftw_malloc(sizeof(fftw_complex) * integrator->size_complex);

        jacobian = new Jacobian<jacobian_type>(2 * integrator->size_real,
                2 * integrator->size_real + 1);

        d2_cv = (fftw_complex*)
            fftw_malloc(sizeof(fftw_complex) * integrator->size_complex);
        d2_v = (double*)
            fftw_malloc(sizeof(fftw_complex) * integrator->size_real);

        d2_cu = (fftw_complex*)
            fftw_malloc(sizeof(fftw_complex) * integrator->size_complex);
        d2_u = (double*)
            fftw_malloc(sizeof(fftw_complex) * integrator->size_real);

        d4_cu = (fftw_complex*)
            fftw_malloc(sizeof(fftw_complex) * integrator->size_complex);
        d4_u = (double*)
            fftw_malloc(sizeof(fftw_complex) * integrator->size_real);

        f_val1 = (double*)
            fftw_malloc(sizeof(double) * 2 * integrator->size_real);
        f_val2 = (double*)
            fftw_malloc(sizeof(double) * 2 * integrator->size_real);

        dx = (double*)
            fftw_malloc(sizeof(double) * 2 * integrator->size_real);

        fftw_import_wisdom_from_filename(".wisdom");
        du_c = fftw_plan_dft_r2c_1d(integrator->size_real, f, d_cu,
                FFTW_PATIENT);
        du_r = fftw_plan_dft_c2r_1d(integrator->size_real, d_cu, du,
                FFTW_PATIENT);

        dv_c = fftw_plan_dft_r2c_1d(integrator->size_real,
                f + integrator->size_real, d_cv, FFTW_PATIENT);
        dv_r = fftw_plan_dft_c2r_1d(integrator->size_real, d_cv, dv,
                FFTW_PATIENT);

        d2v_c = fftw_plan_dft_r2c_1d(integrator->size_real,
                f + integrator->size_real, d2_cv, FFTW_PATIENT);
        d2v_r = fftw_plan_dft_c2r_1d(integrator->size_real, d2_cv, d2_v,
                FFTW_PATIENT);

        d2u_r = fftw_plan_dft_c2r_1d(integrator->size_real, d2_cu, d2_u,
                FFTW_PATIENT);

        d4u_r = fftw_plan_dft_c2r_1d(integrator->size_real, d4_cu, d4_u,
                FFTW_PATIENT);

        fftw_export_wisdom_to_filename(".wisdom");
    }

    void Searcher::init()
    {
        for (size_t i = 0; i < integrator->size_real; ++i) {
            f[i] = integrator->u[i];
            f[i + integrator->size_real] = integrator->v[i];

            // TODO: Implement the commented parts as a test
            /*double x = i * integrator->domain_size / integrator->size_real;
            f[i] = 0.4 * cos(x);
            f[integrator->size_real + i] = 0.5 * sin(x);*/
        }
    }

    std::vector<double> Searcher::run()
    {
        size_t size = integrator->size_real;
        double norm;

        for (size_t i = 0; i < iterations; ++i) {
            // Compute the value of F at current position
            F(f, f_val1);

            if ((norm = l2_norm(f_val1, 2 * size)) < threshold) {
                std::cerr << "Found! " << i << " " << f[0] <<
                    " " << l2_norm(f + size, size) << std::endl;
                return std::vector<double>(f, f + size);
            } else if (norm > overflow) {
                std::cerr << "THROWING " << norm << " " << i << std::endl;
                throw NoResult();
            }
            std::cout << norm << std::endl;

            get_jacobian();
            try {
                gauss(f_val1, dx);
            } catch (const NoResult &e) {
                std::cerr << "No result, wrong at " << e.what() << std::endl;
                continue;
            }

            for (size_t j = 0; j < 2 * size; ++j)
                f[j] += dx[j];
        }

        std::cerr << "Nothing found! " << l2_norm(f_val1, 2 * size) << std::endl;
        throw NoResult();
    }


    void Searcher::F(double *input, double *result)
    {
        double inv_domain_size = 1 / integrator->domain_size;
        double inv_size = 1 / (double) integrator->size_real;
        double tmp[2], tmp2[2];

        // Execute the plan transforming u/v to complex form
        fftw_execute_dft_r2c(du_c, input, d_cu);
        fftw_execute_dft_r2c(dv_c, input + integrator->size_real, d_cv);
        fftw_execute_dft_r2c(d2v_c, input + integrator->size_real, d2_cv);
        for (size_t i = 0; i < integrator->size_complex; ++i) {
            double *d2u = d2_cu[i],
                   *d4u = d4_cu[i],
                   *du = d_cu[i];

            *d2u = *du;
            *(d2_u + 1) = *(du + 1);
            *d4u = *du;
            *(d4u + 1) = *(du + 1);
        }


        // Computing the derivative
        for (size_t i = 0; i < integrator->size_complex; ++i) {
            double k = (double) i * 2 * M_PI * inv_domain_size;
            double *du = d_cu[i];
            double *dv = d_cv[i];
            double *d2v = d2_cv[i];
            double *d2u = d2_cu[i];
            double *d4u = d4_cu[i];

            if (i < integrator->size_complex * 2.0 / 3) {
                tmp[0] = -k * *(du + 1);
                tmp[1] = k * *du;

                tmp2[0] = -k * *(dv + 1);
                tmp2[1] = k * *dv;

                *du = tmp[0];
                *(du + 1) = tmp[1];

                *dv = tmp2[0];
                *(dv + 1) = tmp2[1];

                *d2v *= -pow(k, 2);
                *(d2v + 1) *= -pow(k, 2);
                *d2u *= -pow(k, 2);
                *(d2u + 1) *= -pow(k, 2);
                *d4u *= pow(k, 4);
                *(d4u + 1) *= pow(k, 4);
            } else {
                *du = 0;
                *(du + 1) = 0;
                *dv = 0;
                *(dv + 1) = 0;
                *d2v = 0;
                *(d2v + 1) = 0;
                *d2u = 0;
                *(d2u + 1) = 0;
                *d4u = 0;
                *(d4u + 1) = 0;
            }
        }

        // Transform the derivative back into the real form
        fftw_execute_dft_c2r(du_r, d_cu, du);
        fftw_execute_dft_c2r(dv_r, d_cv, dv);
        fftw_execute_dft_c2r(d2v_r, d2_cv, d2_v);
        fftw_execute_dft_c2r(d2u_r, d2_cu, d2_u);
        fftw_execute_dft_c2r(d4u_r, d4_cu, d4_u);

        for (int i = 0; (unsigned)i < integrator->size_real; ++i) {
            d2_v[i] *= inv_size;
            d2_u[i] *= inv_size;
            d4_u[i] *= inv_size;
            du[i] *= inv_size;
            dv[i] *= inv_size;
        }

        compute_F(input, result);
    }

    void Searcher::compute_F(double *input, double *result) {
        int size = integrator->size_real;
        for (int i = 0; i < size; ++i) {
            double u = input[i];
            double v = input[i + size];

            result[i] = -d4_u[i] - 2 * d2_u[i] +
                (-1 + integrator->e + (1 - integrator->e) *
                (integrator->a * v + integrator->b * pow(v, 2)) +
                du[i]) * u;

            result[i + size] = -v + integrator->D * d2_v[i] + integrator->R * pow(u, 2);
            /*result[i] = d4_u[i] + u * dv[i] - 2 * pow(u, 2) - (1 - u) * u;
            result[i + size] = d2_v[i] + v;*/
        }

    }

    void Searcher::gauss(double *f, double *result)
    {
        int size = jacobian->dims().first;

        // Add F as a suffix to jacobian matrix
        for (int i = 0; i < size; ++i)
            (*jacobian)[i][size] = -f[i];

        // Bring Jacobian to the echelon form
        for (int i = 0; i < size; ++i) {
            sort_jacobian(i, jacobian->dims().first);

            int major_nonzero_i = (*jacobian)[i].prefix();
            for (int j = i + 1; j < size; ++j) {
                int nonzero_i = (*jacobian)[j].prefix();
                if (nonzero_i != major_nonzero_i) break;

                jacobian_type factor = (*jacobian)[j][nonzero_i] /
                                (*jacobian)[i][nonzero_i];

                for (int k = nonzero_i; k < size + 1; ++k)
                    (*jacobian)[j][k] -= factor * (*jacobian)[i][k];
            }
        }

        // Check if Jacobian is upper-triangular
        for (int i = 0; i < size; ++i) {
            if ((*jacobian)[i].prefix() != i) throw NoResult(i);
        }

        // Compute the result
        for (int i = size - 1; i >= 0; --i) {
            result[i] = (*jacobian)[i][size];

            for (int j = size - 1; j >= i; j--) {
                if (j != i)
                    result[i] -= result[j] * (*jacobian)[i][j];
                else
                    result[i] /= (*jacobian)[i][i];
            }
        }
    }

    // Make the partition in-place
    void Searcher::sort_jacobian(int start, int end)
    {
        if (end - start < 2) return;

        int anchor = start + rand()%(end-start);
        int anchor_prefix = (*jacobian)[anchor].prefix();
        int store_index = start;

        // In-place partition
        jacobian->swap_lines(anchor, end - 1);
        for (int i = start; i < end - 1; ++i) {
            if ((*jacobian)[i].prefix() < anchor_prefix) {
                jacobian->swap_lines(i, store_index);
                ++store_index;
            }
        }
        jacobian->swap_lines(store_index, end - 1);

        // Launch the sorts
        sort_jacobian(start, store_index);
        sort_jacobian(store_index + 1, end);
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

            for (size_t j = 0; j < size; ++j)
                (*jacobian)[j][i] =
                    (f_val1[j] - f_val2[j]) * inv_h;
        }
        fftw_free(tmp_f);
    }

    Searcher::~Searcher()
    {
        fftw_free(f);
        fftw_free(du);
        fftw_free(d2_v);
        fftw_free(d2_u);
        fftw_free(d4_u);

        fftw_free(d_cu);
        fftw_free(d2_cv);
        fftw_free(d2_cu);
        fftw_free(d4_cu);

        fftw_free(f_val1);
        fftw_free(f_val2);
        fftw_free(dx);

        delete jacobian;
    }
}
