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
        f = (double*) 
            fftw_malloc(sizeof(double) * 2 * integrator->size_real);
        du = (double*)
            fftw_malloc(sizeof(fftw_complex) * integrator->size_real);
        d_cu = (fftw_complex*)
            fftw_malloc(sizeof(fftw_complex) * integrator->size_complex);
        jacobian = (double*)
            fftw_malloc(sizeof(double) * 
                    2 * integrator->size_real * (2 * integrator->size_real + 1));

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
            double u = input[i];
            double v = input[i + size];
            double k = (double) i * 2 * M_PI * inv_domain_size;

            double u_mult = -pow(k, 4) + 2 * pow(k, 2) - 
                1 + integrator->e + (1 - integrator->e) * 
                (integrator->a * v + integrator->b * pow(v, 2)) +
                du[i];
            double v_mult = -integrator->D * pow(k, 2) - 1;

            result[i] = u * u_mult;
            result[i + size] = v * v_mult + integrator->R * pow(u, 2);
        }
    }

    std::vector<double> Searcher::run()
    {
        size_t size = 2 * integrator->size_real;
        double *tmp_f = (double*) fftw_malloc(sizeof(double) * size);
        double norm;

        for (size_t i = 0; i < iterations; ++i) {
            get_jacobian();

            // Compute the value of F at current position
            F(f, f_val1);
            try {
                gauss(f_val1, dx);
            } catch (const NoResult &e) {
                std::cout << "Caught " << e.what() << std::endl;
                continue;
            }

            for (int j = 0; j < size; ++j)
                std::cout << dx[j] << " ";
            std::cout << std::endl;

            for (size_t j = 0; j < size; ++j)
                f[j] += dx[j];

            if ((norm = l2_norm(f_val1, size)) < threshold) {
                std::cout << "Found! " << norm << std::endl;
                fftw_free(tmp_f);
                return std::vector<double>(f_val1, f_val1 + size);
            } else if (norm > overflow) {
                fftw_free(tmp_f);
                throw NoResult();
            }
        }

        //std::cout << "nothing!" << std::endl;
        fftw_free(tmp_f);
        throw NoResult();
    }

    void Searcher::gauss(double *f, double *result)
    {
        size_t size = 2 * integrator->size_real;
        srand(1234);

        // Add F as a suffix to jacobian matrix
        for (size_t i = 0; i < size; ++i) 
            jacobian[(size + 1) * i + size] = -f[i];

        // Bring jacobian to the echelon form
        for (size_t i = 0; i < size; ++i) {
            sort_jacobian(&jacobian[(size+1)*i], size-i);

            size_t major_nonzero_i = get_prefix(i, jacobian);
            for (size_t j = i + 1; j < size; ++j) {
                size_t nonzero_i = get_prefix(j, jacobian);
                if (nonzero_i != major_nonzero_i) break;

                double factor = jacobian[(size+1) * j + nonzero_i] / 
                                jacobian[(size+1) * i + nonzero_i];

                for (size_t k = nonzero_i; k < size + 1; ++k) 
                    jacobian[(size+1) * j + k] -= 
                        factor * jacobian[(size+1) * i + k];
            }
        }

        for (int i = 0; i < size; ++i)
            //std::cout << jacobian[(size + 1) * i + size - 3] << " ";
            std::cout << i << " " << get_prefix(i, jacobian) << " ";
        std::cout << std::endl;
        /*for (int i = 0; i < size; ++i)
            std::cout << jacobian[(size+1) * i + size] << " ";

        std::cout << std::endl;
        std::cout << size << std::endl;*/

        // Compute the result
        for (int i = size - 1; i >= 0; --i) {
            result[i] = jacobian[(size + 1) * i + size];

            for (int j = size; j >= i; j--) {
                size_t current_idx = (size + 1) * i + j;

                if (j != i) result[i] -= result[j] * jacobian[current_idx];
                else {
                    result[i] /= jacobian[(size+1) * i + i];
                    std::cout << jacobian[(size+1) * i + i] << " ";
                }
            }
            //std::cout << result[i] << " " << fuzzy_eql(result[i], 0) << " ";
            if (fuzzy_eql(result[i], 0) && i != 0) throw NoResult(i);
        }
    }

    int Searcher::get_prefix(size_t line_index, double *where)
    {
        int zero_count = 0;
        int line_size = 2 * integrator->size_real + 1;

        for (int i = 0; i < 2 * integrator->size_real; ++i) {
            if (!fuzzy_eql(where[line_size * line_index + i], 0)) break;
            ++zero_count;
        }
        return zero_count;
    }

    // Make the partition in-place
    void Searcher::sort_jacobian(double *where, int size)
    {
        if (size < 2) return;

        int anchor = rand()%size;
        int smaller_count=0, bigger_count=0;
        int anchor_prefix = get_prefix(anchor, where);
        int line_size = 2 * integrator->size_real + 1;

        // Find out how many bigger and smaller indexes there are
        for (int i=0; i < size; ++i) {
            if (get_prefix(i, where) > anchor_prefix)
                ++bigger_count;
            else if (i != anchor)
                ++smaller_count;
        }

        // Allocate the smaller and bigger arrays
        double *smaller = (double*)
            fftw_malloc(sizeof(double) * smaller_count * line_size);
        double *bigger = (double*)
            fftw_malloc(sizeof(double) * bigger_count * line_size);
        double *anchor_line = (double*)
            fftw_malloc(sizeof(double) * line_size);

        // Backup the anchor line
        for (int i = 0; i < line_size; ++i)
            copy_line(anchor, where, anchor_line, 0);

        // Distribute among arrays
        int smaller_i=0, bigger_i=0;
        for (int i = 0; i < size; ++i) {
            if (get_prefix(i, where) > anchor_prefix)
                copy_line(i, where, bigger, bigger_i++);
            else if (i != anchor)
                copy_line(i, where, smaller, smaller_i++);
        }

        // Launch the sorts
        sort_jacobian(bigger, bigger_count);
        sort_jacobian(smaller, smaller_count);

        // Merge into the starting array
        for (int i = 0; i < size; ++i) {
            if (i < smaller_count)
                copy_line(i, smaller, where, i);
            else if (i == smaller_count) 
                copy_line(0, anchor_line, where, i);
            else 
                copy_line(i - smaller_count - 1, bigger, where, i);
        }

        // Free the memory
        fftw_free(smaller);
        fftw_free(bigger);
        fftw_free(anchor_line);
    }

    void Searcher::copy_line(size_t from_i, double *from, double *to, size_t to_i)
    {
        size_t line_size = 2 * integrator->size_real + 1;
        for (size_t i = 0; i < line_size; ++i)
            to[line_size * to_i + i] = from[line_size * from_i + i];
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
                jacobian[(size + 1) * j + i] =
                    (f_val1[j] - f_val2[j]) * inv_h;
        }
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
