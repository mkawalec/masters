#include "searchers/SimpleSearcher.hpp"
#include "Integrator.hpp"
#include "exceptions.hpp"
#include "helpers.hpp"

#include <fstream>
#include <cstdlib>

namespace turb {

    void SimpleSearcher::allocate(Integrator *integrator)
    {
        this->integrator = integrator;

        srand(1234);

        f = (double*)
            fftw_malloc(sizeof(double) * 2 * integrator->size_real);
        du = (double*)
            fftw_malloc(sizeof(fftw_complex) * integrator->size_real);
        d_cu = (fftw_complex*)
            fftw_malloc(sizeof(fftw_complex) * integrator->size_complex);

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

        dx = (double*)
            fftw_malloc(sizeof(double) * 2 * integrator->size_real);

        fftw_import_wisdom_from_filename(get_wisdom_filename().c_str());
        du_c = fftw_plan_dft_r2c_1d(integrator->size_real, f, d_cu,
                FFTW_MEASURE);
        du_r = fftw_plan_dft_c2r_1d(integrator->size_real, d_cu, du,
                FFTW_MEASURE);

        d2v_c = fftw_plan_dft_r2c_1d(integrator->size_real,
                f + integrator->size_real, d2_cv, FFTW_MEASURE);
        d2v_r = fftw_plan_dft_c2r_1d(integrator->size_real, d2_cv, d2_v,
                FFTW_MEASURE);

        d2u_r = fftw_plan_dft_c2r_1d(integrator->size_real, d2_cu, d2_u,
                FFTW_MEASURE);

        d4u_r = fftw_plan_dft_c2r_1d(integrator->size_real, d4_cu, d4_u,
                FFTW_MEASURE);

        fftw_export_wisdom_to_filename(get_wisdom_filename().c_str());
    }

    void SimpleSearcher::init()
    {
        for (size_t i = 0; i < integrator->size_real; ++i) {
            f[i] = integrator->u[i];
            f[i + integrator->size_real] = integrator->v[i];
        }
    }

    void SimpleSearcher::F(double *__restrict__ input,
                           double *__restrict__ result)
    {
        double inv_domain_size = 1 / integrator->domain_size;
        double inv_size = 1 / (double) integrator->size_real;
        double tmp[2];

        // Execute the plan transforming u/v to complex form
        fftw_execute_dft_r2c(du_c, input, d_cu);
        fftw_execute_dft_r2c(d2v_c, input + integrator->size_real, d2_cv);
        for (size_t i = 0; i < integrator->size_complex; ++i) {
            double *d2u = d2_cu[i],
                   *d4u = d4_cu[i],
                   *du = d_cu[i];

            *d2u = *du;
            *(d2u + 1) = *(du + 1);
            *d4u = *du;
            *(d4u + 1) = *(du + 1);
            std::cout << *d2u << std::endl;
        }


        // Computing the derivative
        for (size_t i = 0; i < integrator->size_complex; ++i) {
            double k = (double) i * 2 * M_PI * inv_domain_size;
            double *du = d_cu[i];
            double *d2v = d2_cv[i];
            double *d2u = d2_cu[i];
            double *d4u = d4_cu[i];

            if (i < integrator->size_complex * 2.0 / 3) {
                tmp[0] = -k * *(du + 1);
                tmp[1] = k * *du;

                *du = tmp[0];
            std::cout << *du << std::endl;
                *(du + 1) = tmp[1];
                *d2v *= -pow(k, 2);
                *(d2v + 1) *= -pow(k, 2);
                *d2u *= -pow(k, 2);
                *(d2u + 1) *= -pow(k, 2);
                *d4u *= pow(k, 4);
                *(d4u + 1) *= pow(k, 4);
            } else {
                *du = 0;
                *(du + 1) = 0;
                *d2v = 0;
                *(d2v + 1) = 0;
                *d2u = 0;
                *(d2u + 1) = 0;
                *d4u = 0;
                *(d4u + 1) = 0;
            }
        }

        // Transform the derivative back into the real form
        fftw_execute(du_r);
        fftw_execute(d2v_r);
        fftw_execute(d2u_r);
        fftw_execute(d4u_r);

        for (int i = 0; (unsigned)i < integrator->size_real; ++i) {
            d2_v[i] *= inv_size;
            d2_u[i] *= inv_size;
            d4_u[i] *= inv_size;
            du[i] *= inv_size;
        }

        compute_F(input, result);
    }

    void SimpleSearcher::compute_F(double *__restrict__ input,
                                   double *__restrict__ result)
    {

        int size = integrator->size_real;
        for (int i = 0; i < size; ++i) {
            double u = input[i];
            double v = input[i + size];

            result[i] = -d4_u[i] - 2 * d2_u[i] +
                (-1 + integrator->e + (1 - integrator->e) *
                (integrator->a * v + integrator->b * pow(v, 2)) +
                du[i]) * u;

            result[i + size] = -v + integrator->D * d2_v[i] + integrator->R * pow(u, 2);
        }
    }

    SimpleSearcher::~SimpleSearcher()
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
        fftw_free(dx);
    }

    void SimpleSearcher::check_verify()
    {
        if (check_filename == "") return;

        std::ifstream input_file(check_filename);
        allocate(integrator);
        if (not input_file.is_open()) {
            throw ProgramDeathRequest("The file " +
                    check_filename + " is not readable!");
        }

        std::string tmp_line;
        for (int i = 0; (unsigned)i < 2 * integrator->size_real; ++i) {
            std::getline(input_file, tmp_line);
            const char *tmp_str = tmp_line.c_str();
            char *first_end;
            strtod(tmp_str, &first_end);
            f[i] = strtod(first_end, NULL);

            std::cout << i << " " << f[i] << std::endl;
        }
        std::cout << std::endl;

        F(f, f_val1);

        double norm = l2_norm(f_val1, integrator->size_real);
        std::cout << norm << std::endl;
        if (norm < 0.03) {
            throw ProgramDeathRequest("The file " +
                    check_filename + " is correct");
        }
        else {
            throw ProgramDeathRequest("The norm for file " +
                    check_filename + " is " + std::to_string(norm) + " and is too big!");
        }
    }

}


