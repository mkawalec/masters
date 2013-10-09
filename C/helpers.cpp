#include <fftw3.h>
#include <complex>
#include <string.h>

#include "helpers.h"


const double e = -0.1;
const double a = 0.125;
const double b = -0.004;
const double D = 40;

const double R = 0.95;
double domain_size = 24 * M_PI;

double l2_norm(double *array, size_t size)
{
        double norm = 0.0;
        for (size_t i = 0; i < size; ++i) 
                norm += pow(array[i], 2);

        return sqrt(norm);
}

double l2_norm_cpx(fftw_complex *array, size_t size)
{
        double norm = 0.0;

        for (size_t i = 0; i < size; ++i) {
                double *value = array[i];
                norm += pow(*value, 2) + pow(*(value + 1), 2);
        }
        return sqrt(norm);
}

Data_pointers allocate_precompute(unsigned long int dim_power, double dt)
{
        /* The 2^(size_t domain_size) is the size of the output array 
         * INCLUDING the padding, since FFTW is reportedly faster 
         * for arrays with dimensions being a power of 2.
         */
        Data_pointers program_data;
        program_data.size_real = pow(2, dim_power);
        program_data.size_complex = pow(2, dim_power) / 2 + 1;

        program_data.u = (double*) fftw_malloc(program_data.size_real * sizeof(double));
        program_data.v = (double*) fftw_malloc(program_data.size_real * sizeof(double));
        program_data.c_u = fftw_alloc_complex(program_data.size_complex);
        program_data.c_v = fftw_alloc_complex(program_data.size_complex);
        // Derivative in Fourier space
        program_data.dc_u = fftw_alloc_complex(program_data.size_complex);
        // Real (non-fourier space) derivative
        program_data.du = (double*) fftw_malloc(program_data.size_real * sizeof(double));

        // The operators acting on u and v
        program_data.Lu = (double*) fftw_malloc(program_data.size_complex * 
                                        sizeof(double) * 2.0 / 3);
        program_data.Lv = (double*) fftw_malloc(program_data.size_complex * 
                                        sizeof(double) * 2.0 / 3);


        // Initial plans
        program_data.i_u = fftw_plan_dft_r2c_1d(program_data.size_real,
                        program_data.u, program_data.c_u,
                        FFTW_MEASURE);
        program_data.i_v = fftw_plan_dft_r2c_1d(program_data.size_real,
                        program_data.v, program_data.c_v,
                        FFTW_MEASURE);
        // End (final) plans
        program_data.e_u = fftw_plan_dft_c2r_1d(program_data.size_real,
                        program_data.c_u, program_data.u,
                        FFTW_MEASURE | FFTW_PRESERVE_INPUT);
        program_data.e_v = fftw_plan_dft_c2r_1d(program_data.size_real,
                        program_data.c_v, program_data.v,
                        FFTW_MEASURE | FFTW_PRESERVE_INPUT);

        // Forwad plans
        program_data.f_u = fftw_plan_dft_c2r_1d(program_data.size_real,
                        program_data.c_u, program_data.u, FFTW_MEASURE);
        program_data.f_v = fftw_plan_dft_c2r_1d(program_data.size_real,
                        program_data.c_v, program_data.v, FFTW_MEASURE);
        program_data.f_du = fftw_plan_dft_c2r_1d(program_data.size_real,
                        program_data.dc_u, program_data.du, FFTW_MEASURE);
        // Backward plans
        program_data.b_u = fftw_plan_dft_r2c_1d(program_data.size_real,
                        program_data.u, program_data.c_u, FFTW_MEASURE);
        program_data.b_v = fftw_plan_dft_r2c_1d(program_data.size_real,
                        program_data.v, program_data.c_v, FFTW_MEASURE);

        compute_linear_operators(&program_data, dt);
        return program_data;
}

void compute_linear_operators(Data_pointers *program_data, double dt)
{
        /* Computes the linear operators acting on
         * u and v. TODO: Change to the C implementation
         * of the complex numbers and check how the performance
         * changes.
         */
        for (size_t i = 0; i < program_data->size_complex * 2.0 / 3; ++i) {
                double k = (double) i * 2 * M_PI / domain_size;
                double Lx, Ly;

                Lx = - pow(k, 4) + pow(k, 2) - (1 - e);
                Ly = - D * pow(k, 2) - 1;

                *(program_data->Lu + i) = (1 + 0.5 * dt * Lx) / (1 - 0.5 * dt * Lx);
                *(program_data->Lv + i) = (1 + 0.5 * dt * Ly) / (1 - 0.5 * dt * Ly);
        }
}

void initialize_modes_outputs(FILE ***outputs, Data_pointers *program_data)
{
        *outputs = (FILE**) malloc(sizeof(FILE*) * program_data->size_complex * 2/3);
        for (size_t j = 0; j < program_data->size_complex * 2.0 / 3; ++j) {
                /* Create the filenames output_n where n is
                 * the number of the wave mode height in that 
                 * file
                 */
                char num[10];
                char filename[20] = "output_";
                sprintf(num, "%ld", j);
                strcat(filename, num);

                *(*outputs + j) = (FILE*) malloc(sizeof(FILE*));
                *(*outputs + j) = fopen(filename, "w");
        }
}

void print_modes(FILE ***outputs, Data_pointers *program_data, double timestamp)
{
        for (size_t j = 0; j < program_data->size_complex * 2.0 / 3; ++j)
                fprintf(*(*outputs + j), "%f %f\n", timestamp, 
                                *(double*) program_data->c_u[j]);
}


