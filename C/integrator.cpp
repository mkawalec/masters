#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include <complex>
#include <math.h>
#include <string.h>

#include "helpers.h"


void compute_nonlinear(Data_pointers *prog_data, double dt) 
{
        size_t i = 0;
        // Compute the derivative
        for (i = 0; i < prog_data->size_complex; ++i) {
                double *du = prog_data->dc_u[i];
                double *c_u = prog_data->c_u[i];

                *du = - (double) i * 2 * M_PI / domain_size * *(c_u + 1);
                *(du + 1) = (double) i * 2 * M_PI / domain_size * *c_u;
        }


        fftw_execute(prog_data->f_u);
        fftw_execute(prog_data->f_v);
        fftw_execute(prog_data->f_du);

        double *u = prog_data->u;
        double *v = prog_data->v;
        double *du = prog_data->du;
        // Multiply the parts according to the nonlinear equation
        for (i = 0; i < prog_data->size_real; ++i) {
                double temp_u;

                /* Calculating the u term.
                 * The u value is saved as a temprorary value
                 * as we want to avoid overwriting u with the
                 * new value before v is computed
                 */
                temp_u = (1 - e) *
                        (a * *(v + i) + b * pow(*(v + i), 2)) * *(u + i) +
                        *(u + i) * *(du + i);

                // Calculating the v term
                *(v + i) += dt * R * pow(*(u + i), 2);

                // Saving the results
                *(u + i) += dt * temp_u;
        }

        fftw_execute(prog_data->b_u);
        fftw_execute(prog_data->b_v);

        // Normalize and nullify the padding
        for (i = 0; i < prog_data->size_complex; ++i) {
                double *u = prog_data->c_u[i];
                double *v = prog_data->c_v[i];

                if (i > 2.0/3 * prog_data->size_complex) {
                        *u = 0;
                        *(u + 1) = 0;
                        *v = 0;
                        *(v + 1) = 0;
                } else {
                        *u /= prog_data->size_real;
                        *(u + 1) /= prog_data->size_real;
                        *v /= prog_data->size_real;
                        *(v + 1) /= prog_data->size_real;
                }
        }
}


void compute_linear(Data_pointers *prog_data)
{
        for (size_t i = 0; i < prog_data->size_complex; ++i) {
                double *c_u = prog_data->c_u[i];
                double *c_v = prog_data->c_v[i];

                /* To compute u and v, multiply both parts of the number
                 * by the (real) value of the corresponding precomputed
                 * linear operator
                 */
                *(c_u) *= *(prog_data->Lu + i);
                *(c_u + 1) *= *(prog_data->Lu + i);
                *(c_v) *= *(prog_data->Lv + i);
                *(c_v + 1) *= *(prog_data->Lv + i);
        }
}

void initialize(Data_pointers *program_data)
{
        /* Initialize the data, according to the starting conditions
         * from the Physica D 238 paper
         */
        for (size_t i = 0; i < program_data->size_real; ++i) {
                double x = (double) i / program_data->size_real * domain_size;
                program_data->u[i] = 2.0 * cos(x) + 0.03 * cos(11.0/12.0 * x);
                program_data->v[i] = 1.0;
        }

        /* On every stage we are interested in working
         * on the Fourier space
         */
        fftw_execute(program_data->i_u);
        fftw_execute(program_data->i_v);
}


int main(int argc, char *argv[]) 
{
        unsigned long int dim_power, steps;
        double end_time, dt;
        FILE *output_u, *output_v;
        char **char_ptr = NULL;

        if (argc < 5) {
                printf("To call this program properly use %s "
                                "dim_power steps domain_size end_time\n", argv[0]);
                return -1;
        }

        // Output files
        output_u = fopen("output_u", "w");
        output_v = fopen("output_v", "w");

        dim_power = strtoul(argv[1], char_ptr, 10);
        steps = strtoul(argv[2], char_ptr, 10);
        domain_size = strtod(argv[3], char_ptr);
        end_time = strtod(argv[4], char_ptr);
        dt = (double) end_time/steps;

        printf("Callibrating FFTW...\n");
        Data_pointers program_data = allocate_precompute(dim_power, dt);
        printf("done\n");

        printf("Applying starting conditions...\n");
        initialize(&program_data);   
        printf("done\n");            

        // Initialize modes outputs and open the files
        FILE **outputs;
        outputs = (FILE**) malloc(sizeof(FILE*) * program_data.size_complex);
        for (size_t j = 0; j < program_data.size_complex; ++j) {
                char num[10];
                char filename[20] = "output_";
                sprintf(num, "%ld", j);
                strcat(filename, num);
                outputs[j] = (FILE*) malloc(sizeof(FILE*));
                outputs[j] = fopen(filename, "w");
        }
                                     
        for (size_t i = 0; i < steps; ++i) {
                compute_nonlinear(&program_data, dt);
                compute_linear(&program_data);
                double timestamp = (double) end_time/steps * i;

                // Print the modes
                for (size_t j = 0; j < program_data.size_complex; ++j){
                        fprintf(outputs[j], "%lf %lf\n", timestamp, 
                                        *(double*) program_data.c_u[j]);
                }


                // Print the results to the output
                if (i%1 == 0) {
                        // Transform to the real basis
                        fftw_execute(program_data.e_u);
                        fftw_execute(program_data.e_v);
                        normalize(program_data.u, program_data.size_real);
                        normalize(program_data.v, program_data.size_real);


                        fprintf(output_u, "%lf %lf\n", timestamp, 
                                        l2_norm(program_data.u, program_data.size_real));
                        fprintf(output_v, "%lf %lf\n", timestamp, 
                                        l2_norm(program_data.v, program_data.size_real));
                }
        }

        fclose(output_u);
        fclose(output_v);
        for (size_t j = 0; j < program_data.size_complex; ++j)
                fclose(outputs[j]);
}
