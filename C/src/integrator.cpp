#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include <complex>
#include <math.h>

#include "helpers.h"


void compute_nonlinear(Data_pointers *prog_data, double dt) 
{
        size_t i = 0;
        double inv_domain_size = 1 / domain_size;
        // Compute the derivative
        for (i = 0; i < prog_data->size_complex; ++i) {
                double *du = prog_data->dc_u[i];
                double *c_u = prog_data->c_u[i];

                *du = - (double) i * 2 * M_PI * inv_domain_size * *(c_u + 1);
                *(du + 1) = (double) i * 2 * M_PI * inv_domain_size * *c_u;
        }


        fftw_execute(prog_data->f_u);
        fftw_execute(prog_data->f_v);
        fftw_execute(prog_data->f_du);

        // Normalize the transfor
        double *u = prog_data->u;
        double *v = prog_data->v;
        double *du = prog_data->du;
        double scale_factor = 1/sqrt(prog_data->size_real);

        // Not using normalize for performance reasons
        for (i = 0; i < prog_data->size_real; ++i) {
                *(u + i) *= scale_factor;
                *(v + i) *= scale_factor;
                *(du + i) *= scale_factor;
        }

        // Multiply the parts according to the nonlinear equation
        for (i = 0; i < prog_data->size_complex; ++i) {
                double temp_u;

                /* Calculating the u term.
                 * The u value is saved as a temprorary value
                 * as we want to avoid overwriting u with the
                 * new value before v is computed
                 */
                temp_u = (1 - e) *
                        (a * v[i] + b * pow(v[i], 2)) * u[i] +
                        u[i] * du[i];

                // Calculating the v term
                v[i] += dt * R * pow(u[i], 2);

                // Saving the results
                u[i] += dt * temp_u;
        }

        fftw_execute(prog_data->b_u);
        fftw_execute(prog_data->b_v);

        // Normalize and nullify the padding
        for (i = 0; i < prog_data->size_complex; ++i) {
                double *u = prog_data->c_u[i];
                double *v = prog_data->c_v[i];

                if (i >= prog_data->size_complex * 2.0 / 3) {
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


void compute_linear(Data_pointers *prog_data)
{
        for (size_t i = 0; i < prog_data->size_complex * 2.0 / 3; ++i) {
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
                program_data->u[i] = 2.0 * cos(x) + 0.03 * cos(11 * x / 12.0);
                program_data->v[i] = 0.0;
        }
        printf("%ld\n", program_data->size_real);
        printf("Initial L2 is: %f\n", l2_norm(program_data->u, program_data->size_real));

        /* On every stage we are interested in working
         * on the Fourier space
         */
        fftw_execute(program_data->i_u);
        fftw_execute(program_data->i_v);

        double scale_factor = 1 / sqrt(program_data->size_real);
        // Make the padding be the padding;)
        for (size_t i = 0; i < program_data->size_complex; ++i) {
                double *c_u = program_data->c_u[i];
                double *c_v = program_data->c_v[i];

                if (i >= program_data->size_complex * 2.0 / 3) {
                        *c_u = 0.0;
                        *(c_u + 1) = 0.0;
                        *c_v = 0.0;
                        *(c_v + 1) = 0.0;
                } else {
                        *c_u *= scale_factor;
                        *(c_u + 1) *= scale_factor;
                        *c_v *= scale_factor;
                        *(c_v + 1) *= scale_factor;
                }
        }

        FILE *start = fopen("output_start", "w");
        for (size_t i = 0; i < program_data->size_complex; ++i) {
                double *c_u = program_data->c_u[i];
                fprintf(start, "%ld %f %f\n", i, *c_u, *(c_u + 1));
        }
        fclose(start);
}


int main(int argc, char *argv[]) 
{
        unsigned long int dim_power;
        double end_time, dt, current_time;
        FILE *output;
#ifdef DEBUG
        FILE **outputs = NULL;
#endif
        char **char_ptr = NULL;

        if (argc < 4) {
                printf("To call this program properly use %s "
                                "dim_power dt end_time\n", argv[0]);
                return -1;
        }

        // Parse the commandline parameters
        dim_power = strtoul(argv[1], char_ptr, 10);
        dt = strtod(argv[2], char_ptr) * 0.5;
        end_time = strtod(argv[3], char_ptr);

        // Bind the output files
        output = fopen("output", "w");

        // Initialize everything
        printf("Callibrating FFTW...\n");
        Data_pointers program_data = allocate_precompute(dim_power, dt);
        printf("done\n");

        printf("Applying starting conditions...\n");
        initialize(&program_data);   
        printf("done\n");            

#ifdef DEBUG
        // Initialize modes outputs and open the files
        initialize_modes_outputs(&outputs, &program_data);
#endif

        current_time = 0.0;
        size_t i = 0;
        while(current_time < end_time) {
                compute_nonlinear(&program_data, dt);
                compute_linear(&program_data);

#ifdef DEBUG
                // Print the modes
                print_modes(&outputs, &program_data, current_time);
#endif

                // Print the L2 norm in complex space 
                fprintf(output, "%f %f %f\n", current_time,
                                l2_norm_cpx(program_data.c_u, program_data.size_complex),
                                l2_norm_cpx(program_data.c_v, program_data.size_complex));

                // Print the results to the output
                if (i%10 == 0) {
                        // Transform to the real basis
                        fftw_execute(program_data.e_u);
                        fftw_execute(program_data.e_v);
                        normalize(program_data.u, program_data.size_real);
                        normalize(program_data.v, program_data.size_real);

                        fprintf(output, "%f %f %f\n", current_time, 
                                        l2_norm(program_data.u, program_data.size_real),
                                        l2_norm(program_data.v, program_data.size_real));
                }

                current_time += 2 * dt;
                ++i;
        }

        fclose(output);

#ifdef DEBUG
        for (size_t j = 0; j < program_data.size_complex * 2.0 / 3; ++j)
                fclose(outputs[j]);
#endif
}
