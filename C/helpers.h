#include <fftw3.h>

extern const double e;
extern const double a;
extern const double b;
extern const double R;
extern const double D;
extern double domain_size;

typedef struct Data_pointers {
        /* The complex arrays are the intermediate arrays
         * used for calculating convolution in the Fourier
         * transformed basis
         */
        double *u, *v, *du, *Lu, *Lv;
        fftw_complex *c_u, *c_v, *dc_u;

        /* Forward and backward FFTW plans for both arrays.
         * They are small, so I will just store them on the stack,
         * why not.
         */
        fftw_plan e_u, e_v, i_u, i_v, 
                  f_u, f_v, f_du, b_u, b_v;

        size_t size_real, size_complex;
} Data_pointers;

double l2_norm(double *array, size_t size);
void normalize(double *array, size_t size);
Data_pointers allocate_precompute(unsigned long int dim_power, double dt);
void compute_linear_operators(Data_pointers *program_data, double dt);
void initialize_modes_outputs(FILE ***outputs, Data_pointers *program_data);
void print_modes(FILE ***outputs, Data_pointers *program_data, double timestamp);
