#ifndef turb_helpers_h
#define turb_helpers_h

#include <fftw3.h>
#include <cmath>

extern const double e;
extern const double a;
extern const double b;
extern const double R;
extern const double D;
extern double domain_size;

double l2_norm(double *array, size_t size);
double l2_norm_cpx(fftw_complex *array, size_t size);

inline void normalize(double *array, size_t size)
{
        double norm_factor = 1 / sqrt(size);
        for (size_t i = 0; i < size; ++i)
                *(array + i) *= norm_factor;
}

#endif
