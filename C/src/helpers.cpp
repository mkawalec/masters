#include <fftw3.h>
#include <string>

#include "helpers.hpp"
#include "exceptions.hpp"

namespace turb {
    const double e = -0.1;
    const double a = 0.125;
    const double b = -0.004;
    const double D = 40;

    const double R = 1.04;

    double l2_norm(double *array, size_t size)
    {
        double norm = 0.0;
        for (size_t i = 0; i < size; ++i) 
            norm += pow(array[i], 2);

        return sqrt(norm);
    }

    double l2_norm_cpx(fftw_complex *array, size_t size, size_t start_i)
    {
        double norm = 0.0;

        for (size_t i = start_i; i < size; ++i) {
            double *value = array[i];
            norm += pow(*value, 2) + pow(*(value + 1), 2);
        }
        return sqrt(norm);
    }
}
