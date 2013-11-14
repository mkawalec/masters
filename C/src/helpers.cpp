#include "helpers.hpp"
#include "exceptions.hpp"

#include <fftw3.h>
#include <string>
#include <ap.h>


namespace turb {

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

    void e_x(const alglib::real_1d_array &c, const alglib::real_1d_array &x, 
            double &func, void *ptr) 
    {
        // this callback calculates f(c,x)=exp(-c0*sqr(x0))
        // where x is a position on X-axis and c is adjustable parameter
        func = exp(-c[0] * x[0]);
        unused(ptr);
    }

}
