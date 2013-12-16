#include "helpers.hpp"
#include "exceptions.hpp"

#include <fftw3.h>
#include <string>
#include <ap.h>
#include <vector>


namespace turb {

    double l2_norm(double *array, size_t size, size_t start_i)
    {
        double norm = 0.0;
        for (size_t i = start_i; i < size; ++i) 
            norm += pow(array[i], 2);

        return sqrt(norm);
    }

    double l2_norm(fftw_complex *array, size_t size, size_t start_i)
    {
        double norm = 0.0;

        for (size_t i = start_i; i < size; ++i) {
            double *value = array[i];
            norm += pow(*value, 2) + pow(*(value + 1), 2);
        }
        return sqrt(norm);
    }

    double l2_norm(std::vector<double>::iterator start,
                   std::vector<double>::iterator end) 
    {
        double norm = 0.0;
        for (; start != end; ++start)
            norm += pow(*start, 2);
        return sqrt(norm);
    }


    void e_x(const alglib::real_1d_array &c, const alglib::real_1d_array &x, 
            double &func, void *ptr) 
    {
        func = exp(-c[0] * x[0]);
        unused(ptr);
    }

}
