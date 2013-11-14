#ifndef turb_helpers_h
#define turb_helpers_h

#include <fftw3.h>
#include <cmath>
#include <ap.h>

namespace turb {
    double l2_norm(double *array, size_t size);
    double l2_norm_cpx(fftw_complex *array, size_t size, size_t start_i=0);

    inline void normalize(double *array, size_t size)
    {
        double norm_factor = 1 / sqrt(size);
        for (size_t i = 0; i < size; ++i)
            *(array + i) *= norm_factor;
    }

    struct history {
        double time, u, v;
    };

    template <typename T>
    void unused(T &&) { }

    void e_x(const alglib::real_1d_array &c, const alglib::real_1d_array &x, 
            double &func, void *ptr);

}

#endif
