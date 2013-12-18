#ifndef turb_helpers_h
#define turb_helpers_h

#include <fftw3.h>
#include <cmath>
#include <ap.h>
#include <vector>

namespace turb {

    /*! \brief Returns L2 norm of a real array.
     *  \param array an array of real numbers
     *  \param size array size, can be smaller than
     *      the size of a complete array
     *  \param start_i start the norm computation
     *      at i-th element
     *
     *  \return L2 norm(array)
     */
    double l2_norm(double *array, size_t size, size_t start_i=0);
    double l2_norm(long double *array, size_t size, size_t start_i=0);
    double l2_norm(__float128 *array, size_t size, size_t start_i=0);

    /*! \brief Returns L2 norm of a complex array
     *  \param array a complex-numbered array
     *  \param size array length
     *  \param start_i start the norm computation at i-th
     *      element
     *
     *  \return L2 norm(array)
     */
    double l2_norm(fftw_complex *array, size_t size, size_t start_i=0);

    double l2_norm(std::vector<double>::iterator start,
                   std::vector<double>::iterator end);

    /*! \brief Normalizes a FFT-transformed array
     *  \param array array of numbers
     *  \param size size of the array
     */
    inline void normalize(double *array, size_t size)
    {
        double norm_factor = 1 / sqrt(size);
        for (size_t i = 0; i < size; ++i)
            *(array + i) *= norm_factor;
    }

    const double EQL_ACCURACY = 1e-14;

    inline bool fuzzy_eql(double first, double second)
    {
        if (fabs(first-second) > EQL_ACCURACY)
            return false;
        return true;
    }


    /*! \brief Holds information about values of
     *      u and v at a given time.
     */
    struct history {
        double time, u, v;
    };

    /*! \brief Used to explicitly mark parameters
     *      as unused and mute compiler warnings.
     */
    template <typename T>
    void unused(T &&) { }

    /*! \brief Returns values of e^-cx
     *  \param c array of parameters
     *  \param x array of values
     *  \param func function value
     *  \param ptr an unused pointer
     */
    void e_x(const alglib::real_1d_array &c, const alglib::real_1d_array &x, 
            double &func, void *ptr);

    double current_time();
}

#endif
