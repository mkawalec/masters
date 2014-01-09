#include "helpers.hpp"
#include "exceptions.hpp"

#include <mpi.h>
#include <fftw3.h>
#include <string>
#include "ap.h"
#include <vector>
#include <sys/time.h>
#include <iostream>


namespace turb {

    double l2_norm(double *array, size_t size, size_t start_i)
    {
        double norm = 0.0;
        for (size_t i = start_i; i < size; ++i)
            norm += pow(array[i], 2);

        return sqrt(norm);
    }

    double l2_norm(long double *array, size_t size, size_t start_i)
    {
        long double norm = 0.0;
        for (size_t i = start_i; i < size; ++i)
            norm += pow(array[i], 2);

        return sqrt(norm);
    }

    double l2_norm(__float128 *array, size_t size, size_t start_i)
    {
        __float128 norm = 0.0;
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

    double current_time()
    {
        timeval time;
        gettimeofday(&time, NULL);
        return time.tv_sec + time.tv_usec * (double)1e-6;
    }

    bool contains(std::vector<std::vector<double> > *collection,
                  std::vector<double> *element)
    {
        for (int i = 0; (unsigned)i < collection->size(); ++i) {
            bool same = true;
            for (int j = 0; (unsigned)j < (*collection)[i].size(); ++j) {
                if (!fuzzy_eql((*collection)[i][j], (*element)[j], 1e-2)) {
                    same = false;
                    break;
                }
            }
            if (same) return true;
        }
        return false;
    }

    std::string get_wisdom_filename()
    {
        int my_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        return "." + std::to_string(my_rank) + "wisdom";
    }

}
