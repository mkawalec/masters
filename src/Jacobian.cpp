#ifndef turb_Jacobian_cpp
#define turb_Jacobian_cpp

#include "exceptions.hpp"
#include "helpers.hpp"
#include "JacobianElement.hpp"

#include <fftw3.h>
#include <iostream>
#include <cmath>

namespace turb {

    template <typename T>
    JacobianElement<T>& Jacobian<T>::operator[] (int index)
    {
        if (index >= y) throw OutOfBounds();
        return elements[index];
    }

    template <typename T>
    Jacobian<T> Jacobian<T>::operator* (Jacobian<T> *second)
    {
        if (x != second->x || y != second->y)
            throw OutOfBounds("The dimensions don't match");

        for (int i = 0; i < y; ++i) {
            for (int j = 0; j < x; ++j) {
                tmp_line[j] = 0;

                for (int k = 0; k < x; k++)
                    tmp_line[j] += (*this)[i][k] * (*second)[k][i];
            }

            // Copy over the tmp_line
            for (int j = 0; j < x; ++j)
                (*this)[i][j] = tmp_line[j];
        }

        return *this;
    }

    template <typename T>
    void Jacobian<T>::swap_lines(int i, int j)
    {
        if (i >= y || j >= y) throw OutOfBounds();
        if (i == j) return;

        elements[i].swap(&elements[j]);
    }

    template <typename T>
    int Jacobian<T>::max_arg(int k)
    {
        double value = -1;
        int max_index = 0;

        for (int i = k; i < y; ++i) {
            if (fabs(elements[i][k]) > value) {
                value = fabs(elements[i][k]);
                max_index = i;
            }
        }

        return max_index;
    }


    template <typename T>
    Jacobian<T>::Jacobian(int m, int n)
    {
        y = m;
        x = n;

        jacobian = (T*) fftw_malloc(sizeof(T) * y * x);
        tmp_line = (T*) fftw_malloc(sizeof(T) * x);
        elements = (JacobianElement<T>*)
            fftw_malloc(sizeof(JacobianElement<T>) * y);

        for (int i = 0; i < y; ++i)
            new (&elements[i]) JacobianElement<T>(&jacobian[i*x], x);
    }

    template <typename T>
    Jacobian<T>::~Jacobian()
    {
        fftw_free(elements);
        fftw_free(jacobian);
    }
}

#endif
