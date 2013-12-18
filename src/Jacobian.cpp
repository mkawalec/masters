#ifndef turb_Jacobian_cpp
#define turb_Jacobian_cpp

#include "exceptions.hpp"
#include "helpers.hpp"
#include "JacobianElement.hpp"

#include <fftw3.h>
#include <iostream>

namespace turb {

    template <typename T>
    JacobianElement<T> Jacobian<T>::operator[](int index)
    {
        if (index >= y) throw OutOfBounds();
        return elements[index];
    }

    template <typename T>
    void Jacobian<T>::swap_lines(int i, int j)
    {
        if (i >= y || j >= y) throw OutOfBounds();
        elements[i].swap(&elements[j]);
    }

    template <typename T>
    Jacobian<T>::Jacobian(int m, int n)
    {
        y = m;
        x = n;

        jacobian = (T*) fftw_malloc(sizeof(T) * y * x);
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
