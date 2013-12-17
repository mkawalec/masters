#include "Jacobian.hpp"
#include "exceptions.hpp"
#include "Searcher.hpp"
#include "helpers.hpp"

#include <fftw3.h>

namespace turb {

    JacobianElement::prefix()
    {
        return get_prefix();
    }

    JacobianElement::swap(JacobianElement *other) 
    {
        double *tmp = line;
        line = other->line;
        other->line = tmp;
    }

    int JacobianElement::get_prefix()
    {
        int zero_count = 0;
        for (int i = 0; i < line_size; ++i) {
            if (!fuzzy_eql(where[line_size * line_index + i], 0)) break;
            ++zero_count;
        }
        return zero_count;
    }

    double JacobianElement::operator[](int index)
    {
        if (index >= line_size) throw OutOfBounds();
        return line[index];
    }

    double* JacobianElement::operator=(double *ptr)
    {
        line = ptr;
    }

    JacobianElement::JacobianElement(int line_size)
    {
        line = (double *) fftw_malloc(sizeof(double) * line_size);
        line_size = line_size;
        free_at_destruction = true;
    }

    JacobianElement::JacobianElement(double *start, size_t size)
    {
        line = start;
        line_size = size;
    }

    JacobianElement::~JacobianElement()
    {
        if (free_at_destruction) fftw_free(line);
    }

    JacobianElement Jacobian::operator[](int index)
    {
        if (index >= y) throw OutOfBounds();
        return elements[index];
    }

    void Jacobian::swap(int i, int j)
    {
        if (i >= y || j >= y) throw OutOfBounds();
        elements[i].swap(&elements[j]);
    }

    Jacobian::Jacobian(int m, int n)
    {
        y = m;
        x = n;

        jacobian = (double*) fftw_malloc(sizeof(double) * y * x);
        elements = (JacobianElement*) 
            fftw_malloc(sizeof(JacobianElement) * y);

        for (int i = 0; i < y; ++i)
            new (&elements[i]) JacobianElement(&jacobian[i*x], x);
    }

    Jacobian::~Jacobian()
    {
        for (int i = 0; i < y; ++i)
            delete &elements[i];
        fftw_free(jacobian);
    }
}


