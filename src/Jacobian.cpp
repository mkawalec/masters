#include "Jacobian.hpp"
#include "exceptions.hpp"
#include "Searcher.hpp"
#include "helpers.hpp"

#include <fftw3.h>
#include <iostream>

namespace turb {

    int JacobianElement::prefix()
    {
        if (prefix_value == -1 || dirty) {
            prefix_value = get_prefix();
            dirty = false;
        }
        return prefix_value;
    }

    double* JacobianElement::swap(JacobianElement *other) 
    {
        double *tmp = line;
        line = other->line;
        other->line = tmp;

        bool b_tmp = dirty;
        dirty = other->dirty;
        other->dirty = b_tmp;

        b_tmp = free_at_destruction;
        free_at_destruction = other->free_at_destruction;
        other->free_at_destruction = b_tmp;

        size_t s_tmp = line_size;
        line_size = other->line_size;
        other->line_size = s_tmp;

        int i_tmp = prefix_value;
        prefix_value = other->prefix_value;
        other->prefix_value = i_tmp;

        return line;
    }

    int JacobianElement::get_prefix()
    {
        int zero_count = 0;
        for (int i = 0; (unsigned)i < line_size; ++i) {
            if (!fuzzy_eql(line[i], 0)) break;
            ++zero_count;
        }
        return zero_count;
    }

    double const& JacobianElement::operator[](int index) const
    {
        if ((unsigned)index >= line_size) throw OutOfBounds();
        return line[index];
    }

    WriteCheck<JacobianElement, double> JacobianElement::operator[](int index)
    {
        return WriteCheck<JacobianElement, double>(this, &line[index], index);
    }

    void JacobianElement::update_state(int index)
    {
        if (prefix_value != -1 && index >= prefix_value) return;
        dirty = true;
    }

    double* JacobianElement::operator=(double *ptr)
    {
        return line = ptr;
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

    void Jacobian::swap_lines(int i, int j)
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
        fftw_free(elements);
        fftw_free(jacobian);
    }
}


