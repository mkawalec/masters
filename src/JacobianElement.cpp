#ifndef turb_JacobianElement_cpp
#define turb_JacobianElement_cpp

#include "exceptions.hpp"
#include "helpers.hpp"
#include "WriteCheck.hpp"

namespace turb {
    template <typename T>
    int JacobianElement<T>::prefix()
    {
        if (prefix_value == -1 || dirty) {
            prefix_value = get_prefix();
            dirty = false;
        }
        return prefix_value;
    }

    template <typename T>
    int JacobianElement<T>::get_prefix()
    {
        int zero_count = 0;
        for (int i = 0; (unsigned)i < line_size; ++i) {
            if (!fuzzy_eql(line[i], 0)) break;
            ++zero_count;
        }
        return zero_count;
    }

    template <typename T>
    T* JacobianElement<T>::swap(JacobianElement *other)
    {
        T *tmp = line;
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

    template <typename T>
    T const& JacobianElement<T>::operator[](int index) const
    {
        if ((unsigned)index >= line_size) throw OutOfBounds();
        return line[index];
    }

    template <typename T>
    WriteCheck<JacobianElement<T>, T> JacobianElement<T>::operator[](int index)
    {
        return WriteCheck<JacobianElement, T>(this, &line[index], index);
    }

    template <typename T>
    void JacobianElement<T>::update_state(int index)
    {
        if (index > prefix_value) return;
        /*if (index == prefix_value && fuzzy_eql(line[index], 0)) {
            ++prefix_value;
            return;
        }*/

        dirty = true;
    }

    template <typename T>
    T* JacobianElement<T>::operator=(T *ptr)
    {
        return line = ptr;
    }

    template <typename T>
    JacobianElement<T>::JacobianElement(int line_size)
    {
        line = (T*) fftw_malloc(sizeof(T) * line_size);
        this->line_size = line_size;
        free_at_destruction = true;
    }

    template <typename T>
    JacobianElement<T>::JacobianElement(T *start, size_t size)
    {
        line = start;
        line_size = size;
    }

    template <typename T>
    JacobianElement<T>::~JacobianElement()
    {
        if (free_at_destruction) fftw_free(line);
    }
}

#endif
