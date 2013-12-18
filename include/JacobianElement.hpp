#ifndef turb_JacobianElement_h
#define turb_JacobianElement_h

#include "WriteCheck.hpp"

#include <utility>
#include <cstdlib>

namespace turb {
    template<typename T>
    class JacobianElement {
    private:
        size_t line_size;
        T *line;
        int get_prefix();
        bool free_at_destruction = false;
        bool dirty = false;
        int prefix_value = -1;

    public:
        int prefix();
        size_t size() { return line_size; }
        T* get_line() { return line;}
        T* swap(JacobianElement *other);

        T const& operator[](int index) const;
        WriteCheck<JacobianElement, T> operator[](int index);
        void update_state(int index);

        T* operator=(T *ptr);

        JacobianElement<T>(int line_size);
        JacobianElement<T>(T *start, size_t size);
        ~JacobianElement<T>();
    };
}

#include "JacobianElement.cpp"

#endif
