#ifndef turb_Jacobian_h
#define turb_Jacobian_h

#include "WriteCheck.hpp"

#include <utility>
#include <cstdlib>

namespace turb {
    class JacobianElement {
    private:
        size_t line_size;
        double *line;
        int get_prefix();
        bool free_at_destruction = false;
        bool dirty = false;
        int prefix_value = -1;

    public:
        int prefix();
        size_t size() { return line_size; }
        double* swap(JacobianElement *other);

        double const& operator[](int index) const;
        WriteCheck<JacobianElement, double> operator[](int index);
        void update_state(int index);

        double* operator=(double *ptr);

        JacobianElement(int line_size);
        JacobianElement(double *start, size_t size);
        ~JacobianElement();
    };

    class Jacobian {
    private:
        double *jacobian;
        JacobianElement *elements;
        int x, y;

    public:
        JacobianElement operator[](int index);
        void swap_lines(int i, int j);

        std::pair<int, int> dims() { return std::pair<int, int>(y, x);}

        Jacobian(int m, int n);
        ~Jacobian();
    };
}

#endif
