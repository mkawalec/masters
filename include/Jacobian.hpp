#ifndef turb_Jacobian_h
#define turb_Jacobian_h

#include "JacobianElement.hpp"

#include <utility>
#include <cstdlib>

namespace turb {

    template <typename T>
    class Jacobian {
    private:
        T *jacobian;
        JacobianElement<T> *elements;

        // A temporary line, used for multiplication
        T *tmp_line;
        int x, y;

    public:
        JacobianElement<T> operator[] (int index);
        Jacobian<T> operator* (Jacobian<T> *second);

        void swap_lines(int i, int j);
        int max_arg(int k);

        std::pair<int, int> dims() { return std::pair<int, int>(y, x);}

        Jacobian<T>(int m, int n);
        ~Jacobian<T>();
    };

    template <typename T>
    using matrix = Jacobian<T>;
}

#include "Jacobian.cpp"

#endif
