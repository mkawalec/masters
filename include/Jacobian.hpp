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
        int x, y;

    public:
        JacobianElement<T> operator[](int index);
        void swap_lines(int i, int j);

        std::pair<int, int> dims() { return std::pair<int, int>(y, x);}

        Jacobian<T>(int m, int n);
        ~Jacobian<T>();
    };
}

#include "Jacobian.cpp"

#endif
