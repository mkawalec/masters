#ifndef turb_helpers_h
#define turb_helpers_h

#include <fftw3.h>
#include <cmath>
#include <list>
#include <string>

namespace turb {
    extern const double e;
    extern const double a;
    extern const double b;
    extern const double R;
    extern const double D;
    extern double domain_size;

    double l2_norm(double *array, size_t size);
    double l2_norm_cpx(fftw_complex *array, size_t size, size_t start_i=0);

    inline void normalize(double *array, size_t size)
    {
        double norm_factor = 1 / sqrt(size);
        for (size_t i = 0; i < size; ++i)
            *(array + i) *= norm_factor;
    }

    template <typename T>
    void unused(T &&) { }

    template <typename T>
    class Base {
    protected:
        void unregister(T *instance);

    public:
        static std::list<T*> available;
        std::string name;
        std::string description;

        static T* choose(std::string name);
    };
}

#endif
