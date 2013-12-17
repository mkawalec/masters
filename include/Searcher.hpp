#ifndef turb_Searcher_h
#define turb_Searcher_h

#include "Integrator.hpp"
#include <vector>

namespace turb {

    class Searcher {
    private:
        Integrator *integrator;
        double *f, *du, *jacobian, 
               *f_val1, *f_val2, *dx;
        fftw_complex *d_cu;

        size_t iterations = 50;
        double threshold = 0.0001;
        double overflow = 1e20;
        double h = 0.0001;

        fftw_plan du_c, du_r;

        void get_jacobian();
        void F(double *input, double *result);

        // Gaussian elimination, must operate on 
        // non-singular matrices, will break badly if
        // a singular matrix is provided
        void gauss(double *f, double *result);

        // An in-place (with additional storage needed) 
        // quicksort based on zeros count. If a more memory-efficient
        // implementation turns out to be needed, it will be added
        void sort_jacobian(double *where, int size);
        int get_prefix(size_t line_index, double *where);
        void copy_line(size_t from_i, double *from, double *to, size_t to_i);

    public:
        Searcher(Integrator *integrator);
        ~Searcher();

        std::vector<double> run();
    };

    class JacobianElement {
    private:
        size_t prefix, line_size;
        double *line;

    public:
        size_t prefix();
        size_t size();
        double* line() { return line;}
        void swap(Element *other);

        Element(int line_size);
        double operator[](int index);
        double* operator=(double *ptr);
    };

    class Jacobian {
    private:
        double *jacobian;

    public:
        Jacobian(int m, int n);
        JacobianElement operator[](int index);
    };
}

#endif
