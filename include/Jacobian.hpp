#ifndef turb_Jacobian_h
#define turb_Jacobian_h

namespace turb {
    class JacobianElement {
    private:
        size_t line_size;
        double *line;
        int get_prefix();
        bool free_at_destruction = false;

    public:
        size_t prefix();
        size_t size() { return line_size; }
        double* line() { return line;}
        void swap(JacobianElement *other);

        double operator[](int index);
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
        void swap(int i, int j);

        Jacobian(int m, int n);
        ~Jacobian();
    };
}


#endif
