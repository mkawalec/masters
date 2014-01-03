#ifndef turb_Integrator_h
#define turb_Integrator_h

#include "Base.hpp"

#include <memory>
#include <boost/program_options/options_description.hpp>
namespace po = boost::program_options;

namespace turb {

    class Integrator : public Base<Integrator> {
    protected:
        virtual void initialize_operators() = 0;

        /// Computes one timestep of the nonlinear transform
        virtual void compute_nonlinear() = 0;
        virtual void compute_linear() = 0;

        virtual void override_initialize() = 0;

        virtual void initialize_function(double x, double *results) = 0;
        /*! \brief The overloadable nonlinear transform applied
         *      in the real (non-Fourier) space
         */
        virtual void nonlinear_transform(size_t i, double *results) = 0;

        // The much needed temporary array of two doubles
        double *temp_array;

        virtual void set_options() { };
        std::shared_ptr<po::options_description> options = NULL;

    public:
        virtual ~Integrator() { };
        virtual void initialize(size_t dim_power,
                double timestep, double domain=2*M_PI) = 0;
        virtual void clear() = 0;

        virtual void apply_step();

        virtual void forward_transform() = 0;
        size_t size_real, size_complex;

        double dt, domain_size, dim_power,
               e=-0.1, a=0.125, b=-0.004, D=40, R=1.04;
        double *u = NULL, *v = NULL, *du = NULL, *Lu = NULL, *Lv = NULL;

        fftw_plan e_u, e_v, i_u, i_v,
                  f_u, f_v, f_du, b_u, b_v;
        fftw_complex *c_u = NULL, *c_v = NULL, *dc_u = NULL;

        void parse_params(int argc, const char *argv[]);
        virtual Integrator* clone() const = 0;
    };
}

#endif
