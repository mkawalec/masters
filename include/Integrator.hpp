#ifndef turb_Integrator_h
#define turb_Integrator_h

#include "Base.hpp"

#include <memory>
#include <boost/program_options/options_description.hpp>
namespace po = boost::program_options;

namespace turb {

    class Integrator : public Base<Integrator> {
    protected:
        void initialize_operators() = 0;

        /// Computes one timestep of the nonlinear transform
        virtual void compute_nonlinear();
        virtual void compute_linear();

        virtual void override_initialize();

        virtual void initialize_function(double x, double *results) = 0;
        /*! \brief The overloadable nonlinear transform applied
         *      in the real (non-Fourier) space
         */
        virtual void nonlinear_transform(size_t i, double *results) = 0;

        // The much needed temprorary array of two doubles
        double *temp_array;

        virtual void set_options() { };
        std::shared_ptr<po::options_description> options = NULL;

    public:
        virtual ~Integrator();
        virtual void initialize() = 0;

        virtual void apply_step();
        virtual void forward_transform() = 0;
        size_t size_real, size_complex;

        double dt, domain_size;

        fftw_plan e_u, e_v, i_u, i_v,
                  f_u, f_v, f_du, b_u, b_v;

        virtual void parse_params(int argc, const char *argv[]);
        virtual Integrator* clone() const = 0;
    };
}
