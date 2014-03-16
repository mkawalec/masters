#ifndef turb_Integrator_h
#define turb_Integrator_h

#include "Base.hpp"
#include "Searcher.hpp"

#include <memory>
#include <string>
#include <boost/program_options/options_description.hpp>
namespace po = boost::program_options;

namespace turb {
    class Searcher;

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
        virtual void initialize(size_t dim_power,
                double timestep, double domain=2*M_PI) = 0;


    public:
        Integrator() { Integrator::available.push_back(this);}
        virtual ~Integrator() { };
        virtual void clear(size_t dim_power, double dt, double domain_size=2*M_PI) = 0;

        virtual void apply_step();

        virtual void forward_transform() = 0;
        size_t size_real, size_complex;

        double dt, domain_size, dim_power,
               e, a, b, D, R, start_mult;
        double *u = NULL, *v = NULL, *du = NULL, *Lu = NULL, *Lv = NULL,
               *tau = NULL, *dtau = NULL;

        fftw_plan e_u, e_v, i_u, i_v,
                  f_u, f_v, f_du, b_u, b_v;
        fftw_complex *c_u = NULL, *c_v = NULL, *dc_u = NULL;

        Searcher *searcher = NULL;
        void setup_searcher();
        void setup_searcher(int argc, char *argv[]);

        std::string selected_searcher;

        virtual std::vector<double> get_norms(std::vector<double> coords);

        void parse_params(int argc, char *argv[]);
        std::string additional_info();

        virtual Integrator* clone() const = 0;

        // If set to true, static points will be searched for
        bool search;
    };
}

#endif
