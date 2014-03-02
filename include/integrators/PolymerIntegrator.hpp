#ifndef turb_PolymerIntegrator_h
#define turb_PolymerIntegrator_h

#include "helpers.hpp"
#include "Integrator.hpp"
#include "Jacobian.hpp"

#include <fftw3.h>
#include <fstream>
#include <complex>


namespace turb {

    /*! \brief The class doing actual computation through
     *      integrating a non-linear PDE
     *
     *  It is NOT thread-safe, don't try to use the same
     *  Integrator instance in multiple threads.
     */
    class PolymerIntegrator : Integrator {

    private:
        void allocate(size_t dim_power, double timestep, double domain);
        double lambda;

        matrix<std::complex<double> > *L_u_tau = NULL;
        fftw_complex *c_tau, *c_dtau, *tmp_u, *tmp_tau;

        fftw_plan e_tau, i_tau, f_tau, f_dtau, b_tau;

    protected:
        void initialize_operators();

        /// Computes one timestep of the nonlinear transform
        void compute_nonlinear();
        void compute_linear();

        virtual void override_initialize();

        virtual void initialize_function(double x, double *results);
        /*! \brief The overloadable nonlinear transform applied
         *      in the real (non-Fourier) space
         */
        virtual void nonlinear_transform(size_t i, double *results);

        // The much needed temporary array of two doubles
        double *temp_array = NULL;

        void set_options();
        void initialize(size_t dim_power, double timestep, double domain=2*M_PI);

    public:
        PolymerIntegrator();
        ~PolymerIntegrator();
        void clear(size_t dim_power, double dt, double domain_size=2*M_PI);

        Integrator* clone() const { return new PolymerIntegrator(*this); }

        void forward_transform();
    };
}

#endif
