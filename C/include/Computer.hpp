#ifndef turb_Computer_h
#define turb_Computer_h

#include "Serializer.hpp"
#include "Integrator.hpp"
#include "Base.hpp"

#include <thread>
#include <memory>
#include <boost/program_options/options_description.hpp>
namespace po = boost::program_options;

namespace turb {

    /*! \brief Base of Computer-type workers that implement
     *      computation logic
     *
     *  These workers manage binding the Integrators,
     *  setting up the simulation and then performing
     *  the computation. 
     *
     *  They also provide some additional bookkeeping
     *  functionality like parsing of additional command
     *  line parameters
     */
    class Computer : public Base<Computer> {
    protected:
        /// Currently used Integrator
        Integrator *integrator = NULL;

        /// \brief Function doing the actual computation.
        virtual void compute() = 0;

        /*! \brief Initializes integrator with constants with
         *      proper values.
         */ 
        void set_constants();

        /*! \brief Implements class-specific options setting logic.
         */
        virtual void set_options() { };

        /*! \brief If set points to a description of 
         *      the class-specific command line options
         */
        std::shared_ptr<po::options_description> options = NULL;

        /*! \brief If set and no serializer is set at 
         *      set-up stage this serializer will be used.
         *
         *  Note that if the specified serializer is 
         *  unavailable at runtime, InstanceNotFound will
         *  still be thrown.
         */
        std::string suggested_serializer;

        /*! \brief If called and the serializer is unset,
         *      set it to a sensible value.
         */
        void set_serializer();

    public:
        /*! \brief Execute to start the simulation.
         *  \return an std::thread object that runs the integration
         *
         *  Note that you will want to call join() on the returned
         *  thread to actually do any computation whatsoever
         */
        std::thread run() { return std::thread(&Computer::compute, this); }

        /// \brief Currently used Serializer
        Serializer *serializer = NULL;

        /// \brief Parameters needing to be set before run
        size_t print_every=100, samples=7, runs=2000;
        double end_time=2000, dt=0.0005, e, a, b, D, 
               R, domain_size=20*M_PI, threshold=300;
        std::string output_filename;
        bool split_files, fit=true;

        /*! \brief Produces additional information about
         *      current instance
         */
        std::string additional_info();

        /*! \brief Call to parse the class-specific command line
         *      parameters.
         */
        virtual void parse_params(int argc, const char *argv[]);

        /*! \brief Overloads the clone method so that 
         *      a correct pointer is returned.
         */
        virtual Computer* clone() const = 0;
    };
}

#endif
