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

    class Computer : public Base<Computer> {
    protected:
        Integrator *integrator;
        virtual void compute() = 0;
        void set_constants();

        std::shared_ptr<po::options_description> options = NULL;

    public:
        typedef Base<Computer> super;

        std::thread run() { return std::thread(&Computer::compute, this); }
        Serializer *serializer;

        size_t print_every, samples, runs;
        double end_time, dt, e, a, b, D, 
               R, domain_size, threshold;
        std::string output_filename;
        bool split_files, fit;

        std::string suggested_serializer;
        std::string additional_info();
        void set_serializer();

        virtual void parse_params(int argc, const char *argv[]);
        virtual void set_options() { };

        virtual Computer* clone() const = 0;
    };
}

#endif