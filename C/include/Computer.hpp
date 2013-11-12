#ifndef turb_Computer_h
#define turb_Computer_h

#include "Serializer.hpp"
#include "Integrator.hpp"
#include "Base.hpp"

#include <thread>

namespace turb {

    class Computer : public Base<Computer> {
    protected:
        Serializer *serializer;
        Integrator *integrator;
        virtual void compute() = 0;

    public:
        std::thread run() { return std::thread(&Computer::compute, this); }

        size_t print_every;
        double end_time;
        double dt;
        size_t samples;
        size_t domain_size;
        std::string output_filename;
    };
}

#endif
