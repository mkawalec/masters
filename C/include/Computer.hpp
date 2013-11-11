#ifndef turb_Computer_h
#define turb_Computer_h

#include "Serializer.hpp"
#include "helpers.hpp"

#include <thread>

namespace turb {

    class Computer : public Base<Computer> {
    protected:
        Serializer *serializer;

    public:
        virtual std::thread run() = 0;

        size_t print_every;
    };
}

#endif
