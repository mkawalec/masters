#include <memory>

#include "Computer.hpp"

namespace turb {

    class Serializer {
    protected:
        std::shared_ptr<Computer> bound_instance;

    public:
        Serializer(const std::shared_ptr<Computer> instance);
        Serializer(const Serializer &instance);

        virtual void serialize(double time) = 0;
    };
}

