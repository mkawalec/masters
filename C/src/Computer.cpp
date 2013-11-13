#include "Computer.hpp"
#include "exceptions.hpp"
#include "Serializer.hpp"

#include <list>

namespace turb {

    void Computer::set_constants()
    {
        integrator->e = e;
        integrator->a = a;
        integrator->b = b;
        integrator->D = D;
        integrator->R = R;
        integrator->initialize();
    }

    std::string Computer::additional_info()
    {
        return "The default Serializer is '" + suggested_serializer + "'.";
    }

    void Computer::set_serializer()
    {
        if (serializer == NULL)
            serializer = Serializer::choose(suggested_serializer);
    }
}

