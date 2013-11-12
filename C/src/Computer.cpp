#include "Computer.hpp"
#include "exceptions.hpp"

#include <list>

namespace turb {

    void Computer::set_constants()
    {
        integrator->e = e;
        integrator->a = a;
        integrator->b = b;
        integrator->D = D;
        integrator->R = R;
    }
}

