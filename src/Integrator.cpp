#include "Integrator.hpp"

namespace turb {

        void Integrator::apply_step()
        {
            compute_nonlinear();
            compute_linear();
        }
}
