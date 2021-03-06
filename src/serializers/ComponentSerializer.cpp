#include "serializers/ComponentSerializer.hpp"
#include "helpers.hpp"

#include <fftw3.h>
#include <fstream>
#include <list>
#include <cmath>

namespace turb {

    ComponentSerializer::ComponentSerializer()
    {
        name = "cpt";
        class_name = "Serializer";
        description = "Prints all the nonzero components of "
                      "u in Fourier space";
    }

    void ComponentSerializer::serialize(Integrator *instance,
            std::ofstream *output, void *time)
    {
        unused(time);

        for (size_t i = 0; i < instance->size_complex * 2/3.0; ++i) {
            double *value = instance->c_u[i];
            *output << sqrt(pow(*value, 2) + pow(*(value + 1), 2)) << " ";
        }

        *output << std::endl;
    }
}

DECLARE_TURB_PLUGIN(ComponentSerializer);

