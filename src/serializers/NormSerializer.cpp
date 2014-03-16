#include "serializers/NormSerializer.hpp"
#include "Integrator.hpp"
#include "helpers.hpp"

#include <fftw3.h>
#include <fstream>
#include <list>


namespace turb {

    NormSerializer::NormSerializer()
    {
        name = "norm";
        class_name = "Serializer";
        description = "Prints current time and then norm "
                      "real(u) and norm real(v)";
    }

    void NormSerializer::serialize(Integrator *instance,
            std::ofstream *output, void *time)
    {
        instance->forward_transform();

        *output << *static_cast<double*>(time) << " " << l2_norm(instance->u,
                instance->size_real) << " " <<
            l2_norm(instance->v, instance->size_real) <<
            std::endl;
    }

    NormSerializer *norm_serializer_instance = new NormSerializer();
}


