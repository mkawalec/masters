#include "serializers/NormSerializer.hpp"
#include "helpers.hpp"

#include <fftw3.h>
#include <fstream>
#include <list>

namespace turb {

    NormSerializer::NormSerializer()
    {
        name = "norm";
        description = "Prints current time and then norm "
                      "real(u) and norm real(v)";
        Serializer::available.push_back(this);
    }

    void NormSerializer::serialize(Integrator *instance, 
            std::ofstream *output, double time)
    {
        instance->forward_transform();

        *output << time << " " << l2_norm(instance->u,
                instance->size_real) << " " <<
            l2_norm(instance->v, instance->size_real) <<
            std::endl;
    }

    NormSerializer norm_serializer_instance();
}

        
