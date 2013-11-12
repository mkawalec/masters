#include "serializers/NormSerializer.hpp"
#include "Integrator.cpp"
#include "helpers.hpp"

#include <fftw3.h>
#include <fstream>
#include <list>
#include <iostream>

namespace turb {
    template <typename T>
    std::list<T*> Base<T>::available;

    NormSerializer::NormSerializer()
    {
        name = "norm";
        class_name = "Serializer";
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

    NormSerializer *norm_serializer_instance = new NormSerializer();
}

        
