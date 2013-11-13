#include "serializers/GenericSerializer.hpp"
#include "helpers.hpp"
#include "Integrator.hpp"

#include <fstream>
#include <string>

namespace turb {
    
    GenericSerializer::GenericSerializer()
    {
        name = "generic";
        class_name = "Serializer";
        description = "Prints a string passed as a third argument."
            " Should only be used with Computers capable of"
            " utilizing it properly";
        
        Serializer::available.push_back(this);
    }

    void GenericSerializer::serialize(Integrator *nothing,
            std::ofstream *output, void *output_data)
    {
        unused(nothing);

        *output << *static_cast<std::string*>(output_data) << std::endl;
    }

    GenericSerializer *generic_serializer_instance = 
        new GenericSerializer();
}
