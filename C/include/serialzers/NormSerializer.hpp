#include "Serializer.hpp"
#include "Computer.hpp"
#include "helpers.hpp"

#include <fftw3.h>

namespace turb {

    class NormSerializer {
    public:
        using Serializer::Serializer;

        void serialize(double time);
    };
}

        
