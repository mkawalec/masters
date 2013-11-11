#include "helpers.hpp"
#include "Serializer.hpp"
#include "computers/SimpleComputer.hpp"


namespace turb {

    template <typename T>
    std::list<T*> Base<T>::available;

    SimpleComputer::SimpleComputer()
    {
        name = "simple";
        description = "Runs Integrator a certain number of times "
                      "and serializes the results every few frames";
        Computer::available.push_back(this);
    }

    SimpleComputer::SimpleComputer(Serializer *serializer) :
        SimpleComputer::SimpleComputer() 
    {
        serializer = serializer;
    }

    std::thread SimpleComputer::run()
    {
    }
}

