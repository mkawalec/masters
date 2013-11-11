#include "Serializer.hpp"
#include "exceptions.hpp"

#include <list>

namespace turb {
    std::list<Serializer*> Serializer::available;

    Serializer* Serializer::choose_serializer(std::string name)
    {
        for(std::list<Serializer*>::iterator it = Serializer::available.begin();
                it != Serializer::available.end(); ++it) {
            if ((*it)->name == name)
                return *it;
        }

        throw SerializerNotFound("Serializer " + name + 
                                 "is not currently available");
    }

    void Serializer::unregister(Serializer *instance)
    {
        Serializer::available.remove(instance);
    }
}
