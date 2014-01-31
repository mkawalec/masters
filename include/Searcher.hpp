#ifndef turb_Searcher_h
#define turb_Searcher_h

#include "Integrator.hpp"
#include "Base.hpp"

#include <vector>
#include <string>

namespace turb {

    class Integrator;

    class Searcher : public Base<Searcher> {
    public:
        virtual ~Searcher() { };

        virtual std::vector<double> run() = 0;
        virtual void init() = 0;
        virtual void allocate(Integrator *integrator) = 0;

        virtual Searcher* clone() const = 0;

        // Contains a list of integrators that
        // are compatible with this searcher
        std::vector<std::string> compatible_integrators;
    };
}

#endif
