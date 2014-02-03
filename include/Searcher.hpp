#ifndef turb_Searcher_h
#define turb_Searcher_h

#include "Integrator.hpp"
#include "Base.hpp"

#include <vector>
#include <string>
#include <boost/program_options/options_description.hpp>
namespace po = boost::program_options;

namespace turb {

    class Integrator;

    class Searcher : public Base<Searcher> {
    protected:
        std::shared_ptr<po::options_description> options = NULL;

        virtual void set_options();
        virtual void check_verify() = 0;

    public:
        virtual ~Searcher() { };

        virtual std::vector<double> run() = 0;
        virtual void init() = 0;
        virtual void allocate(Integrator *integrator) = 0;

        virtual Searcher* clone() const = 0;

        // Contains a list of integrators that
        // are compatible with this searcher
        std::vector<std::string> compatible_integrators;

        virtual void parse_params(int argc, char *argv[]);

        std::string check_filename;
        std::string additional_info();
    };
}

#endif
