#include "Integrator.hpp"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

namespace turb {

    void Integrator::apply_step()
    {
        compute_nonlinear();
        compute_linear();
    }

    void Integrator::parse_params(int argc, char *argv[])
    {
        set_options();
        if (options) {
            po::variables_map vm;
            po::store(po::command_line_parser(argc, argv).
                options(*options).allow_unregistered().run(),
                vm);
            po::notify(vm);
        }
    }

    std::string Integrator::additional_info()
    {
        std::string message;

        if (options) {
            std::stringstream params_desc;
            params_desc << *options;
            message += params_desc.str();
        }
        return message;

    }

    void Integrator::setup_searcher()
    {
        if (searcher) delete searcher;
        if (selected_searcher != "no_searcher") {
            searcher = Searcher::choose(selected_searcher)->clone();
            searcher->allocate(this);
        }
    }

    void Integrator::setup_searcher(int argc, char *argv[])
    {
        setup_searcher();

        if (search)
            searcher->parse_params(argc, argv);
    }

    std::vector<double> Integrator::get_norms(std::vector<double> coords)
    {
        std::vector<double> norms;
        norms.push_back(l2_norm(coords.begin(),
                    coords.begin() + coords.size() / 2));
        norms.push_back(l2_norm(coords.begin() + coords.size() / 2,
                    coords.end()));
        return norms;
    }
}
