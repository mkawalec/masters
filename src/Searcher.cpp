#include "Searcher.hpp"

#include <vector>
#include <string>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

namespace turb {

    void Searcher::set_options()
    {
        options.reset(new po::options_description(name + " options"));
        options->add_options()
            ("check-stationary", po::value<std::string>(&check_filename)->default_value(""),
             "If specified only a stationary file check will be ran");
    }

    void Searcher::parse_params(int argc, char *argv[])
    {
        set_options();
        if (options) {
            po::variables_map vm;
            po::store(po::command_line_parser(argc, argv).
                options(*options).allow_unregistered().run(),
                vm);
            po::notify(vm);
        }

        check_verify();
    }

    std::string Searcher::additional_info()
    {
        std::string message;

        if (options) {
            std::stringstream params_desc;
            params_desc << *options;
            message += params_desc.str();
        }
        return message;
    }
}
