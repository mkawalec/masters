#include "Computer.hpp"
#include "exceptions.hpp"
#include "Serializer.hpp"

#include <list>
#include <sstream>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

namespace turb {

    std::string Computer::additional_info()
    {
        std::string message;
        message = "The default Serializer is '" + suggested_serializer + "'.\n";

        if (options) {
            std::stringstream params_desc;
            params_desc << *options;
            message += params_desc.str();
        }
        return message;

    }

    void Computer::set_serializer()
    {
        if (serializer == NULL)
            serializer = Serializer::choose(suggested_serializer);
    }

    void Computer::parse_params(int argc, const char *argv[])
    {
        set_options();
        if (options) {
            po::variables_map vm;
            po::store(po::command_line_parser(argc, argv).
                options(*options).allow_unregistered().run(),
                vm);
            po::notify(vm);
        }

        integrator->parse_params(argc, argv);
    }

}

