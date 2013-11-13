#include "Computer.hpp"
#include "exceptions.hpp"
#include "Serializer.hpp"

#include <list>
#include <sstream>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

namespace turb {

    void Computer::set_constants()
    {
        integrator->e = e;
        integrator->a = a;
        integrator->b = b;
        integrator->D = D;
        integrator->R = R;
        integrator->initialize();
    }

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
        std::cout << "inhere" << std::endl;
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
            options(*options).allow_unregistered().run(), 
            vm);
        po::notify(vm);
    }

}
