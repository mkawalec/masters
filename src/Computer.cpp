#include "Computer.hpp"
#include "exceptions.hpp"
#include "Serializer.hpp"

#include <mpi.h>
#include <list>
#include <sstream>
#include <string>
#include <iomanip>
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

    void Computer::parse_params(int argc, char *argv[])
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

    void Computer::set_filename(std::string *output)
    {
        std::string prepend = prefix;

        if (output->compare(0, prepend.length(), prepend) == 0)
            return;

        // Getting process id from MPI
        int rank, process_count;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &process_count);
        std::cerr << "Hi from rank " << rank << std::endl;

        std::stringstream process_number;
        process_number.width(log(process_count - 1)/log(10) + 1);
        process_number << std::setfill('0') << rank;

        *output = prepend + process_number.str() + "-" + *output;
    }



}

