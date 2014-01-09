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


    void Computer::fold(std::vector<std::vector<double> > *folded, int target_rank)
    {
        int my_rank, process_count,
            point_size=2*integrator->size_real;
        int max_size = point_size *
            ceil(end_time / static_interval);

        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &process_count);
        MPI_Request send_request;
        MPI_Status tmp_status;
        //std::cout << my_rank << " here1" << std::endl;

        // TODO: MPI_Barrier
        MPI_Barrier(MPI_COMM_WORLD);
        //std::cout << my_rank << " here2" << std::endl;

        if (my_rank == target_rank)
            std::cerr << "About to fold stationary points..." << std::flush;

        // Gather
        // Prepare send array
        //std::cout << my_rank << " here3" << std::endl;
        double *send_array = new double[max_size];
        for (int i = 0; (unsigned)i < folded->size(); ++i) {
            for (int j = 0; j < point_size; ++j)
                send_array[i * point_size + j] = (*folded)[i][j];
        }

        //std::cout << my_rank << " here4" << std::endl;
        MPI_Issend(send_array, folded->size() * point_size, MPI_DOUBLE,
                   target_rank, 0, MPI_COMM_WORLD, &send_request);

        //std::cout << my_rank << " here5" << std::endl;
        if (my_rank != 0) {
            MPI_Wait(&send_request, &tmp_status);
            delete[] send_array;
            return;
        }

        //std::cout << my_rank << " here6" << std::endl;
        // Receive and align
        folded->clear();

        double *recv_array = new double[max_size];

        //std::cout << my_rank << " here7" << std::endl;
        for (int i = 0; i < process_count; ++i) {
            int received_count, j = 0;
            MPI_Status recv_status;

            std::cerr << "before recv" << std::endl;
            MPI_Recv(recv_array, max_size, MPI_DOUBLE,
                    i, 0, MPI_COMM_WORLD, &recv_status);
            std::cerr << "after recv" << std::endl;
            MPI_Get_count(&recv_status, MPI_DOUBLE, &received_count);
            std::cerr << "after count " << received_count << std::endl;

            while (j < received_count) {
                std::cerr << j << " " << received_count << std::endl;
                std::vector<double> to_add(&recv_array[j], &recv_array[j] + point_size);
                folded->push_back(to_add);

                j += point_size;
            }
        }

        //std::cout << my_rank << " here8" << std::endl;
        MPI_Wait(&send_request, &tmp_status);
        //std::cout << my_rank << " here9" << std::endl;
        delete[] send_array;
        delete[] recv_array;

        if (my_rank == target_rank) std::cerr << " done" << std::endl;
    }


}

