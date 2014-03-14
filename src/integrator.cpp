#include "Serializer.hpp"
#include "Computer.hpp"
#include "Searcher.hpp"

#include "helpers.hpp"
#include "exceptions.hpp"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <sstream>

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;


turb::Computer* initialize(int argc, char *argv[])
{
    std::string output_filename, config_filename,
        serializer_name, computer_name, integrator_name,
        searcher, output_prefix;
    double end_time, dt, domain_size, threshold, start_mult;
    size_t print_every, samples, runs;
    bool split_files, fit, find_zeros, use_output;

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    po::options_description generic_opts("Generic options");
    generic_opts.add_options()
        ("help,h", "produce help message")
        ;

    po::options_description file_opts("IO options");
    file_opts.add_options()
        ("config", po::value<std::string>(&config_filename),
         "config file filename")
        ("output,o", po::value<std::string>(&output_filename)->default_value("output"),
         "output file, will be overwritten")
        ("split-files", po::value<bool>(&split_files)->default_value(false),
         "put each output round in a separate file if set")
        ("print-every,p", po::value<size_t>(&print_every)->default_value(100),
         "output every p iterations")
        ("prefix", po::value<std::string>(&output_prefix)->default_value(""),
         "the prefix of all output files")
        ("use-output", po::value<bool>(&use_output)->default_value(true),
         "if specified, the standard output files will be used")
        ;

    po::options_description modules_opts("Modules options");
    modules_opts.add_options()
        ("serializer,s", po::value<std::string>(&serializer_name),
        ("name of a selected Serializer. If not specified, "
         "the default Serializer specified by a selected "
         "Computer will be used. Available are:\n\n" +
         turb::Serializer::list_available()).c_str())
        ("computer,c", po::value<std::string>(&computer_name)->default_value("simple"),
         ("name of a selected Computer. Available are:\n\n" +
         turb::Computer::list_available()).c_str())
        ("integrator,i", po::value<std::string>(&integrator_name)->default_value("paper"),
         ("name of a selected Integrator. Available are:\n\n" +
          turb::Integrator::list_available()).c_str())
        ;

    po::options_description simulation_opts("Simulation options");
    simulation_opts.add_options()
        ("dt", po::value<double>(&dt)->default_value(0.0005),
         "timestep size")
        ("end-time", po::value<double>(&end_time)->default_value(2000),
         "local simulation time when integration of a single run ends")
        ("samples", po::value<size_t>(&samples)->default_value(7),
         "log[base 2] of the number of samples in real space")
        ("domain-size,d", po::value<double>(&domain_size)->default_value(24 * M_PI),
         "size of integration domain")
        ("threshold,t", po::value<double>(&threshold)->default_value(300),
         "threshold value. Action depends on Computer used")
        ("runs,r", po::value<size_t>(&runs)->default_value(2000),
         "total number of runs, used with multirun computers")
        ("fit,f", po::value<bool>(&fit)->default_value(true),
         "if specified, the MultirunComputers will try to fit"
         " returned values of u to an exponential curve and print"
         " results to stdout")
        ("find-zeros", po::value<bool>(&find_zeros)->default_value(true),
         "If true, the search for a stationary point will be performed")
        ("searcher",
         po::value<std::string>(&searcher)->default_value("no-tau-simple"),
         ("Specifies a searcher to use. Available are:\n\n" +
         turb::Searcher::list_available()).c_str())
        ("start-mult", po::value<double>(&start_mult)->default_value(1),
         "Specifies the multiplier of starting conditions")
        ;

    po::options_description cmdline_opts;
    cmdline_opts.add(generic_opts).add(file_opts).
        add(modules_opts).add(simulation_opts);

    po::options_description config_file_opts;
    config_file_opts.add(file_opts).add(modules_opts).
        add(simulation_opts);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
            options(cmdline_opts).allow_unregistered().run(), vm);
    po::notify(vm);

    if (vm.count("config")) {
        std::ifstream config_file(config_filename);
        if (!config_file) {
            throw turb::ProgramDeathRequest(
                    "Could not open config file " + config_filename);
        }
        po::store(po::parse_config_file(config_file, config_file_opts),
                vm);
        po::notify(vm);
    }

    if (vm.count("help")) {
        std::stringstream output;

        if (my_rank == 0)
            output << "Called as: " << argv[0] << " [params]"
                << std::endl << cmdline_opts;

        throw turb::ProgramDeathRequest(&output);
    }


    turb::Computer *computer = turb::Computer::choose(computer_name)->clone();
    computer->integrator = turb::Integrator::choose(integrator_name)->clone();

    std::vector<std::string> comp_ints =
        turb::Searcher::choose(searcher)->compatible_integrators;

    bool compatible = false;
    for (int i = 0; (unsigned)i < comp_ints.size(); ++i) {
        if (comp_ints[i] == computer->integrator->name) {
            compatible = true;
            break;
        }
    }

    if (!compatible) {
        std::stringstream output;
        if (my_rank == 0)
            std::cerr << "Either no searcher is selected, "
                "or the selected searcher is not compatible "
                "with selected integrator. Disabling searching."
                << std::endl;

        computer->integrator->selected_searcher = "no_searcher";
        computer->integrator->search = false;
    } else {
        computer->integrator->selected_searcher = searcher;
    }

    computer->serializer = turb::Serializer::choose(serializer_name);


    // Setting the params
    computer->print_every = print_every;
    computer->end_time = end_time;
    computer->dt = dt;
    computer->samples = samples;
    computer->domain_size = domain_size;
    computer->output_filename = output_filename;
    computer->output_prefix = output_prefix;
    computer->find_zeros = find_zeros;
    computer->use_output = use_output;
    computer->integrator->start_mult = start_mult;
    if (!computer->split_files)
        computer->split_files = split_files;

    computer->fit = fit;
    computer->threshold = threshold;
    computer->runs = runs;
    computer->parse_params(argc, argv);

    return computer;
}

void segv_handler(int sig)
{
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (my_rank == 0) {
        void *array[10];
        size_t size;

        // get void*'s for all entries on the stack
        size = backtrace(array, 10);

        // print out all the frames to stderr
        fprintf(stderr, "Error: signal %d:\n", sig);
        backtrace_symbols_fd(array, size, STDERR_FILENO);
    }

    MPI_Finalize();
    exit(1);
}


int main(int argc, char *argv[])
{
    // Gracefully catching a segfault
    signal(SIGSEGV, segv_handler);

    MPI_Init(&argc, &argv);
    turb::Computer *computer = NULL;

    // Test the segfault
    /*int *foo = (int*)-1;
    printf("%d\n", *foo);*/

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    try {
        computer = initialize(argc, argv);
    } catch (const turb::ProgramDeathRequest& e) {
        std::cerr << e.what() << std::endl;
        MPI_Finalize();
        return 0;
    } catch (const turb::InstanceNotFound& e) {
        std::cerr << e.what() << std::endl;
        MPI_Finalize();
        return 0;
    }

    try {
        computer->run();
    } catch (const turb::ProgramDeathRequest& e) {
        if (my_rank == 0) std::cerr << e.what() << std::endl;
    }

    MPI_Finalize();
    return 0;
}
