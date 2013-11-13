#include "Serializer.hpp"
#include "Computer.hpp"
#include "helpers.hpp"
#include "exceptions.hpp"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <sstream>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>
namespace po = boost::program_options;


turb::Computer* initialize(int argc, char *argv[])
{
    std::string output_filename, config_filename, 
        serializer_name, computer_name;
    double end_time, dt, e, a, b, D, 
           R, domain_size, threshold;
    size_t print_every, samples, runs;
    bool split_files;

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
        ;

    po::options_description simulation_opts("Simulation options");
    simulation_opts.add_options()
        ("dt,t", po::value<double>(&dt)->default_value(0.0005),
         "timestep size")
        ("end-time", po::value<double>(&end_time)->default_value(2000),
         "time when simulation ends")
        ("samples", po::value<size_t>(&samples)->default_value(7),
         "log[base 2] of the number of samples in real space")
        ("domain-size,d", po::value<double>(&domain_size)->default_value(24 * M_PI),
         "size of integration domain")
        ("threshold,t", po::value<double>(&threshold)->default_value(300),
         "threshold value. Action depends on Computer used")
        ("runs,r", po::value<size_t>(&runs)->default_value(2000),
         "total number of runs, used with multirun computers")
        ("e,e", po::value<double>(&e)->default_value(-0.1),
         "an integration parameter")
        ("a,a", po::value<double>(&a)->default_value(0.125),
         "an integration parameter")
        ("b,b", po::value<double>(&b)->default_value(-0.004),
         "an integration parameter")
        ("D,D", po::value<double>(&D)->default_value(40),
         "an integration parameter")
        ("R,R", po::value<double>(&R)->default_value(1.04),
         "integration parameter corresponding to a Reynolds number")
        ;

    po::options_description cmdline_opts;
    cmdline_opts.add(generic_opts).add(file_opts).
        add(modules_opts).add(simulation_opts);

    po::options_description config_file_opts;
    config_file_opts.add(file_opts).add(modules_opts).
        add(simulation_opts);

    po::positional_options_description p;
    p.add("config", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
            options(cmdline_opts).positional(p).run(), vm);
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
        output << cmdline_opts;
        throw turb::ProgramDeathRequest(&output);
    }
    

    turb::Computer *computer = turb::Computer::choose(computer_name)->clone();
    computer->serializer = turb::Serializer::choose(serializer_name);

    // Setting the params
    computer->print_every = print_every;
    computer->end_time = end_time;
    computer->dt = dt;
    computer->samples = samples;
    computer->domain_size = domain_size;
    computer->output_filename = output_filename;
    computer->split_files = split_files;
    computer->e = e;
    computer->a = a;
    computer->b = b;
    computer->D = D;
    computer->R = R;
    computer->threshold = threshold;
    computer->runs = runs;

    return computer;
}


int main(int argc, char *argv[])
{
    // Speeds up the performance for large
    // inputs and outputs. Makes using printf and scanf
    // very unpredicatable
    std::ios_base::sync_with_stdio(false);
    turb::Computer *computer = NULL;

    try{
        computer = initialize(argc, argv);
    } catch (const turb::ProgramDeathRequest& e) {
        std::cerr << e.what() << std::endl;
        return 0;
    }

    std::thread comp = computer->run();
    comp.join();
    return 0;
}
