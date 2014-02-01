#include "computers/DecayPathComputer.hpp"
#include "Computer.hpp"
#include "exceptions.hpp"
#include "Searcher.hpp"

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

namespace turb {

    DecayPathComputer::DecayPathComputer()
    {
        name = "decay-path";
        class_name = "Computer";
        description = "Prints decay paths for fast-decaying runs."
            " Will break your program if ran with different Serializer!";
        serializer = NULL;
        suggested_serializer = "generic";

        split_files = true;

        set_options();
        Computer::available.push_back(this);
    }

    void DecayPathComputer::set_options()
    {

        options.reset(new po::options_description(name + " options"));
        options->add_options()
            ("decay-threshold", po::value<double>(&decay_threshold)->default_value(5.0),
             "Specifies a value which, when reached by u, will"
             " indicate that u had decayed")
            ("fast-threshold", po::value<double>(&fast_threshold)->default_value(50.0),
             "If a run decays before reaching t = fast_threshold"
             " it will be counted as fast-decaying")
            ("static-interval", po::value<double>(&static_interval)->default_value(1.0),
             "Tries to find a static point every static-interval"
             " of time")
            ("find-zeros", po::value<bool>(&search)->default_value(true),
             "If true, the search for a stationary point will be performed")
            ;
    }

    double DecayPathComputer::compute_single(std::ofstream *output, MultirunComputer *base)
    {
        set_serializer();
        integrator->clear(samples, dt, domain_size);
        integrator->setup_searcher();

        history *run_history =
            new history[(size_t)(fast_threshold / dt / (double) print_every)];

        for (size_t i = 0; i * dt < fast_threshold; ++i) {
            integrator->apply_step();

            if (i%print_every == 0) {
                integrator->forward_transform();
                double norm_u = l2_norm(integrator->u, integrator->size_real);
                double norm_v = l2_norm(integrator->v, integrator->size_real);
                size_t index = i / (double) print_every;

                run_history[index].time = i * dt;
                run_history[index].u    = norm_u;
                run_history[index].v    = norm_v;

                if (integrator->tau != NULL) {
                    double norm_tau = l2_norm(integrator->tau,
                            integrator->size_real);
                    run_history[index].tau    = norm_tau;
                }

                if (norm_u < decay_threshold || (i + 2 * print_every) * dt  > fast_threshold) {
                    std::string output_data;
                    output_data.reserve(30 * index);
                    for (size_t j = 0; j < index + 1; ++j) {
                        output_data.append(std::to_string(run_history[j].time));
                        output_data += " ";
                        output_data.append(std::to_string(run_history[j].u));
                        output_data += " ";
                        output_data.append(std::to_string(run_history[j].v));

                        if (integrator->tau != NULL) {
                            output_data += " ";
                            output_data.append(std::to_string(run_history[j].tau));
                        }

                        output_data += "\n";
                    }
                    delete run_history;

                    serializer->serialize(NULL, output, &output_data);
                    std::cerr << "Run decayed fast at t = " << i * dt
                        << std::endl;
                    return i * dt;
                }
            }

            // Try to find a static point at the current position
            if (i != 0 && search && integrator->search
                    && i%(int)(static_interval / dt) == 0) {
                try {
                    integrator->searcher->init();
                    base->add_stationary(integrator->searcher->run());
                } catch(NoResult e) {}
            }
        }

        delete run_history;
        throw RemoveOutput();
    }


    DecayPathComputer *decay_path_instance = new DecayPathComputer();
}

