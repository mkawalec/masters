#include "computers/DecayPathComputer.hpp"
#include "Computer.hpp"
#include "exceptions.hpp"

#include <fstream>
#include <string>
#include <iostream>

namespace turb {

    DecayPathComputer::DecayPathComputer()
    {
        name = "decay-path";
        class_name = "Computer";
        description = "Prints decay paths for fast-decaying runs."
            "Will break your program if ran with different Serializer!";
        serializer = NULL;
        suggested_serializer = "generic";

        decay_threshold = 5.0;
        fast_threshold = 50.0;
        split_files = true;

        Computer::available.push_back(this);
    }

    double DecayPathComputer::compute_single(std::ofstream *output)
    {
        set_serializer();
        integrator = new Integrator(samples, dt, domain_size);
        set_constants();

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

                if (norm_u < decay_threshold) {
                    std::string output_data;
                    output_data.reserve(30 * index);
                    for (size_t j = 0; j < index + 1; ++j) {
                        output_data.append(std::to_string(run_history[j].time));
                        output_data += " ";
                        output_data.append(std::to_string(run_history[j].u));
                        output_data += " ";
                        output_data.append(std::to_string(run_history[j].v));
                        output_data += "\n";
                    }
                    delete run_history;

                    serializer->serialize(NULL, output, &output_data);
                    std::cerr << "Run decayed fast at t = " << i * dt 
                        << std::endl;
                    return i * dt;
                }
            }
        }

        delete run_history;
        throw RemoveOutput();
    }


    DecayPathComputer *decay_path_instance = new DecayPathComputer();
}

