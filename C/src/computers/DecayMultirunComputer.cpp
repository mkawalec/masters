#include "computers/DecayMultirunComputer.hpp"
#include "Computer.hpp"

#include <fstream>
#include <string>
#include <iostream>

namespace turb {

    DecayMultirunComputer::DecayMultirunComputer()
    {
        name = "decay_mult";
        class_name = "Computer";
        description = "Runs the integrator loop many times"
            " and outputs the times to decay for each run";
        serializer = NULL;
        suggested_serializer = "generic";
        
        decay_threshold = 5.0;

        Computer::available.push_back(this);
    }

    void DecayMultirunComputer::compute_single(std::ofstream *output)
    {
        set_serializer();
        integrator = new Integrator(samples, dt, domain_size);
        set_constants();

        for (size_t i = 0; i * dt < end_time; ++i) {
            integrator->apply_step();

            if (i%print_every == 0) {
                integrator->forward_transform();
                double norm = l2_norm(integrator->u, integrator->size_real);

                if (norm < decay_threshold) {
                    std::string current_time = std::to_string(i * dt);
                    serializer->serialize(integrator, output, &current_time);

                    std::cerr << "Finished run at time " << i * dt << 
                        " with norm " << norm << std::endl;
                    return;
                }
            }
        }

        delete integrator;
    }

    DecayMultirunComputer *decay_mult_instance = new DecayMultirunComputer();
}


               
