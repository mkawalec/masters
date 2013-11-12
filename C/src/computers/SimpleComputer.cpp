#include "helpers.hpp"
#include "computers/SimpleComputer.hpp"
#include "Integrator.hpp"
#include "helpers.hpp"

#include <fstream>


namespace turb {

    SimpleComputer::SimpleComputer()
    {
        name = "simple";
        class_name = "Computer";
        description = "Runs Integrator a certain number of times "
                      "and serializes the results every few frames";
        serializer = NULL;

        Computer::available.push_back(this);
    }

    void SimpleComputer::compute()
    {
        integrator = new Integrator(samples, dt, domain_size);
        set_constants();

        std::ofstream output(output_filename);
        for (size_t i = 0; i * dt < end_time; ++i) {
           integrator->apply_step();
           if (i%print_every == 0)
               serializer->serialize(integrator, &output,
                       i * dt);
        }

        output.close();
        delete integrator;
    }

    SimpleComputer *simple_computer_instance = new SimpleComputer();
}

