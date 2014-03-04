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
        suggested_serializer = "norm";

        Computer::available.push_back(this);
    }

    void SimpleComputer::compute()
    {
        set_serializer();
        integrator->clear(samples, dt, domain_size);

        set_filename(&output_filename);

        std::ofstream output(output_prefix + output_filename);
        for (size_t i = 0; i * dt < end_time; ++i) {
            double current_time = i * dt;
            integrator->apply_step();
            if (i%print_every == 0)
                serializer->serialize(integrator, &output,
                       &current_time);
        }

        output.close();
        delete integrator;
    }

    SimpleComputer *simple_computer_instance = new SimpleComputer();
}

