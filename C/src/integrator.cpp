#include "Integrator.hpp"

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    if (argc < 3) {
        std::cout << "The correct way to run this program is " <<
            argv[0] << " samples_power dt end_time" << std::endl;
        return -1;
    }

    char **null_p = NULL;
    size_t samples = strtoul(argv[1], null_p, 10);
    double dt = strtod(argv[2], null_p) * 0.5;
    double end_time = strtod(argv[3], null_p);

    turb::Integrator main_structure(samples, dt);

    double current_time = 0.0;
    std::ofstream output;
    output.open("output");

    while (current_time < end_time) {
        main_structure.apply_step();
        main_structure.serialize(&output, current_time);

        current_time += 2 * dt;
    }

    output.close();
    return 0;
}
