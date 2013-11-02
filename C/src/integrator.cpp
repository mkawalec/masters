#include "Integrator.hpp"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

int main(int argc, char *argv[])
{
    // Speeds up the performance for large
    // inputs and outputs. Makes using printf and scanf
    // very unpredicatable
    std::ios_base::sync_with_stdio(false);

    if (argc < 4) {
        std::cout << "The correct way to run this program is " <<
            argv[0] << " samples_power dt end_time runs_number" << std::endl;
        return -1;
    }

    char **null_p = NULL;
    size_t samples = strtoul(argv[1], null_p, 10);
    double dt = strtod(argv[2], null_p);
    double end_time = strtod(argv[3], null_p);
    size_t runs = strtoul(argv[4], null_p, 10);
    std::ofstream output;
    output.open("output");
    
    for (size_t j = 0; j < runs; ++j) {
        std::cout << "Run " << j << std::endl;
        turb::Integrator main_structure(samples, dt, 24 * M_PI);
        main_structure.initialize();

        double current_time = 0.0;
        size_t i = 0;

        while (current_time < end_time) {
            main_structure.apply_step();
            main_structure.forward_transform();
            
            if (l2_norm(main_structure.u, main_structure.size_real) < 5) {
                output << current_time << std::endl;
                break;
            }

            current_time += dt;
            ++i;
        }
    }

    output.close();
    return 0;
}
