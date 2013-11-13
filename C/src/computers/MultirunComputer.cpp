#ifndef turb_MultirunComputer_cpp
#define turb_MultirunComputer_cpp

#include "exceptions.hpp"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <iostream>

namespace turb {

    template <typename T>
    void MultirunComputer<T>::compute()
    {
        for (size_t i = 0; i < runs; ++i) {
            std::string current_filename = output_filename;
            if (split_files) {
                std::ostringstream output_number;
                output_number.width(log(runs)/log(10) + 1);
                output_number << std::setfill('0') << i;

                current_filename += output_number.str();
            }

            std::ofstream output(current_filename, std::ios::app);

            T* instance = static_cast<T*>(clone());
            try {
                instance->compute_single(&output);
                output.close();
            } catch (const RemoveOutput &e) {
                output.close();
                remove(current_filename.c_str());
            }

            delete instance;
        }
    }
}

#endif
