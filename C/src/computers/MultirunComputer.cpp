#ifndef turb_MultirunComputer_cpp
#define turb_MultirunComputer_cpp

#include "exceptions.hpp"
#include "helpers.hpp"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

// ALGLIB
#include <stdafx.h>
#include <interpolation.h>
#include <ap.h>


namespace turb {

    template <typename T>
    void MultirunComputer<T>::compute()
    {
        std::vector<double> decay_times;
        decay_times.reserve(runs);

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
                decay_times.push_back(instance->compute_single(&output));
                output << std::endl << std::endl;

                output.close();
            } catch (const RemoveOutput &e) {
                output.close();
                remove(current_filename.c_str());
            } catch (const NoResult &e) {
                output.close();
            }

            delete instance;
        }

        if (fit) fit_it(&decay_times);
    }

    template <typename T>
    void MultirunComputer<T>::fit_it(std::vector<double> *decay_times) 
    {
        std::cerr << "Fitting starting... " << std::flush;

        std::sort(decay_times->begin(), decay_times->end());

        alglib::real_2d_array x;
        alglib::real_1d_array y, c = "[0]";

        size_t start = (1 - fit_part) * decay_times->size();
        size_t points_no = decay_times->size() - start;
        double *values = new double[points_no];
        double *points = new double[points_no];

        for (size_t i = start; i < decay_times->size(); ++i) {
            values[(int)(i - start)] = (decay_times->size() - i) / (double) decay_times->size();
            points[(int)(i - start)] = decay_times->at(i);
        }

        x.setcontent(points_no, 1, points);
        y.setcontent(points_no, values);

        double epsf = 0, epsx = 0.000001, diffstep = 0.0001;
        alglib::ae_int_t maxits = 0, info;
        alglib::lsfitstate state;
        alglib::lsfitreport rep;

        alglib::lsfitcreatef(x, y, c, diffstep, state);
        alglib::lsfitsetcond(state, epsf, epsx, maxits);
        alglib::lsfitfit(state, e_x);
        alglib::lsfitresults(state, info, c, rep);
        std::cerr << "done" << std::endl;
        std::cout << c[0] << " " << rep.rmserror << std::endl;

        delete[] values;
        delete[] points;

    }
}

#endif
