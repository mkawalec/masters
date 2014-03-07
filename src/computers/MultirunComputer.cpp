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
        stationary_pts.reserve(runs);

        set_filename(&output_filename);

        std::ofstream *output = NULL;
        if (!split_files && use_output)
            output = new std::ofstream(output_prefix + output_filename);
        if (!use_output)
            output = new std::ofstream("/dev/null");



        for (size_t i = 0; i < runs; ++i) {
            std::string current_filename = output_filename;
            if (split_files && use_output) {
                std::ostringstream output_number;
                output_number.width(log(runs)/log(10) + 1);
                output_number << std::setfill('0') << i;

                current_filename += output_number.str();
                output = new std::ofstream(output_prefix + current_filename);
            }

            // Create a new instance, run and save the
            // output values
            T* instance = static_cast<T*>(clone());
            try {
                decay_times.push_back(instance->compute_single(output, this));
                *output << std::endl << std::endl;
            } catch (const RemoveOutput &e) {
                output->close();
                remove(current_filename.c_str());
                std::cerr << "Removing the file" << std::endl;
            } catch (const NoResult &e) {}

            delete instance;
            if (split_files) {
                if (output->is_open()) output->close();
                delete output;
            }
        }

        if (find_zeros)
            print_stationary();

        if (fit)
            fit_it(&decay_times);

        if (!split_files) {
            if (output->is_open()) output->close();
            delete output;
        }
    }


    template <typename T>
    void MultirunComputer<T>::print_stationary()
    {
        int my_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

        // TODO: Remove local duplicates before folding
        std::vector<std::vector<double> > temp_vec;
        for (int i = 0; (unsigned)i < stationary_pts.size(); ++i) {
            if (!contains(&temp_vec, &stationary_pts[i])) {
                temp_vec.push_back(stationary_pts[i]);
            }
        }

        fold(&temp_vec, 0);
        if (my_rank != 0) return;

        std::vector<std::vector<double> > single_pts;
        std::cerr << "Starting with number of points: " << temp_vec.size() << std::endl;

        for (int i = 0; (unsigned)i < temp_vec.size(); ++i) {
            if (!contains(&single_pts, &temp_vec[i]))
                    single_pts.push_back(temp_vec[i]);
        }

        if (single_pts.size() == 0) {
            std::cerr << "No points found" << std::endl;
            return;
        }

        std::cerr << "------------------------------------" << std::endl;
        std::cerr << "Amount of stationary points is " << single_pts.size() << std::endl;

        for (size_t i = 0; i < single_pts.size(); ++i) {
            std::vector<double> norms = integrator->get_norms(single_pts[i]);
            for (int j = 0; j < norms.size(); ++j) {
                if (fabs(norms[j]) > 300) continue;
            }

            std::string current_filename = "stationary";
            std::ostringstream output_number;
            output_number.width(log(single_pts.size())/log(10) + 1);

            output_number << std::setfill('0') << i;

            current_filename += output_number.str();

            for (int j = 0; j < norms.size(); ++j) {
                current_filename += "-" + std::to_string(norms[j]);
                std::cerr << norms[j] << " ";
            }
            std::cerr << std::endl;

            std::ofstream output(output_prefix + current_filename);
            for (size_t j = 0; j < single_pts[i].size(); ++j)
                output << j << " " << single_pts[i][j] << std::endl;
            output.close();

        }
    }

    template <typename T>
    void MultirunComputer<T>::fit_it(std::vector<double> *decay_times)
    {
        fold(decay_times, 0);

        int my_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        if (my_rank != 0) return;

        std::cerr << "Fitting starting... " << std::flush;
        std::sort(decay_times->begin(), decay_times->end());

        // All the survival probabilities are printed to another file
        std::ofstream probs_out(output_prefix + "survival_probability");
        for (int i = 0; (unsigned)i < decay_times->size(); ++i) {
            double prob = (decay_times->size() - i) / (double)decay_times->size();
            probs_out << decay_times->at(i) << " " << prob << std::endl;
        }
        probs_out.close();

        // The actual fitting code starts here
        alglib::real_2d_array x;
        alglib::real_1d_array y, c = "[1, 0]";

        int start = (1 - fit_part) * decay_times->size(),
            points_no = decay_times->size() - start;

        double *values = new double[points_no],
               *points = new double[points_no];

        for (int i = start; (unsigned)i < decay_times->size(); ++i) {
            values[i - start] = (decay_times->size() - i) / (double) decay_times->size();
            points[i - start] = decay_times->at(i);

            std::cout << values[i-start] << " " << points[i-start] << std::endl;
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
        std::cout << c[0] << " " << c[1] << " " << rep.rmserror
                  << " " << decay_times.size() << std::endl;

        std::ofstream fit_out(output_prefix + "fit");
        fit_out << c[0] << " " << c[1] << " " << rep.rmserror
                << " " << decay_times.size() << std::endl;
        fit_out.close();

        delete[] values;
        delete[] points;

    }

    template <typename T>
    void MultirunComputer<T>::add_stationary(std::vector<double> stationary)
    {
        stationary_pts.push_back(stationary);
    }
}

#endif
