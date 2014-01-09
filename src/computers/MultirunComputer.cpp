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
        if (!split_files)
            output = new std::ofstream(output_filename);


        for (size_t i = 0; i < runs; ++i) {
            std::string current_filename = output_filename;
            if (split_files) {
                std::ostringstream output_number;
                output_number.width(log(runs)/log(10) + 1);
                output_number << std::setfill('0') << i;

                current_filename += output_number.str();
                output = new std::ofstream(current_filename);
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
            } catch (const NoResult &e) {}

            delete instance;
            if (split_files) {
                if (output->is_open()) output->close();
                delete output;
            }
        }

        print_stationary();
        if (fit) fit_it(&decay_times);

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

        fold(&stationary_pts, 0);
        if (my_rank != 0) return;

        std::vector<std::vector<double> > single_pts;

        for (int i = 0; (unsigned)i < stationary_pts.size(); ++i) {
            if (!contains(&single_pts, &stationary_pts[i]))
                    single_pts.push_back(stationary_pts[i]);
        }

        if (single_pts.size() == 0) {
            std::cerr << "No points found" << std::endl;
            return;
        }

        std::cerr << "------------------------------------" << std::endl;
        std::cerr << "Amount of stationary points is " << single_pts.size() << std::endl;

        for (size_t i = 0; i < single_pts.size(); ++i) {
            std::string current_filename = "stationary";
            std::ostringstream output_number;
            output_number.width(log(stationary_pts.size())/log(10) + 1);
            output_number << std::setfill('0') << i;

            current_filename += output_number.str();
            std::ofstream output(current_filename);

            double norm_u = l2_norm(single_pts[i].begin(),
                   single_pts[i].begin() + single_pts[i].size() / 2);
            double norm_v = l2_norm(single_pts[i].begin() +
                   single_pts[i].size() / 2, single_pts[i].end());

            std::cout << norm_u << " " << norm_v <<std::endl;
            for (size_t j = 0; j < single_pts[i].size(); ++j)
                output << single_pts[i][j] << std::endl;
            output.close();
        }
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

    template <typename T>
    void MultirunComputer<T>::add_stationary(std::vector<double> stationary)
    {
        stationary_pts.push_back(stationary);
    }
}

#endif
