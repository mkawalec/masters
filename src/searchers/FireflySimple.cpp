#include "searchers/FireflySimple.hpp"
#include "Jacobian.hpp"
#include "Searcher.hpp"
#include "Integrator.hpp"
#include "helpers.hpp"
#include "exceptions.hpp"

#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <climits>
#include <armadillo>
using namespace arma;


namespace turb {

    FireflySimple::FireflySimple()
    {
        name = "firefly-simple";
        class_name = "Searcher";
        description = "Searcher using the Firefly Algorithm";

        compatible_integrators.push_back("paper");
        set_options();
        Searcher::available.push_back(this);
    }

    void FireflySimple::allocate(Integrator *integrator)
    {
        this->integrator = integrator;
        N = 2 * integrator->size_real;

        gamma = 0.1;

        old_points  = new mat(N, points_n);
        //old_points  = new mat(points_n, N);
        new_points  = new mat(N, points_n);
        //new_points  = new mat(points_n, N);
        fitness     = new vec(points_n);
        tmp_value   = new vec(N);

        SimpleSearcher::allocate(integrator);
    }

    FireflySimple::~FireflySimple()
    {
        delete old_points;
        delete new_points;
        delete fitness;
        delete tmp_value;
    }

    std::vector<double> FireflySimple::run()
    {
        alpha = 0.001 * problem_scale;

        double last_fitness = INT_MAX;
        int last_fitness_change = 0;

        // Generate the initial points
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> d(0.0, 1.0);

        for (int i = 0; i < points_n; ++i) {
            for (int j = 0; j < N; ++j)
                (*old_points)(j, i) = f[j] + 0.01 * (double)d(gen);

            F(old_points->colptr(i), tmp_value->memptr());
            (*fitness)(i) = norm(*tmp_value, 2);
            std::cerr << (*fitness)(i) << " ";
        }
        std::cerr << std::endl;

        std::ofstream out("Firefly-" + std::to_string(fabs(d(gen))));
        out.precision(12);

        for (int i = 0; i < stop_iters; ++i) {
            // Generate the new points
            for (int j = 0; j < points_n; ++j) {
                // Copy the old point
                for (int k = 0; k < N; ++k)
                    new_points->col(j) = old_points->col(j);

                for (int k = 0; k < points_n; ++k) {
                    if (k == j || (*fitness)(k) > (*fitness)(j)) continue;

                    // Calculate the distance vector
                    *tmp_value = old_points->col(k) - old_points->col(j);

                    double distance = norm(*tmp_value, 2);

                    // Set the new coordinate
                    new_points->col(j) += beta *
                        exp(-gamma * pow(distance, 2)) * *tmp_value +
                        alpha * randn(N);
                }

                // Calculate the fitness for the point
                F(new_points->colptr(j), tmp_value->memptr());
                (*fitness)(j) = norm(*tmp_value, 2);
            }

            // Find min fitness
            int min_index = 0;
            double min_value = INT_MAX;
            for (int j = 0; j < points_n; ++j) {
                if ((*fitness)(j) < min_value) {
                    min_index = j;
                    min_value = (*fitness)(j);
                }
            }

            double norm_u = l2_norm(new_points->colptr(min_index), N/2);
            double norm_v = l2_norm(new_points->colptr(min_index) + N/2, N/2);
            out << norm_u << " " << norm_v << std::endl;
            if (i%100 == 0) std::cerr << min_index << " " << min_value << std::endl;


            if (min_value < threshold) {
                out.close();
                std::cerr << "Found!" << std::endl;
                return std::vector<double>(new_points->colptr(min_index),
                                          new_points->colptr(min_index) + N);
            } else if (min_value > 1e5) {
                out.close();
                std::cerr << "Out of bounds" << std::endl;
                throw NoResult();
            }

            if (fabs(last_fitness - min_value) > 0.01) {
                last_fitness = min_value;
                last_fitness_change = i;
            } else if (i - last_fitness_change > 200) {
                std::cerr << "Flat, ending..." << std::endl;
                out.close();
                throw NoResult();
            }

            alpha *= delta;
            mat *tmp = old_points;
            old_points = new_points;
            new_points = tmp;
        }

        out.close();
        std::cerr << "finished here" << std::endl;
        throw NoResult();
    }

    FireflySimple *firefly_instance = new FireflySimple();
}






