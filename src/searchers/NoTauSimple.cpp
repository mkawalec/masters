#include "Integrator.hpp"
#include "Jacobian.hpp"
#include "searchers/SimpleSearcher.hpp"
#include "searchers/NoTauSimple.hpp"

#include <mpi.h>


namespace turb {

    std::vector<double> NoTauSimple::run()
    {
        size_t size = integrator->size_real;
        double norm;
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        for (size_t i = 0; i < iterations; ++i) {
            get_jacobian();

            // Compute the value of F at current position
            F(f, f_val1);

            if ((norm = l2_norm(f_val1, 2 * size)) < threshold) {
                std::cerr << "Found! " << i << " " << f[0] <<
                    " " << l2_norm(f + size, size) << std::endl;
                return std::vector<double>(f, f + 2 * size);
            } else if (norm > overflow) {
                throw NoResult();
            }

            try {
                gauss(f_val1, dx);
            } catch (const NoResult &e) {
                int rank;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                std::cerr << "No result, wrong at " << e.what()
                          << " at rank " << rank <<
                          " at iteration " << i << std::endl;
                throw e;
            }

            for (size_t j = 0; j < 2 * size; ++j)
                f[j] += dx[j];
        }

        std::cerr << "Nothing found! " << l2_norm(f_val1, 2 * size) << std::endl;
        throw NoResult();
    }


    void NoTauSimple::gauss(double *f, double *result)
    {
        int size = jacobian->dims().first;

        // Add F as a suffix to jacobian matrix
        for (int i = 0; i < size; ++i)
            (*jacobian)[i][size] = -f[i];

        // Bring Jacobian to the echelon form
        for (int k = 0; k < size; ++k) {
            int pivot = jacobian->max_arg(k);
            if (fuzzy_eql((*jacobian)[pivot][k], 0, 1e-10))
                throw NoResult("singular");

            jacobian->swap_lines(k, pivot);

            for (int i = k + 1; i < size; ++i) {
                jacobian_type factor = (*jacobian)[i][k] / (*jacobian)[k][k];

                for (int j = k; j < size + 1; ++j)
                    (*jacobian)[i][j] -= (*jacobian)[k][j] * factor;

                (*jacobian)[i][k] = 0;
            }
        }

        // Compute the result
        for (int i = size - 1; i >= 0; --i) {
            result[i] = (*jacobian)[i][size];

            for (int j = size - 1; j >= i; j--) {
                if (j != i)
                    result[i] -= result[j] * (*jacobian)[i][j];
                else
                    result[i] /= (*jacobian)[i][i];
            }
        }
    }

    void NoTauSimple::get_jacobian()
    {
        size_t size = 2 * integrator->size_real;
        double *tmp_f1 = (double*) fftw_malloc(sizeof(double) * size),
               *tmp_f2 = (double*) fftw_malloc(sizeof(double) * size),
               *tmp_f3 = (double*) fftw_malloc(sizeof(double) * size);
        double inv_h = 1 / h;

        for (size_t i = 0; i < size; ++i) {
            // Initialize the tmp array
            for (size_t j = 0; j < size; ++j)
                tmp_f1[j] = f[j];
            tmp_f1[i] += h;

            // Calculate the values of F
            F(tmp_f1, tmp_f2);
            F(f, tmp_f3);

            for (size_t j = 0; j < size; ++j)
                (*jacobian)[j][i] =
                    (tmp_f2[j] - tmp_f3[j]) * inv_h;
        }
        fftw_free(tmp_f1);
        fftw_free(tmp_f2);
        fftw_free(tmp_f3);
    }

    void NoTauSimple::allocate(Integrator *integrator)
    {
        this->integrator = integrator;

        jacobian = new Jacobian<jacobian_type>(2 * integrator->size_real,
                                2 * integrator->size_real + 1);
        SimpleSearcher::allocate(integrator);
    }

    NoTauSimple::NoTauSimple()
    {
        name = "no-tau-simple";
        class_name = "Searcher";
        description = "Searches for stationary points using "
            "a simple linar descent strategy";
        compatible_integrators.push_back("paper");

        Searcher::available.push_back(this);
    }

    NoTauSimple *simple_instance = new NoTauSimple;
}
