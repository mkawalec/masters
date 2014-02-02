#include "searchers/CMASimple.hpp"
#include "Searcher.hpp"
#include "Integrator.hpp"
#include "helpers.hpp"

#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <armadillo>
using namespace arma;


namespace turb {

    CMASimple::CMASimple()
    {
        name = "CMA-simple";
        class_name = "Searcher";
        description = "Searches using the Covariant Matrix "
            "Adaptation";

        compatible_integrators.push_back("paper");

        Searcher::available.push_back(this);
    }

    std::vector<double> CMASimple::run()
    {
        set_params();

        std::vector<std::pair<double, int> > fitness;
        double *result = new double[N];

        for (int i = 0; i < stop_iters; ++i) {
            fitness.clear();
            fitness.reserve(lambda);

            for (int j = 0; j < lambda; ++j) {
                *values[j] = *xmean + sigma * (*B * (*D % randn(N)));
                F(values[j]->memptr(), result);
                fitness.push_back(std::make_pair(l2_norm(result, N), j));
            }

            // Recombination of the mean
            std::sort(fitness.begin(), fitness.end());

            *xold = *xmean;
            xmean->fill(0);
            for (int j = 0; j < mu; ++j)
                *xmean += *values[fitness[j].second] * (*weights)(j);

            // Updating evolution paths
            *ps = (1 - cs) * *ps
                + sqrt(cs * (2 - cs) * mueff) * *invsqrtC * (*xmean - *xold) / sigma;
            int hsig = norm(*ps, 2) / sqrt(1 - pow(1 - cs, 2 * (i + 1) )) / chiN
                < 1.4 + 2 / (N + 1) ? 1 : 0;
            *pc = (1 - cc) * *pc
                + hsig * sqrt(cc * (2 - cc) * mueff) * (*xmean - *xold) / sigma;

            // Adapting covariance matrix
            mat arrtmp(N, mu, fill::zeros);
            for (int j = 0; j < mu; ++j)
                arrtmp.col(j) = 1 / sigma * *values[fitness[j].second] - *xold;

            *C = (1 - c1 - cmu) * *C
                + c1 * (*pc * pc->t()
                        + (1 - hsig) * cc * (2 - cc) * *C)
                + cmu * arrtmp * diagmat(*weights) * arrtmp.t();

            // Adapt the step size
            sigma = sigma * exp((cs / damps) * (norm(*ps, 2) / chiN - 1));

            // Diagonalization of C
            if ((i + 1) - eigenval > 1 / (c1 + cmu) / N / 10) {
                eigenval = i + 1;
                *C = symmatu(*C);
                eig_sym(*D, *B, *C);
                *D = sqrt(sign(*D) % *D);
                for (int j = 0; (unsigned)j < B->n_cols; ++j)
                    B->col(j) /= norm(B->col(j), 2);

                *invsqrtC = *B * diagmat(pow(*D, -1)) * B->t();
            }

            if (fabs(fitness[0].first) > 1e7 || isnan(fitness[0].first)) {
                delete[] result;
                throw NoResult();
            } else if (fitness[0].first < stop_fitness) {
                std::cout << "Found " << fitness[0].first << endl;

                delete[] result;
                return std::vector<double>(xmean->memptr(), xmean->memptr() + N);
            }

            if (fitness[floor(0.7 * lambda)].first - fitness[0].first < 0.1)
                sigma *= exp(0.2 * cs / damps);

            if (i%1 == 0) {
                for (int j = 0; j < 5; ++j)
                    std::cout << fitness[j].first << " ";
                std::cout << std::endl;
            }

        }
        std::cout << "finished" << std::endl;

        delete[] result;
        throw NoResult();
    }

    void CMASimple::set_params()
    {
        int i = 1;
        // Setting the weighted recombination
        for (mat::col_iterator iter = weights->begin_col(0);
             iter < weights->end_col(0); ++iter, ++i)
            *iter = (log(mu + 0.5) - log(i)) / log(10);

        // Setting the initial point
        i = 0;
        for (mat::col_iterator iter = xmean->begin_col(0);
             iter < xmean->end_col(0); ++iter, ++i)
            *iter = f[i];

        *weights /= accu(*weights);
        mueff = pow(accu(*weights), 2) / accu(*weights % *weights);

        // Adaptation parameters
        cc =  (4 + mueff / N) / (N + 4 + 2 * mueff / N);
        cs = (mueff + 2) / (N + mueff + 5);
        c1 = 2 / (pow(N + 1.3, 2) + mueff);
        cmu = std::min(1 - c1, 2 * (mueff - 2 + 1/mueff) / (pow(N + 2, 2) + mueff));
        damps = 1 + 2 * std::max(0.0, sqrt((mueff - 1) / (N + 1)) - 1) + cs;

        // Initialize dynamic strategy parameters
        pc->fill(0.0);
        ps->fill(0.0);
        B->fill(fill::eye);
        D->fill(1.0);
        *C = *B * diagmat(*D % *D) * B->t();
        *invsqrtC = *B * diagmat(pow(*D, -1)) * B->t();
        eigenval = 0;
        chiN = pow(N, 0.5) * (1 - 1 / (4 * N) + 1 / (21 * pow(N, 2)));
    }

    void CMASimple::allocate(Integrator *integrator)
    {
        this->integrator = integrator;

        // Parameter setting
        N = 2 * integrator->size_real;
        lambda = 100 + floor(3 * log(N)/log(10));
        mu = lambda / 2;
        stop_iters = 1e3 * N;

        weights = new vec(mu);
        xmean   = new vec(N);
        xold    = new vec(N);
        pc      = new vec(N);
        ps      = new vec(N);
        D       = new vec(N);
        B       = new mat(N, N);
        C       = new mat(N, N);
        invsqrtC= new mat(N, N);

        values = new vec *[lambda];
        for (int i = 0; i < lambda; ++i)
            values[i] = new vec(N);

        SimpleSearcher::allocate(integrator);
    }

    CMASimple::~CMASimple()
    {
        delete weights;
        delete xmean;
        delete xold;
        delete pc;
        delete ps;
        delete B;
        delete D;
        delete C;
        delete invsqrtC;

        for (int i = 0; i < lambda; ++i)
            delete values[i];
        //delete values;
    }


    CMASimple *cma_instance = new CMASimple();
}
