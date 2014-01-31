#include "searchers/CMASimple.hpp"
#include "Searcher.hpp"
#include "Integrator.hpp"

#include <vector>
#include <cmath>
#include <algorithm>
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

    CMASimple::~CMASimple()
    {
        delete weights;
        delete xmean;
        delete pc;
        delete ps;
        delete B;
        delete D;
        delete C;
        delete invsqrtC;
    }

    std::vector<double> CMASimple::run()
    {
        set_params();

        for (int i = 0; i < stop_iters; ++i) {

        }
    }

    void CMASimple::set_params()
    {
        int i = 1;
        // Setting the weighted recombination
        for (mat::col_iterator iter = weights->begin_col(0);
             iter < weights->end_col(0); ++iter, ++i)
            *iter = log(mu + 0.5) - log(i) / log(10);

        // Setting the initial point
        i = 0;
        for (mat::col_iterator iter = xmean->begin_col(0);
             iter < weights->end_col(0); ++iter, ++i)
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
        lambda = 4 + floor(3 * log(N)/log(10));
        mu = lambda / 2;

        weights = new mat(mu, 1);
        xmean   = new mat(N, 1);
        pc      = new mat(N, 1);
        ps      = new mat(N, 1);
        B       = new mat(N, N);
        D       = new mat(N, 1);
        C       = new mat(N, N);
        invsqrtC= new mat(N, N);

        SimpleSearcher::allocate(integrator);

    }

}
