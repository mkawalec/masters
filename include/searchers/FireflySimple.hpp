#ifndef turb_Firefly_Simple_hpp
#define turb_Firefly_Simple_hpp

#include "searchers/SimpleSearcher.hpp"
#include "Integrator.hpp"
#include "Jacobian.hpp"

#include <vector>
#include <armadillo>
using namespace arma;

namespace turb {
    class FireflySimple : public SimpleSearcher {
    private:
        double beta = 1,
               gamma,
               alpha,
               delta = 0.99,
               problem_scale = 4,
               threshold = 3e-1;
        int N, points_n = 1000, stop_iters = 1e4;

        mat *old_points, *new_points;
        vec *fitness, *tmp_value;

    public:
        FireflySimple();
        ~FireflySimple();

        std::vector<double> run();
        void allocate(Integrator *integrator);

        Searcher* clone() const { return new FireflySimple(*this); }
    };
}

#endif
