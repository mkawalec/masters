#ifndef turb_MultirunComputer_h
#define turb_MultirunComputer_h

#include "Computer.hpp"

#include <fstream>
#include <vector>

namespace turb {

    /*! \brief Base class for computers that need
     *      to run more than one run to get useful data.
     */
    template <typename T>
    class MultirunComputer : public Computer {
    private:
        /*! \brief Fit a decay profile of a 
         *      vector of decay times to an e^-ax function
         */
        void fit_it(std::vector<double> *decay_times);

        /*! \brief The proportion of tail of the 
         *      decay times that will be used for the fit.
         */
        double fit_part = 0.3;

    protected:
        /*! \brief Manages multiple runs, reduces their
         *      output data etc.
         */
        void compute();

    public:
        /*! \brief Implemented by derived classes to execute
         *      a single run.
         */
        virtual double compute_single(std::ofstream *output) = 0;

    };
}

#include "computers/MultirunComputer.cpp"

#endif
