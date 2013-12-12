#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "test_integrators.hpp"
#include <fstream>
#include <string>
#include <sstream>
#include <cfloat>
#include <cmath>
#include <stdio.h>

using namespace boost::unit_test;
using namespace turb;
using namespace std;

const double dt = 0.001;
const string filename = "test_output";

template <typename T>
void setup_test(T *instance)
{
    instance->initialize();

    ofstream output;
    output.open(filename);

    for (size_t i = 0; i < 100; ++i) {
        instance->apply_step();
        instance->serialize(&output, dt * i);
    }
    output.close();
}

/** Checks if in successive iterations the value 
 *  of a 0-th component of u decays, as it should,
 *  exponentially with a = 1 + e
 */
BOOST_AUTO_TEST_CASE(test2)
{
    // Setup
    TestDecayIntegrator tested(7, dt);
    setup_test(&tested);

    ifstream input;
    input.open(filename);
    double value, prev_value = DBL_MAX;
    string temp1, temp2;
    stringstream temp;
    BOOST_CHECK(false);

    for (size_t i = 0; i < 100; ++i) {
        input >> temp1 >> temp2;

        temp.clear();
        temp.str(temp2);
        temp >> value;

        if (prev_value != DBL_MAX) 
            BOOST_CHECK(abs(value - (1 - (1 - tested.e)*dt/2)*prev_value) < 1e-05);

        prev_value = value;
    }

    remove(filename.c_str());
} 

/** Checks if the L2 norm is stable over long periods,
 *  does NOT check correctness in any way
 */
BOOST_AUTO_TEST_CASE(test3)
{
    TestStabilityIntegrator tested(7, dt);
    setup_test(&tested);

    ifstream input;
    input.open(filename);
    double value1, value2, prev_value1 = DBL_MAX, prev_value2 = DBL_MAX;
    string temp1, temp2, temp3;
    stringstream temp;

    for (size_t i = 0; i < 100; ++i) {
        input >> temp1 >> temp2 >> temp3;

        temp.clear();
        temp.str(temp2);
        temp >> value1;

        temp.clear();
        temp.str(temp3);
        temp >> value2;

        if (prev_value1 != DBL_MAX) {
            BOOST_CHECK(abs(value1 - prev_value1) < 1e-05);
            BOOST_CHECK(abs(value2 - prev_value2) < 1e-05);
        }

        prev_value1 = value1;
        prev_value2 = value2;
    }

    remove(filename.c_str());
}

/** Checks if the functions get combined to create a correct
 *  function, which mostly checks the correctness of the normalization
 *  and of Fourier transforms.
 */
// TODO: Add a templated setup function
BOOST_AUTO_TEST_CASE(test4)
{
    TestMultIntegrator tested(7, dt);
    tested.initialize();

    ofstream output;
    output.open(filename);

    tested.apply_step();
    tested.serialize(&output);
    output.close();

    // The actual testing
    ifstream input;
    input.open(filename);
    double value1, value2;
    string temp1, temp2, temp3;
    size_t iteration;
    stringstream temp;

    for (size_t i = 0; i < pow(2, 7); ++i) {
        input >> temp1 >> temp2 >> temp3;

        temp.clear();
        temp.str(temp1);
        temp >> iteration;

        temp.clear();
        temp.str(temp2);
        temp >> value1;

        temp.clear();
        temp.str(temp3);
        temp >> value2;

        double x = iteration * tested.domain_size / pow(2, 7);
        double u = cos(x) * cos(2 * x);
        double v = cos(2 * x);

        BOOST_CHECK(abs(u - value1) < 1e-05);
        BOOST_CHECK(abs(v - value2) < 1e-05);
    }

    remove(filename.c_str());
}

BOOST_AUTO_TEST_CASE(step)
{
    TestNonLinear tested(7, dt, 24 * M_PI);
    tested.initialize();

    ofstream output;
    output.open("step_output");

    tested.serialize(&output, dt);
    tested.apply_step();
    tested.serialize(&output, dt);
    output.close();

    BOOST_CHECK(true);
}
    
    

test_suite *init_unit_test_suite(int, char *[])
{
    framework::master_test_suite().p_name.value = "PUMAS unit test";

    return 0;
}
