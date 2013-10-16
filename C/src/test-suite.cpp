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

BOOST_AUTO_TEST_CASE(test1)
{
    BOOST_CHECK_EQUAL(1, 1);
}


/** Checks if in successive iterations the value 
 *  of a 0-th component of u decays, as it should,
 *  exponentially with a = 1 + e
 */
BOOST_AUTO_TEST_CASE(test2)
{
    // Setup
    TestDecayIntegrator tested(7, dt);
    tested.initialize();

    ofstream output;
    ifstream input;
    output.open(filename);

    for (size_t i = 0; i < 100; ++i) {
        tested.apply_step(); 
        tested.serialize(&output, dt * i);
    }
    output.close();

    input.open(filename);
    double value, prev_value = DBL_MAX;
    string temp1, temp2;
    stringstream temp;

    for (size_t i = 0; i < 100; ++i) {
        input >> temp1 >> temp2;

        temp.clear();
        temp.str(temp2);
        temp >> value;

        if (prev_value != DBL_MAX) 
            BOOST_CHECK(abs(value - (1 - (1 + e)*dt/2)*prev_value) < 1e-05);

        prev_value = value;
        ++i;
    }

    remove(filename.c_str());
} 

/** Checks if the L2 norm is stable over long periods,
 *  does NOT check correctness in any way
 */
BOOST_AUTO_TEST_CASE(test3)
{
    TestStabilityIntegrator tested(7, dt);
    tested.initialize();

    ofstream output;
    ifstream input;
    output.open(filename);

    for (size_t i = 0; i < 100; ++i) {
        tested.apply_step();
        tested.serialize(&output, dt * i);
    }
    output.close();

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
        ++i;
    }

    remove(filename.c_str());
}

/** Checks if the functions get combined to create a correct
 *  function, which mostly checks the correctness of the normalization
 *  and of Fourier transforms.
 */
BOOST_AUTO_TEST_CASE(test4)
{
}
    
    

test_suite *init_unit_test_suite(int, char *[])
{
    framework::master_test_suite().p_name.value = "PUMAS unit test";

    return 0;
}
