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
    double dt = 0.001;
    string filename = "test_output";

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
    double dt = 0.001;
    string filename = "test_output";

    TestStabilityIntegrator tested(7, dt);
    tested.initialize();

}
    
    

test_suite *init_unit_test_suite(int, char *[])
{
    framework::master_test_suite().p_name.value = "PUMAS unit test";

    return 0;
}
