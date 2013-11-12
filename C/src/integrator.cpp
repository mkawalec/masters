#include "Serializer.hpp"
#include "Computer.hpp"
#include "helpers.hpp"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <list>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>
namespace po = boost::program_options;


void initialize()
{
    std::string output_methods_desc = turb::Base<turb::Computer>::list_available();
    std::cout << "methods: " << output_methods_desc << std::endl;
}


int main(int argc, char *argv[])
{
    // Speeds up the performance for large
    // inputs and outputs. Makes using printf and scanf
    // very unpredicatable
    std::ios_base::sync_with_stdio(false);
    initialize();
    return 0;
}
