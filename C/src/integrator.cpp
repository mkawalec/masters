#include "helpers.hpp"
#include "Serializer.hpp"
#include "Computer.hpp"
#include "computers/SimpleComputer.hpp"

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
    std::string output_methods_desc;
    std::list<turb::Serializer*>::iterator it;
    for (it = turb::Serializer::available.begin();
         it != turb::Serializer::available.end(); ++it) {
        output_methods_desc += (*it)->name + ": \t" +
            (*it)->description + "\n\n";
    }
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
