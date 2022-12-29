#include <iostream>
#include <string>
#include <limits>
#include <sstream>

#include "xigoi.h"
#include "xigoi_optims.h"
#include "tests.h"


int main()
{
    
    std::ostringstream original_stream;
    std::ostringstream alloc_reduced_stream;


    alloc_reduced_stream << xigoi_optims::factorial(100);
    std::cout << "alloc_reduced: " << alloc_reduced_stream.str() << std::endl;

    original_stream << xigoi::factorial(100);
    std::cout << "original:      " << original_stream.str() << std::endl;

    assert(original_stream.str() == alloc_reduced_stream.str());

    //Timings (in ns) to factorial(1000):
    //original custom vec: 167232267500
    //original std    vec: 161497220300
    //reduced alloc par:   150307045800
    //reduced alloc seq:   130714175200

    auto original_time = ellapsed_time([&]{
        original_stream << xigoi::factorial(1000);
    });

    auto alloc_reduced_time = ellapsed_time([&]{
        alloc_reduced_stream << xigoi_optims::factorial(1000);
    });

    std::cout << "original time:      " << original_time << std::endl;
    std::cout << "alloc_reduced time: " << alloc_reduced_time << std::endl;
    std::cout << "ori / red ratio:     " << (double) original_time / (double) alloc_reduced_time << std::endl;

    //test::run_tests();
    std::cout << "\nOK" << std::endl;
}
