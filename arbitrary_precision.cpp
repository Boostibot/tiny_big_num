#include <iostream>
#include <string>
#include <limits>
#include <sstream>

#include "xigoi.h"
#include "xigoi2.h"
#include "tests.h"


int main()
{
    
    std::ostringstream original_stream;
    std::ostringstream alloc_reduced_stream;

    original_stream << xigoi::factorial(20);
    alloc_reduced_stream << xigoi2::factorial(20);

    std::cout << "original:      " << original_stream.str() << std::endl;
    std::cout << "alloc_reduced: " << alloc_reduced_stream.str() << std::endl;
    assert(original_stream.str() == alloc_reduced_stream.str());

    auto original_iters = count_iters(2000, 200, [&]{
        original_stream << xigoi::factorial(20);
    });

    auto alloc_reduced_iters = count_iters(2000, 200, [&]{
        alloc_reduced_stream << xigoi2::factorial(20);
    });

    std::cout << "original iters:      " << original_iters << std::endl;
    std::cout << "alloc_reduced iters: " << alloc_reduced_iters << std::endl;
    std::cout << "red / ori ratio:     " << (double) alloc_reduced_iters / (double) original_iters << std::endl;

    //test::run_tests();
    std::cout << "\nOK" << std::endl;
}
