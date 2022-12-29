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


    {
        using namespace xigoi2;
        //std::ostringstream normal_stream;
        //std::ostringstream weird_stream;

        Integer int1 = factorial(20);
        Integer int2 = factorial(10);

        Integer diff_normal = int1 - int2;
        Digits carry;

        Integer diff_weird1;
        naf::add_or_sub(&diff_weird1.digits, &carry, int1.digits, int2.digits, false);

        Integer diff_weird2 = int1;
        naf::incr_or_decr(&carry, &diff_weird2.digits, int2.digits, false);


        std::cout << "normal: " << diff_normal << std::endl;
        std::cout << "weird1: " << diff_weird1 << std::endl;
        std::cout << "weird2: " << diff_weird2 << std::endl;
    }


    original_stream << xigoi::factorial(100);
    std::cout << "original:      " << original_stream.str() << std::endl;

    alloc_reduced_stream << xigoi2::factorial(100);
    std::cout << "alloc_reduced: " << alloc_reduced_stream.str() << std::endl;
    assert(original_stream.str() == alloc_reduced_stream.str());


    auto original_iters = count_iters(20000, 600, [&]{
        original_stream << xigoi::factorial(100);
    });

    auto alloc_reduced_iters = count_iters(20000, 600, [&]{
        alloc_reduced_stream << xigoi2::factorial(100);
    });

    std::cout << "original iters:      " << original_iters << std::endl;
    std::cout << "alloc_reduced iters: " << alloc_reduced_iters << std::endl;
    std::cout << "red / ori ratio:     " << (double) alloc_reduced_iters / (double) original_iters << std::endl;

    //test::run_tests();
    std::cout << "\nOK" << std::endl;
}
