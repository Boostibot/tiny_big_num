#include <iostream>
#include <string>
#include <limits>

#include "benchmark.h"
#include "big_int.h"
#include "big_int_ops.h"
#include "tests.h"

int main()
{
    run_tests();
    std::cout << "\nOK" << std::endl;
}
