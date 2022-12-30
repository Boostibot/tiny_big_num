#include <iostream>
#include <string>
#include <limits>
#include <sstream>

#include "tests.h"
#include "xigoi.h"


int main()
{

    auto time = ellapsed_time([&]{
        xigoi::factorial(1000);
    });
    //xigoi: 33541100
    //tiny     173000
    std::cout << "FACTORIAL 1000: " << time << std::endl; //173000

    test::run_tests();
    std::cout << "\nOK" << std::endl;
}
