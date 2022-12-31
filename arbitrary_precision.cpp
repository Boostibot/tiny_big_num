#include <iostream>
#include <string>
#include <limits>
#include <sstream>

#include "tests.h"
#include "xigoi.h"


int main()
{
    size_t fact_of = 30;
    {
        using namespace test;

        jot::Allocator* def = jot::memory_globals::default_allocator();
        Vector<umax> digits = make_sized_vector<umax>(fact_of, def);
        Slice<umax> digits_s = to_slice(&digits);
        Vector<char> converted;

        auto time = ellapsed_time([&]{
            Slice<umax> output = factorial<umax>(&digits_s, fact_of, Optim_Info{});
            size_t size = required_to_base_out_size<umax>(output.size, 10);

            converted = to_base<umax>(output, 10, def);
        });

        Slice<char> converted_s = to_slice(&converted);
        std::cout.write(converted_s.data, converted_s.size);
        std::cout << "\nFACTORIAL: " << time << std::endl; //173000
    }
    {
        auto time = ellapsed_time([&]{
            std::cout << xigoi::factorial(fact_of);
        });
        //xigoi: 33541100
        //tiny     173000
        std::cout << "\nFACTORIAL: " << time << std::endl; //173000
    }
    test::run_tests();
    std::cout << "\nOK" << std::endl;
}
