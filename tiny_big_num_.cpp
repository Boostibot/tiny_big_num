#include <iostream>
#include <string>

#include "tiny_big_num.h"
//#include "tests.h"

int main()
{
    //tiny_num_test::run_tests();
}

/*
#include "tests.h"
#include "interpret.h"

int main()
{
    std::cout << 
        "=== BIG NUM INTERPRET ===\n"
        "type: expression to evaluate\n"
        "    : RUN_TESTS to run all tests\n"
        "    : AST to toggle ast printing\n"
        "    : Q to quit\n\n";

    int has_ast = -1;
    while(true)
    {
        std::string str;
        std::getline(std::cin, str);

        if(str == "RUN_TESTS")
        {
            tiny_num::test::run_tests();
            big_int::test_big_int();
            interpret::test_interpret();
            std::cout << "ALL TESTS OK" << std::endl;
            continue;
        }
        else if(str == "AST")
        {
            has_ast *= -1;
            std::cout << "AST printing toggled" << std::endl;
            continue;
        }
        else if(str == "Q")
        {
            break;
        }
        
        using namespace interpret;
        try
        {
            Big_Int result;
            if(has_ast == 1)
            {
                auto tokens = interpret::tokenize(str.data());
                auto parsed = interpret::parse_or_throw(tokens);
                std::cout << "> " << parsed.expr << std::endl;
                std::cout << "> " << interpret::Parentesized_Prefix_Print{&parsed.expr} << std::endl;
                std::cout << "> " << interpret::eval(parsed.expr) << std::endl;
            }
            else
            {
                std::cout << "> " << interpret::eval(str.data()) << std::endl;
            }
        }
        catch(Tokenizer_Exception ex)
        {
            std::cout << "Tokenizer error at column " << ex.from << std::endl;
            std::cout << "comment: \"" << ex.comment << "\"" << std::endl;
            std::cout << "text: \"" << ex.text << "\"" << std::endl;
        }
        catch(Parser_Exception ex)
        {
            std::cout << ex.msg << std::endl;
        }
        catch(std::domain_error ex)
        {
            std::cout << ex.what() << std::endl;
        }
        catch(...)
        {
            std::cout << "unspecified error occured" <<  std::endl;
        }
    }

}
*/