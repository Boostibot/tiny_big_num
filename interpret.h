#pragma once

#include <vector>
#include <iostream>
#include <sstream>
#include "big_int.h"

#define cast(...) (__VA_ARGS__)

namespace interpret
{
    using tiny_num::Slice;
    using cstring = const char*;

    template <typename T>
    using Vector = std::vector<T>;
    using String_Builder = std::string;
    using String = Slice<const char>;
    using Mutable_String = Slice<char>;

    using std::move;

    constexpr String slice(String_Builder const& str) noexcept 
    {
        return Slice{str.data(),str.size()};
    }

    constexpr String slice(cstring str)
    {
        size_t size = 0;
        while(str[size] != '\0')
            size ++;

        return String{str, size};
    }

    enum class Token_Error
    {
        NONE = 0,
        UNSUPPORTED_CHARACTER,
        NUM_PREFIX_STRANDED,
        NUM_INVALID_BASE_CHAR,
    };

    enum class Token_Type
    {
        SPACE,
        NUMBER,
        OPERATOR,
        BRACKET_L,
        BRACKET_R,
        NONE,
    };

    struct Token
    {
        Token_Type type;
        size_t from;
        size_t to;
        Big_Int num_value;
        char text_value; //for this simple parser we only ever use single char operators
    };

    struct Tokenizer_Exception
    {
        Token_Error type;
        String_Builder text;
        String_Builder comment;
        size_t from;
        size_t to;
    };

    bool is_space(char character) noexcept
    {
        switch(character)
        {
            case ' ': 
            case '\v': 
            case '\t': 
            case '\r': 
            case '\f': 
            case '\n': 
                return true;
            default: 
                return false;
        }
        return false;
    }

    bool is_digit(char character) noexcept
    {
        switch(character)
        {
            case '0': 
            case '1': 
            case '2': 
            case '3': 
            case '4': 
            case '5': 
            case '6': 
            case '7': 
            case '8': 
            case '9': 
                return true;
            default: 
                return false;
        }
        return false;
    }

    bool is_number(char character) noexcept
    {
        return is_digit(character)
            || character == '_';
    }

    bool is_bracket(char character) noexcept
    {
        return character == '(' || character == ')';
    }

    bool is_operator(char character) noexcept
    {
        switch(character)
        {
            case '+':
            case '-':
            case '*':
            case '/':
            case '%':
            case '^':
            case '\\':
            case '!':
                return true;
            default: 
                return false;
        }
    }

    char to_char(Digit value)
    {
        constexpr char val_to_char_table[36] = {
            '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
            'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
            'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
            'U', 'V', 'W', 'X', 'Y', 'Z'
        };

        if(value > 36)
            return -1;
        else
            return val_to_char_table[value];
    }

    Digit to_digit(char value)
    {
        switch(value)
        {
            default: return -1;
            case '0': return 0;
            case '1': return 1;
            case '2': return 2;
            case '3': return 3;
            case '4': return 4;
            case '5': return 5;
            case '6': return 6;
            case '7': return 7;
            case '8': return 8;
            case '9': return 9;
            case 'A': return 10;
            case 'a': return 10;
            case 'B': return 11;
            case 'b': return 11;
            case 'C': return 12;
            case 'c': return 12;
            case 'D': return 13;
            case 'd': return 13;
            case 'E': return 14;
            case 'e': return 14;
            case 'F': return 15;
            case 'f': return 15;
            case 'G': return 16;
            case 'g': return 16;
            case 'H': return 17;
            case 'h': return 17;
            case 'I': return 18;
            case 'i': return 18;
            case 'J': return 19;
            case 'j': return 19;
            case 'K': return 20;
            case 'k': return 20;
            case 'L': return 21;
            case 'l': return 21;
            case 'M': return 22;
            case 'm': return 22;
            case 'N': return 23;
            case 'n': return 23;
            case 'O': return 24;
            case 'o': return 24;
            case 'P': return 25;
            case 'p': return 25;
            case 'Q': return 26;
            case 'q': return 26;
            case 'R': return 27;
            case 'r': return 27;
            case 'S': return 28;
            case 's': return 28;
            case 'T': return 29;
            case 't': return 29;
            case 'U': return 30;
            case 'u': return 30;
            case 'V': return 31;
            case 'v': return 32;
            case 'W': return 32;
            case 'w': return 32;
            case 'X': return 33;
            case 'x': return 33;
            case 'Y': return 34;
            case 'y': return 34;
            case 'Z': return 35;
            case 'z': return 35;
        }
    }

    String to_string(Token_Type type)
    {
        switch(type)
        {
            case Token_Type::SPACE:     return slice("SPACE");
            case Token_Type::NUMBER:    return slice("NUMBER");
            case Token_Type::OPERATOR:  return slice("OPERATOR");
            case Token_Type::BRACKET_L: return slice("BRACKET_L");
            case Token_Type::BRACKET_R: return slice("BRACKET_R");
        };

        return slice("");
    }

    Token tokenize_number_value(String input, size_t from, size_t i, size_t base = 10)
    {
        size_t to = i;
        Big_Int_Pusher pusher;
        for(; to < input.size; to++)
        {
            char current = input[to];
            if(current == '_')
                continue;

            Digit converted = to_digit(current);
            if(converted == -1)
                break;
           
            pusher.push(converted, base);
        }

        assert(to > from && "must move!");

        Token out = {Token_Type::NUMBER, from, to};
        out.num_value = Big_Int(move(pusher));

        assert(tiny_num::is_striped_number(out.num_value.slice()));

        return out;
    }

    Token tokenize_number(String input, size_t from)
    {
        size_t i = from;
        if(i + 2 <= input.size)
        {
            if(input[i] == '0')
            {
                switch(input[i + 1])
                {
                    case 'b': return tokenize_number_value(input, from, i + 2, 2);
                    case 'o': return tokenize_number_value(input, from, i + 2, 8);
                    case 'x': return tokenize_number_value(input, from, i + 2, 16);
                }
            }
        }

        return tokenize_number_value(input, from, i, 10);
    }

    Token tokenize_operator(String input, size_t from)
    {
        return Token{Token_Type::OPERATOR, from, from + 1, 0, input[from]};
    }

    Token tokenize_space(String input, size_t from)
    {
        size_t to = from;
        for(; to < input.size; to++)
        {
            if(is_space(input[to]) == false)
                break;
        }

        assert(to > from && "must move!");

        return Token{Token_Type::SPACE, from, to};
    }

    Token tokenize_bracket(String input, size_t from)
    {
        Token_Type type = (input[from] == '(')
            ? Token_Type::BRACKET_L
            : Token_Type::BRACKET_R;

        return Token{type, from, from + 1, 0, input[from]};
    }

    Vector<Token> tokenize(String input, size_t from)
    {
        Vector<Token> output;
        for(size_t pos = from; pos < input.size; )
        {
            char current = input[pos];
            Token retrieved;
            bool is_important = true;

            if(is_space(current))
            {
                retrieved = tokenize_space(input, pos);
                is_important = false;
            }
            else if(is_bracket(current))
                retrieved = tokenize_bracket(input, pos);
            else if(is_digit(current))
                retrieved = tokenize_number(input, pos);
            else if(is_operator(current))
                retrieved = tokenize_operator(input, pos);
            else
            {
                Tokenizer_Exception exception;
                exception.type = Token_Error::UNSUPPORTED_CHARACTER;
                exception.from = pos;
                exception.to = pos + 1;
                exception.text = current;
                exception.comment = "Invalid character!";

                throw exception;
            }

            assert(pos < retrieved.to && "must move!");
            pos = retrieved.to;
            if(is_important)
                output.push_back(move(retrieved));
        }

        return output;
    }

    Vector<Token> tokenize(cstring input)
    {
        return tokenize(slice(input), 0);
    }

    String_Builder format(Big_Int const& num, size_t base = 10)
    {
        using umax = uint64_t;
        using namespace tiny_num;
        if((2 <= base && base < 36) == false)
            throw std::domain_error("format into invalid base! (must be between 2 and 35 inclusive)");

        if(num.digit_count() == 0)
            return String_Builder("0");


        size_t required = required_to_base_out_size(num.digit_count(), base, BIT_SIZE<uint64_t>);
        POD_Vector<umax> temp;
        String_Builder formatted;

        formatted.resize(required);
        temp.resize_for_overwrite(num.digit_count());

        Slice temp_slice = {temp.data(), temp.capacity()};
        tiny_num::Slice<char> output = {formatted.data(), formatted.capacity()};

        output = to_base<umax, char>(&output, &temp_slice, num.slice(), base, to_char, Optims{});
        formatted.resize(output.size);

        if(num.sign() == -1)
            formatted.insert(0, "-");

        return formatted;
    } 

    enum class Expr_Type
    {
        PRIMARY,
        BRACKET,
        PREFIX_UNARY,
        POSTFIX_UNARY,
        POW,
        MUL_DIV,
        ADD_SUB,
        PROGRAM,
    };


    struct Expr
    {
        Expr_Type type;
        char operator_text; 
        Big_Int value;
        Vector<Expr> children;    
    };

    struct Parse_Result
    {
        bool ok;
        Expr expr;
        size_t from;
        size_t to;
    };

    String to_string(Expr_Type type)
    {
        switch (type)
        {
        case Expr_Type::PRIMARY:        return slice("PRIMARY");
        case Expr_Type::BRACKET:        return slice("BRACKET");
        case Expr_Type::PREFIX_UNARY:   return slice("PREFIX_UNARY");
        case Expr_Type::POSTFIX_UNARY:  return slice("POSTFIX_UNARY");
        case Expr_Type::POW:            return slice("POW");
        case Expr_Type::MUL_DIV:        return slice("MUL_DIV");
        case Expr_Type::ADD_SUB:        return slice("ADD_SUB");
        case Expr_Type::PROGRAM:        return slice("PROGRAM");

        default: return slice("ERROR?!");
        }
    }

    Token_Type type_at(Vector<Token> const& tokens, size_t at)
    {
        if(at >= tokens.size())
            return Token_Type::NONE;

        return tokens[at].type;
    }

    bool is_op_at(Vector<Token> const& tokens, size_t at, std::initializer_list<char> op_texts)
    {
        if(type_at(tokens, at) != Token_Type::OPERATOR)
            return false;

        Token const& token = tokens[at];
        for(auto i = op_texts.begin(); i != op_texts.end(); ++i)
        {
            if(tokens[at].text_value == *i)
                return true;
        }
         
        return false;
    }

    using Parsing_Func = Parse_Result(*)(Vector<Token> const& tokens, size_t at);

    //for now we assume left associativity for everything
    Parse_Result parse_binary(Vector<Token> const& tokens, size_t at, std::initializer_list<char> op_texts, Expr_Type type, Parsing_Func above)
    {
        if(at >= tokens.size())
            return {false, {type}};

        Parse_Result left = above(tokens, at);
        if(left.ok == false)
            return left;

        while(true)
        {
            size_t op_at = left.to;
            if(is_op_at(tokens, op_at, op_texts) == false)
                break;

            size_t right_at = op_at + 1;
            Parse_Result right = parse_binary(tokens, right_at, op_texts, type, above);
            if(right.ok == false)
                break;

            Parse_Result result = {true};
            result.from = at;
            result.to = right.to;
            result.expr = Expr{type, tokens[op_at].text_value};
            result.expr.children.push_back(move(left.expr));
            result.expr.children.push_back(move(right.expr));

            left = move(result);
        }

        return left;
    }

    Parse_Result parse_expr(Vector<Token> const& tokens, size_t at);

    Parse_Result parse_primary(Vector<Token> const& tokens, size_t at)
    {
        if(type_at(tokens, at) == Token_Type::NUMBER)
        {
            Parse_Result out = {true};
            out.from = at;
            out.to = at + 1;
            out.expr = Expr{Expr_Type::PRIMARY, 0, tokens[at].num_value};
            return out;
        }

        if(type_at(tokens, at) == Token_Type::BRACKET_L)
        {
            Parse_Result inner = parse_expr(tokens, at + 1);
            if(inner.ok == false || type_at(tokens, inner.to) != Token_Type::BRACKET_R)
                return {false, {Expr_Type::BRACKET}};

            Parse_Result out = {true};
            out.from = at;
            out.to = inner.to + 1;
            out.expr = Expr{Expr_Type::BRACKET, 0};
            out.expr.children.push_back(move(inner.expr));

            return out;
        }

        return {false, {Expr_Type::PRIMARY}};
    }

    Parse_Result parse_prefix_unary(Vector<Token> const& tokens, size_t at)
    {
        if(is_op_at(tokens, at, {'+', '-', '!'}))
        {
            Parse_Result result = parse_prefix_unary(tokens, at + 1);
            if(result.ok)
            {
                Parse_Result out = {true};
                out.from = at;
                out.to = result.to;
                out.expr = Expr{Expr_Type::PREFIX_UNARY, tokens[at].text_value};
                out.expr.children.push_back(move(result.expr));

                return out;
            }
        }
        return parse_primary(tokens, at);
    }

    Parse_Result parse_pow(Vector<Token> const& tokens, size_t at)
    {
        return parse_binary(tokens, at, {'^', '\\'}, Expr_Type::MUL_DIV, parse_prefix_unary);
    }

    Parse_Result parse_mul_div(Vector<Token> const& tokens, size_t at)
    {
        return parse_binary(tokens, at, {'*', '/', '%'}, Expr_Type::MUL_DIV, parse_pow);
    }

    Parse_Result parse_add_sub(Vector<Token> const& tokens, size_t at)
    {
        return parse_binary(tokens, at, {'+', '-'}, Expr_Type::MUL_DIV, parse_mul_div);
    }

    Parse_Result parse_expr(Vector<Token> const& tokens, size_t at = 0)
    {
        return parse_add_sub(tokens, at);
    }

    std::ostream& operator <<(std::ostream& out, String str)
    {
        out.write(str.data, str.size);
        return out;
    }

    std::ostream& operator <<(std::ostream& out, Big_Int const& num)
    {
        return out << format(num);
    }

    std::ostream& operator <<(std::ostream& out, Token const& token)
    {
        out << "Token{"; 

        if(token.type == Token_Type::NUMBER)
            out << token.num_value;
        else
            out << token.text_value;

        out << ", " << to_string(token.type) << ", " << token.from << ":" << token.to;

        return out << "}";
    }

    std::ostream& operator <<(std::ostream& out, Expr_Type const& type)
    {
        return out << to_string(type);
    }

    struct Parentesized_Prefix_Print
    {
        const Expr* expr;
    };

    std::ostream& operator <<(std::ostream& out, Parentesized_Prefix_Print const& prefix)
    {
        if(prefix.expr->type == Expr_Type::PRIMARY)
            return out << prefix.expr->value;

        Parentesized_Prefix_Print left = {&prefix.expr->children[0]};
        if(prefix.expr->type == Expr_Type::PREFIX_UNARY)
            return out << '(' << prefix.expr->operator_text << ' ' << left << ')';

        if(prefix.expr->type == Expr_Type::POSTFIX_UNARY)
            return out << '(' << '`' << prefix.expr->operator_text << '`' << ' ' << left << ')';

        if(prefix.expr->type == Expr_Type::BRACKET)
            return out << left;

        Parentesized_Prefix_Print right = {&prefix.expr->children[1]};
        out << "(" << prefix.expr->operator_text << " " << left << " " << right << ")";
    }
    
    std::ostream& operator <<(std::ostream& out, Expr const& expr)
    {
        if(expr.type == Expr_Type::PRIMARY)
            return out << expr.value;

        Expr left = expr.children[0];
        if(expr.type == Expr_Type::POSTFIX_UNARY)
            return out << left << expr.operator_text;

        if(expr.type == Expr_Type::PREFIX_UNARY)
            return out << expr.operator_text << left;

        if(expr.type == Expr_Type::BRACKET)
            return out << "(" << left << ")";

        Expr right = expr.children[1];
        return out << left << " " << expr.operator_text << " " << right;
    }

    std::ostream& operator <<(std::ostream& out, Parse_Result const& result)
    {
        if(result.ok == false)
            return out << "Parse_Result{ERROR}";

        return out << "Parse_Result{OK, " << result.from << ", " << result.to << ", " << result.expr << "}";
    }


    Big_Int eval(Expr const& expr)
    {
        if(expr.type == Expr_Type::PRIMARY)
            return expr.value;

        assert(expr.children.size() >= 1);
        Expr left = expr.children[0];
        if(expr.type == Expr_Type::BRACKET)
            return eval(left);

        if(expr.type == Expr_Type::POSTFIX_UNARY || expr.type == Expr_Type::PREFIX_UNARY)
        {
            Big_Int left_num = eval(left);
            switch(expr.operator_text)
            {
                case '+': 
                    return left_num; 
                
                case '-': 
                    left_num.flip_sign();
                    return left_num;

                case '!':
                    big_int::imax single = left_num.to_single_digit();
                    return factorial(single);
            }
        }


        assert(expr.children.size() >= 2);
        Expr right = expr.children[1];

        Big_Int left_num = eval(left);
        Big_Int right_num = eval(right);

        switch(expr.operator_text)
        {
            case '+': return left_num + right_num; 
            case '-': return left_num - right_num;
            case '*': return left_num * right_num;
            case '/': return left_num / right_num;
            case '%': return left_num % right_num;

            case '^': {
                big_int::imax single = right_num.to_single_digit();
                return pow(left_num, single);
            }
            
            case '\\': {
                big_int::imax single = right_num.to_single_digit();
                return root(left_num, single);
            }
        }

        assert(false && "unreachable");
        return 0;
    }

    struct Parser_Exception
    {
        String_Builder msg;
    };


    Parse_Result parse_or_throw(Vector<Token> const& tokens, size_t at = 0)
    {
        Parse_Result parsed = parse_expr(tokens, 0);
        if(parsed.ok == false || parsed.to != tokens.size())
        {
            std::stringstream stream;
            stream << "Couldnt parse the entire program\n";
            stream << "Last parsed result: " << parsed;
            throw Parser_Exception{stream.str()};
        }

        return parsed;
    }

    Big_Int eval(String str)
    {
        auto tokens = tokenize(str, 0);
        Parse_Result parsed = parse_or_throw(tokens, 0);
        return eval(parsed.expr);
    }

    Big_Int eval(cstring str)
    {
        String string = slice(str);
        return eval(string);
    }

    void test_tokenizer()
    {
        auto tokens = tokenize("400 \\ 2 + 0xFF % 3");
        assert(tokens.size() == 7);
        assert(tokens[0].type == Token_Type::NUMBER && tokens[0].num_value == 400);
        assert(tokens[1].type == Token_Type::OPERATOR && tokens[1].text_value == '\\');
        assert(tokens[4].type == Token_Type::NUMBER && tokens[4].num_value == 255);

        try { auto tokens = tokenize("400 \\ 2 + 0x % 3"); }
        catch(Tokenizer_Exception ex)
        {
            assert(ex.type == Token_Error::UNSUPPORTED_CHARACTER);
        }

        try { auto tokens = tokenize("u 400 \\ 2 ++ % 3"); }
        catch(Tokenizer_Exception ex)
        {
            assert(ex.type == Token_Error::UNSUPPORTED_CHARACTER);
            assert(ex.from == 0);
        }
    }

    void test_parser()
    {
        auto tokens = tokenize("400 \\ 2 + 0xFF%3");
        std::stringstream stream;
        Parse_Result parsed = parse_expr(tokens);
        stream << parsed.expr;
        assert(stream.str() == "400 \\ 2 + 255 % 3");
        assert(parsed.ok == true);
        assert(parsed.to == tokens.size());
    }

    void test_eval()
    {
        //assert(Big_Int(10) - 20 == -10);
        std::cout << parse_expr(tokenize("10 - 20"));

        assert(eval("10 - 20") == -20);
        //assert(format())
        assert(eval("10*   10+ 7*10") == 170);
        assert(eval("10^9 + 123") == 1000'000'123);
        assert(eval("10^9 \\ 9") == 10);
        assert(eval("(!3^2) + --10*10 + 400 \\ 2 + 10%3") == 36 + 100 + 20 + 1);
    }

    void test_interpret()
    {
        std::cout << "TESTING INTERPRET "<< std::endl;

        std::cout << "tokenizer...";
        test_tokenizer();
        std::cout << "ok" << std::endl;

        std::cout << "parser...   ";
        test_parser();
        std::cout << "ok" << std::endl;

        std::cout << "eval...     ";
        test_eval();
        std::cout << "ok" << std::endl;

        std::cout << "=== OK ===\n" << std::endl;        
    }
}