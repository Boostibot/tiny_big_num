#include <cstdint>
#include <cstddef>
#include <cassert>
#include <cstring>

#if !defined(TINY_NUM_FIRST_TIME_INCLUDED) && !defined(INCLUDED_NEW_TYPE)
    #define INCLUDED_NEW_TYPE size_t
    #define TINY_NUM_FIRST_TIME_INCLUDED
#endif

#if !defined(MAX_UNSIGNED_TYPE)
    #define MAX_UNSIGNED_TYPE size_t
#endif

#if defined(INCLUDED_NEW_TYPE)

//type used for each digit of numbers
typedef INCLUDED_NEW_TYPE Digit;

//maximum avalible unsigned type
typedef MAX_UNSIGNED_TYPE umax;

enum State
{
    OK = 0,                 
    OK_EXEPTIONAL,          //State where the result is well defined yet not normal
                            // when substarcting bigger number from smaller number the result is well defined 2s complement
                            //  but maybe not what the user expected 
                            // when we correctly stop pop operation (no more digits to pop)

    UNDEFINED_VALUE,        //division by zero, root of 0
    OUT_SIZE_TOO_SMALL,
    AUX_SIZE_TOO_SMALL,

    //The following two states could be merged into one but we want to give at least some hints as to what went wrong
    // (because invalid arguments code is not very helpful)
    ARGUMENTS_TOO_BIG,      //some ops require number to be only low bits
    ARGUMENTS_TOO_SMALL,    //some ops (to_chars) requires arguments (base) to be bigger than some value (2)


    //NUMBER_NOT_STRIPPED,  //considered hard errors that should be prevented! 
                            // (requires outside manipulation to get into invalid state)
    //INVALID_ALIASING,     //considered hard error
};

typedef struct Optims
{
    bool mul_shift;
    bool mul_half_bits;
    bool div_shift;

    size_t max_recursion_depth;
    size_t mul_quadratic_both_below_size;
    size_t mul_quadratic_single_below_size;
    size_t pow_trivial_below_power;
} Optims;  

Optims make_def_optims();


//Structure used for storing big numbers - just pointer and size
//Is desigend in such a way to make it trivial to create from / store in
//  custom data structures (for example std::vector)
typedef struct Slice
{
    Digit* data;
    size_t size;
} Slice;

typedef struct CSlice
{
    const Digit* data;
    size_t size;
} CSlice;

typedef struct Char_Slice
{
    char* data;
    size_t size;
} Char_Slice;

typedef struct Char_CSlice
{
    const char* data;
    size_t size;
} Char_CSlice;


//Structure returned from most operations
// contains state and resulting Slice containing the output.
// The result can only be used if the state is OK (or OVERFLOWN but the value is complemented)!
// The output is only defined inside of the `output` Slice! 
//  (ie even though we are passing overside buffer the function is only required to fill
//   the `output` portion of the buffer - the rest migh contain junk values!)
typedef struct Result
{
    State state;
    Slice output;
} Result;

typedef struct Div_Short_Result
{
    State state;
    Slice quotient;
    Digit remainder;
} Div_Short_Result;

typedef struct Div_Result
{
    State state;
    Slice quotient;
    Slice remainder;
} Div_Result;

typedef struct Char_Result
{
    State state;
    Char_Slice output;
} Char_Result;


//Algorhitms operating on single digit as one of their arguments
Result add_short(Slice* to, CSlice left, Digit right);
Result sub_short(Slice* to, CSlice left, Digit right);
Result sub_short_in_place(Slice* left, Digit right);
Result mul_short(Slice* to, CSlice left, Digit right, const Optims* optims);
Div_Short_Result div_short(Slice* to, CSlice left, Digit right, const Optims* optims);
Div_Short_Result div_short_in_place(Slice* left, Digit right, const Optims* optims);

//general purpose algorhitms
// the most general name refers to a version that should be most commonly used and
// is designed for optimal performance - however more specific cases might benefit from more specific algorhitm
Result sub_in_place(Slice* left, CSlice right);
Result add(Slice* to, CSlice left, CSlice right);
Result sub(Slice* to, CSlice left, CSlice right);

Div_Result div_bit_by_bit    (Slice* quotient, Slice* remainder, CSlice num, CSlice den, const Optims* optims);
Div_Result div               (Slice* quotient, Slice* remainder, CSlice num, CSlice den, const Optims* optims);
Result rem                   (Slice* to, CSlice num, CSlice den, const Optims* optims);
Result mul_quadratic_fused   (Slice* to, CSlice left, CSlice right, const Optims* optims);
Result mul_quadratic         (Slice* to, Slice* aux, CSlice left, CSlice right, const Optims* optims);
Result mul_karatsuba         (Slice* to, Slice* aux, CSlice left, CSlice right, const Optims* optims);
Result mul                   (Slice* to, Slice* aux, CSlice left, CSlice right, const Optims* optims);
Result pow_by_squaring       (Slice* to, Slice* aux, CSlice num, umax power, const Optims* optims);
Result pow_trivial           (Slice* to, Slice* aux, CSlice num, umax power, const Optims* optims);
Result pow                   (Slice* to, Slice* aux, CSlice num, umax power, const Optims* optims);
Result root                  (Slice* to, Slice* aux, CSlice num, umax root, const Optims* optims);
size_t log2(CSlice left);
size_t bit_size(CSlice left); //same as log2 only different name
Result factorial(Slice* out, size_t of_value, const Optims* optims);

int compare(CSlice left, CSlice right);
bool is_equal(CSlice left, CSlice right);

//Very fast minimally formated dump of numbers into a char buffer
// for more elaborate formatting / loading use Pop_Digit / Push_Digit to write your
// own functions
Char_Result to_chars(Char_Slice* rep, Slice* temp, CSlice num, Digit base, const Optims* optims);
Result from_chars(Slice* num, Char_CSlice rep, Digit base, const Optims* optims);

//Struct used for popping diggits in base off number using auxiliary buffer
// is optimized using chunking and is intended for easy writing of custom
// printing functions
struct Pop_Digit;
typedef struct Pop_Digit Pop_Digit;

State       pop_digit_can_continue(const Pop_Digit* pop);
Pop_Digit   pop_digit_init(Slice* out, CSlice num, Digit base);
Digit       pop_digit(Pop_Digit* pop, Slice* out, const Optims* optims);
Result      pop_digit_result(Pop_Digit* pop, Slice* out, const Optims* optims);

//Struct used for pushing diggits in base to number
// is optimized using chunking and is intended for easy writing of custom
// reading functions
struct Push_Digit;
typedef struct Push_Digit Push_Digit;

State       push_digit_can_continue(const Push_Digit* push);
Push_Digit  push_digit_init(Slice* out, CSlice num, Digit base);
State       push_digit(Push_Digit* push, Slice* out, Digit digit, const Optims* optims);
Result      push_digit_result(Push_Digit* push, Slice* out, const Optims* optims);

//Functions used for querying the necessary output / auxiliary size for a given algorhitm.
// In case the specific algorhitm isnt among these use the more general option
// (ie. add_short_in_place --> required_add_out_size(left_size, 1) )
size_t required_add_out_size(size_t left_size, size_t right_size);
size_t required_sub_out_size(size_t left_size, size_t right_size);
size_t required_factorial_out_size(umax of_value);

size_t required_mul_out_size(size_t left_size, size_t right_size);
size_t required_mul_quadratic_out_size(size_t left_size, size_t right_size);
size_t required_mul_quadratic_aux_size(size_t left_size, size_t right_size);
size_t required_mul_quadratic_fused_out_size(size_t left_size, size_t right_size);
size_t required_mul_karatsuba_out_size(size_t left_size, size_t right_size);
size_t required_mul_karatsuba_aux_size(size_t left_size, size_t right_size, size_t recursion_depth);

size_t required_div_quotient_size(size_t num_size, size_t den_size);
size_t required_div_remainder_size(size_t num_size, size_t den_size);

size_t required_pow_out_size(size_t num_bit_size, umax power, size_t digit_bit_size);
size_t required_pow_by_squaring_aux_size(size_t num_bit_size, umax power, size_t digit_bit_size);
size_t required_pow_trivial_aux_size(size_t num_bit_size, umax power, size_t digit_bit_size);
size_t required_pow_aux_size(size_t num_bit_size, umax power, size_t digit_bit_size);

size_t required_root_out_size(size_t num_bit_size, umax root, size_t digit_bit_size);
size_t required_root_aux_size(size_t num_bit_size, umax root, size_t digit_bit_size);

size_t required_to_chars_out_size(size_t num_size, umax to_base, size_t digit_bit_size);
size_t required_to_chars_aux_size(size_t num_size, umax to_base, size_t digit_bit_size);
size_t required_from_chars_out_size(size_t rep_size, umax from_base, size_t digit_bit_size);

size_t required_pop_digit_out_size(size_t num_size, umax to_base, size_t digit_bit_size);
size_t required_pop_digit_out_size(size_t num_size, umax to_base, size_t digit_bit_size);
size_t required_push_digit_out_size(size_t rep_size, umax from_base, size_t digit_bit_size);


//Some algorhimts can benefit from extra storage. This storage is not required but will speed up
// the execution considerably
size_t optimal_mul_aux_size(size_t left_size, size_t right_size, size_t recursion_depth);
size_t optimal_pow_aux_size(size_t num_bit_size, umax power, size_t digit_bit_size);
size_t optimal_root_aux_size(size_t num_bit_size, umax root, size_t digit_bit_size);


#endif

//Optim switches:
#ifndef DO_OPTIM_MUL_SHIFT
#define DO_OPTIM_MUL_SHIFT false
#endif

#ifndef DO_OPTIM_MUL_HALF_BITS
#define DO_OPTIM_MUL_HALF_BITS false
#endif

#ifndef DO_OPTIM_DIV_SHIFT
#define DO_OPTIM_DIV_SHIFT false
#endif

#ifndef MAX_RECURSION_DEPTH 
#define MAX_RECURSION_DEPTH 20
#endif

#ifndef MUL_QUADRATIC_SINGLE_BELOW_SIZE 
#define MUL_QUADRATIC_SINGLE_BELOW_SIZE 3
#endif

#ifndef MUL_QUADRATIC_BOTH_BELOW_SIZE 
#define MUL_QUADRATIC_BOTH_BELOW_SIZE 4
#endif

#ifndef TRIVIAL_POW_BELOW_POWER 
#define TRIVIAL_POW_BELOW_POWER 4
#endif