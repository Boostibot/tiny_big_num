#include <cstdint>
#include <cstddef>
#include <cassert>
#include <cstring>

//Defaults:
#ifndef TINY_NUM_MAX_UNSIGNED_TYPE
#define TINY_NUM_MAX_UNSIGNED_TYPE uint64_t
#endif 

#ifndef TINY_NUM_TYPE_DEF
#define TINY_NUM_TYPE_DEF uint64_t
#endif 

#ifdef INCLUDED_NEW_TYPE
#undef INCLUDED_NEW_TYPE
#undef INCLUDED_TINY_NUM_TYPE_DEF
#endif

//if TINY_NUM_TYPE_1 is first time included (TINY_NUM_TYPE_1_INCLUDED wasnt defined yet casuse is defined after)
// => define unlock so that the library code gets compiledw with the type
#if defined(TINY_NUM_TYPE_1) && !defined(TINY_NUM_TYPE_1_INCLUDED)
#define INCLUDED_NEW_TYPE TINY_NUM_TYPE_1

#elif defined(TINY_NUM_TYPE_2) && !defined(TINY_NUM_TYPE_2_INCLUDED)
#define INCLUDED_NEW_TYPE TINY_NUM_TYPE_2

#elif defined(TINY_NUM_TYPE_3) && !defined(TINY_NUM_TYPE_3_INCLUDED)
#define INCLUDED_NEW_TYPE TINY_NUM_TYPE_3

#elif defined(TINY_NUM_TYPE_4) && !defined(TINY_NUM_TYPE_4_INCLUDED)
#define INCLUDED_NEW_TYPE TINY_NUM_TYPE_4

#elif defined(TINY_NUM_TYPE_DEF) && !defined(TINY_NUM_TYPE_DEF_INCLUDED) && !defined(INCLUDED_NEW_TYPE)
#define INCLUDED_NEW_TYPE TINY_NUM_TYPE_DEF
#define INCLUDED_TINY_NUM_TYPE_DEF
#endif

//define as included all that were included so far
#ifdef TINY_NUM_TYPE_DEF
#define TINY_NUM_TYPE_DEF_INCLUDED
#endif

#ifdef TINY_NUM_TYPE_1
#define TINY_NUM_TYPE_1_INCLUDED
#endif

#ifdef TINY_NUM_TYPE_2
#define TINY_NUM_TYPE_2_INCLUDED
#endif

#ifdef TINY_NUM_TYPE_3
#define TINY_NUM_TYPE_3_INCLUDED
#endif

#ifdef TINY_NUM_TYPE_4
#define TINY_NUM_TYPE_4_INCLUDED
#endif




#ifndef TINY_NUM_FIRST_TIME_INCLUDED
#define TINY_NUM_FIRST_TIME_INCLUDED
enum State
{
    OK = 0,
    MATH_ERROR,
    OUT_SIZE_TOO_SMALL,
    AUX_SIZE_TOO_SMALL,
    NUMBER_TOO_BIG, //some ops require number to be only low bits
    //NUMBER_NOT_STRIPPED, //considered hard errors that should be prevented! (requires outside manipulation to get into invalid state)
    //INVALID_ALIASING,
};

struct Optims
{
    bool mul_shift;
    bool mul_consts;
    bool mul_half_bits;
    bool div_shift;
    bool div_consts;
    bool rem_optims;

    size_t max_recursion_depth;
    size_t mul_quadratic_both_below_size;
    size_t mul_quadratic_single_below_size;
    size_t pow_trivial_below_power;
};  
typedef struct Optims Optims;

Optims make_def_optims();

using umax = uint64_t; //@TODO: make switch
#endif

#ifdef INCLUDED_NEW_TYPE
typedef INCLUDED_NEW_TYPE Digit;

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


typedef struct Result
{
    State state;
    Slice slice;
} Result;

typedef struct Char_Result
{
    State state;
    Char_Slice slice;
} Result;

//@TODO: move to shared - or maybe remove shared completely
//     : add inline qualifiers to all of these
size_t required_add_out_size(size_t left_size, size_t right_size);
size_t required_sub_out_size(size_t left_size, size_t right_size);
size_t required_mul_out_size(size_t left_size, size_t right_size);
size_t required_mul_quadratic_out_size(size_t left_size, size_t right_size);
size_t required_mul_quadratic_aux_size(size_t left_size, size_t right_size);
size_t required_mul_quadratic_fused_out_size(size_t left_size, size_t right_size);
size_t required_mul_karatsuba_out_size(size_t left_size, size_t right_size);
size_t required_mul_karatsuba_aux_size(size_t left_size, size_t right_size, size_t recursion_depth);
size_t required_pow_out_size(size_t num_bit_size, umax power, size_t single_digit_bit_size);
size_t required_pow_by_squaring_aux_size(size_t num_bit_size, umax power, size_t single_digit_bit_size);
size_t required_div_quotient_size(size_t num_size, size_t den_size);
size_t required_div_remainder_size(size_t num_size, size_t den_size);
size_t required_pow_trivial_aux_size(size_t num_bit_size, umax power, size_t single_digit_bit_size);
size_t required_pow_aux_size(size_t num_bit_size, umax power, size_t single_digit_bit_size);
size_t required_root_out_size(size_t num_bit_size, size_t root, size_t single_digit_bit_size);
size_t required_root_aux_size(size_t num_bit_size, size_t root, size_t single_digit_bit_size);
size_t required_to_base_out_size(size_t num_size, size_t to_base, size_t single_digit_bit_size);
size_t required_from_base_out_size(size_t rep_size, size_t from_base, size_t single_digit_bit_size);
size_t required_factorial_size(size_t of_value);

size_t optimal_mul_aux_size(size_t left_size, size_t right_size, size_t recursion_depth);
size_t optimal_pow_aux_size(size_t num_bit_size, umax power, size_t single_digit_bit_size);
size_t optimal_root_aux_size(size_t num_bit_size, size_t root, size_t single_digit_bit_size);

int compare(CSlice left, CSlice right);
bool is_equal(CSlice left, CSlice right);

Result add_short_in_place(Slice* left, CSlice right);
Result sub_short_in_place(Slice* left, CSlice right);
Result add_short(Slice* to, CSlice left, CSlice right);
Result sub_short(Slice* to, CSlice left, CSlice right);
Result mul_short(Slice* to, CSlice left, Digit right, const Optims* optims);
Result div_short(Slice* to, CSlice left, Digit right, const Optims* optims);

Div_Result div_bit_by_bit   (Slice* quotient, Slice* remainder, CSlice num, CSlice den, const Optims* optims, bool DO_QUOTIENT = true);
Div_Result div              (Slice* quotient, Slice* remainder, CSlice num, CSlice den, const Optims* optims);
Result rem   (Slice* remainder, CSlice num, CSlice den, const Optims* optims);
Result mul_quadratic_fused   (Slice* to, CSlice left, CSlice right, const Optims* optims);
Result mul_quadratic         (Slice* to, Slice* aux, CSlice left, CSlice right, const Optims* optims);
Result mul_karatsuba         (Slice* to, Slice* aux, CSlice left, CSlice right, const Optims* optims, size_t depth = 0, bool is_run_alone = true);
Result mul                   (Slice* to, Slice* aux, CSlice left, CSlice right, const Optims* optims, size_t depth = 0);
Result pow_by_squaring       (Slice* to, Slice* aux, CSlice num, umax power, const Optims* optims);
Result pow_trivial           (Slice* to, Slice* aux, CSlice num, umax power, const Optims* optims);
Result pow                   (Slice* to, Slice* aux, CSlice num, umax power, const Optims* optims);
Result root                  (Slice* to, Slice* aux, CSlice num, size_t root, const Optims* optims);
size_t log2(CSlice left);
Result factorial(Slice* out, size_t of_value, const Optims* optims);

Result from_chars(Slice* num, Char_CSlice rep, Digit base, const Optims* optims);
Char_Result to_chars(Char_Slice* rep, Slice* temp, CSlice num, Digit base, const Optims* optims);

//@TODO fuze with the other state - figure out names generic enough to apply to both but specific enough so that they are helpful
enum Pop_Digit_Error
{
    OK,
    BASE_TOO_LARGE,
    BASE_TOO_SMALL,
    SMALL_BUFFER,
    NON_STRIPPED_NUM
};

typedef struct Pop_Digit
{
    Slice buffer;
    CSlice number; //can also point into a different buffer from buffer
    size_t stashed_base = 0;
    size_t stashed_value = 0;
    size_t stashed_power = 0;
    Power_And_Powed highest;

    Pop_Digit_Error error = OK;
    bool ok = true;
} Pop_Digit;

Pop_Digit pop_digit_init(Slice* buffer, Slice num, Digit base);
Digit pop_digit(Pop_Digit* state, const Optims* optims);

typedef struct Push_Digit
{
    Slice buffer;
    size_t num_size = 0;
    umax stashed_base = 0;
    umax stashed_digit = 0;
    umax stashed_power = 0;
    Power_And_Powed highest;
} Push_Digit;

Push_Digit push_digit_init(Slice* buffer, size_t num_size = 0);
bool push_digit(Push_Digit* state, Digit digit, Digit base, const Optims* optims, bool flush = false);
Result push_digit_result(Push_Digit* state, const Optims* optims);
#endif

//Optim switches:
#ifndef DO_OPTIM_MUL_SHIFT
#define DO_OPTIM_MUL_SHIFT false
#endif

#ifndef DO_OPTIM_MUL_CONSTS
#define DO_OPTIM_MUL_CONSTS false
#endif

#ifndef DO_OPTIM_MUL_HALF_BITS
#define DO_OPTIM_MUL_HALF_BITS false
#endif

#ifndef DO_OPTIM_DIV_SHIFT
#define DO_OPTIM_DIV_SHIFT false
#endif

#ifndef DO_OPTIM_DIV_CONSTS
#define DO_OPTIM_DIV_CONSTS false
#endif

#ifndef DO_OPTIM_REM_OPTIMS
#define DO_OPTIM_REM_OPTIMS false
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