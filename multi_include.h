#pragma once

//Defaults:
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
