#ifndef INC_PRECISION
#define INC_PRECISION

/*
 *******************************************************************************
 *
 *
 *                       Copyright (c) 2002-2021
 *                       Henrik Vestermark
 *                       Denmark
 *
 *                       All Rights Reserved
 *
 *   This source file is subject to the terms and conditions of the
 *   Henrik Vestermark Software License Agreement which restricts the manner
 *   in which it may be used.
 *   Mail: hve@hvks.com
 *
 *******************************************************************************
*/
/*
 *******************************************************************************
 *
 *
 * Module name     :   iprecision.h
 * Module ID Nbr   :
 * Description     :   Arbitrary integer precision class
 * --------------------------------------------------------------------------
 * Change Record   :
 *
 * Version	Author/Date		Description of changes
 * -------  ---------------	----------------------
 * 03.01	HVE/14-Aug-2021	Switch to new internal binary format iptype
 * 03.02	HVE/19-Oct-2021	Passed all float testing
 * 03.03	HVE/3-Nov-2021	Added support for !,&&,|| operator working on int_precision class
 * 03.04	HVE/19-Nov-2021	Fixed compiler bugs reported by GNU cmpler on Mac
 * 03.05	HVE/20-Nov-2021 A few bugs fixed and change to avoid compiler warnings
 * 03.06	HVE/21-Nov-2021 More cleaning up and improvements
 * 03.07	HVE/11-Dec-2021	Change atoip() to be a string reference
 * 03.08	HVE/19-Jan-2022 Rename karatsuba & schonhage_strassen to followed the name convention for function call.
 *							Renamed _int_precision_umul to _int_precision_umul_school and added _int_precision_umul as a common entry for multiplcation that
 *							branch out to the most optimal multiplication algorithm to call. e.g. school, karatsuba, FFT or Schohage_strassen
 * 03.09	HVE/21-Mar-2022 Added fixed size arbitrry integers by introduce the mLimit field in the int_precision class. Added the method .precision()
 * 03.10	HVE/25-Mar-2022	Fixed a bug in the _int_precision_unegate()
 * 03.11	HVE/7-Aug-2022	Handle a carry bug in the _int_precision_umul_add()
 * 03.12	HVE/24-Aug-2022	Cleaning up code
 * 03.13	HVE/26-Aug-2022	Added int_precision constructor for float and double for completeness
 * 03.14	HVE/5-Sep-2022	Restore the _int_precision_compare() to the previous pointer arguments
 *
 * End of Change Record
 * --------------------------------------------------------------------------
*/

/* define version string */
static char _VI_[] = "@(#)iprecision.h 03.14 -- Copyright (C) Henrik Vestermark";

// If _INT_PRECESION_FAST_DIV_REM is defined it will use a magnitude faster div and rem integer operation.
//#define _INT_PRECISSION_FAST_DIV_REM

#include <climits>
#include <cstdint>
#include <string>
#include <vector>
#include <complex>   // Need <complex> to support FFT functions for fast multiplications
#include <cstdlib>

// THIS is the only configuration parameter to set or change.
typedef uintmax_t iptype;	// The default size of the internal binary vector type, an unsigned 64bit. It should ALWAYS be set to the 'biggest' integer type.
							// However performance will suffer if iptype < 64bit or not set to the maxium the enviroment can handle.
const unsigned int Bitsiptype = sizeof(iptype) * 8;  // Const use throughtout the source which is the number of bits the iptype can hold.

// Definining som default base numbers
static const int BASE_2	  = 2;
static const int BASE_8   = 8;
static const int BASE_10  = 10;  // Default
static const int BASE_16  = 16;

static const int RADIX = BASE_10;			// Set internal base for the arbitrary precision. NOT USED ANYMORE

// this is only for a few instance where we still use decimal arithmetic
inline int CHAR_SIGN( char x )            { return x == '-' ? -1 : 1; }
inline unsigned char IDIGIT10( char x )   { return (unsigned char)( x - '0'); }
inline unsigned char ICHARACTER10( char x){ return (unsigned char)( x + '0'); }

#undef TEMPLIFY

//
// @class int_precision
// @author Henrik Vestermark (hve@hvks.com)
// @date	14/Aug/2021
// @brief  This is an arbitrary integer class
//
// @todo
//
// Precision class
// Also number is always strip for leading zeros
// Since we always initiate to a valid int_precision number, it will never be a empty number e.g. mNumber.size()>=1
// The least significan number is at mNumber[0], the most signidicant number is a mNumber[n-1]
// For iptype = uint64_t the number in mNumber is stored as:
//		mNumber=mNumber[0]+mNumber[1]*2^64+mNumber[2]*(2^64)^2...mNumber[n-1]*(2^64)^(n-1)
// if iptype=uint3_t the number in mNumber is stored as:
//		mNumber=mNumber[0]+mNumber[1]*2^32+mNumber[2]*(2^32)^2...mNumber[n-1]*(2^32)^(n-1)
// for short we will use the notation that the least significant part of the number is stored in mNumber[0] and will be denoted a0, the most significant of the number in mNumber[n-1] as an-1
// the radix, R will be 2^64 for iptype=uint64_t and 2^32 for iptype=uint3_t etc.
// the number in mNumber can be written as:
//		mNumber=a0*R^0+a1*R^1+a2*R^2...an-1*R^n-1
//
class int_precision
	{
	int mSign;						// Sign of the int_precision. Version 2+ only. In version 2 sign has been separated from mNumber to avoid many uncessary copies and string.substr() calls
									// mSign is either +1 or -1. For mNumber==0 then sign is always +1
	size_t mLimit;					// By default int_precision i ulimited precision but it can be limit to force a certain size of an integer. e.g. 128bit has limit=2, 512bit has limit=8 etc.
									// unlimit preciion has a limit on UINTMAX_MAX or SIZE_MAX
	std::vector<iptype> mNumber;	// The binary vector of iptype that holds the integer. Per definition the vector when the constructor is invoked will always be initialized to zero if no argument is provided.


	public:
	// Constructor
	int_precision();										// No initialization
 	int_precision( char, const size_t=SIZE_MAX);					// When initialized through a char
    int_precision( unsigned char, const size_t=SIZE_MAX);			// When initialized through a unsigned char
    int_precision( short, const size_t=SIZE_MAX);					// When initialized through an short
    int_precision( unsigned short, const size_t=SIZE_MAX);			// When initialized through an unsigned short
    int_precision( int, const size_t=SIZE_MAX);						// When initialized through an int
    int_precision( unsigned int, const size_t=SIZE_MAX);				// When initialized through an unsigned int
    int_precision( long, const size_t=SIZE_MAX);						// When initialized through an long
    int_precision( unsigned long, const size_t=SIZE_MAX);			// When initialized through an unsigned long
	int_precision( long long, const size_t=SIZE_MAX);				// When initialized through an long. Same as int64_t
	int_precision( unsigned long long, const size_t=SIZE_MAX);		// When initialized through an unsigned long. Same as uintmax_t
    int_precision( const char *, const size_t=SIZE_MAX);			// When initialized through a char string
	int_precision(const float, const size_t=SIZE_MAX);				// When initialized through a float
	int_precision(const double, const size_t=SIZE_MAX);				// When initialized through a double
	int_precision( const std::string&, const size_t=SIZE_MAX);		// When initialized through a std::string
	int_precision( const std::vector<iptype>&, const size_t=SIZE_MAX);	// When initialized through a std::vector<iptype>. Notice sign will be 1 since vector<iptype> is unsigned
	int_precision( const int_precision&, const size_t=SIZE_MAX);		// When initialized through another int_precision

	//template<class _TY> inline int_precision(_TY c);		// Not working as intended

    // Coordinate memebr functions
	std::vector<iptype> *pointer();							// Return a pointer to mNumber
	std::vector<iptype> number() const;						// Return a copy of mNumber
	std::vector<iptype> number(std::vector<iptype> &mb);	// Set mNumber and return a copy of mNumber 
	iptype index(const size_t inx)	const;
	int sign() const;			// Return current sign
	int sign(int s);			// Set and return sign
	int change_sign();			// Toggle and return sign 
	size_t size() const;		// Return number of iptype digits. iptype is the allocation unit of typicall 64bit?
	size_t precision() const;	// Return the maximum fixed integer precision the variable can hold in number of iptype units
	size_t precision(const size_t p);	// Set a new fixed integer precision or arbitrary precision
	int_precision& abs();		// Change sign to + and return number
	// Start of Bit Methods
	void setbit(size_t i);		// Set bit at bi position i
	void resetbit(size_t i);	// Reset bit  bit position i
	void flipbit(size_t i);		// Flip bit at bit position i
	bool testbit(size_t i);		// Test bit at bit position i
	size_t ctz();				// Count trailing zeros
	size_t clz();				// Count leading zeros
	bool even() const;			// Test for even number
	bool odd() const;			// Test for odd number
	bool iszero()	const;		// Test for zero and return true or false
	// End of Bit Methods
	// Conversion methods. Safer and less ambiguous than overloading implicit/explicit conversion operators
	std::string toString(const int);	//  Convert number to Decimal String with an optional base parameter
	
	// Implicit/explicit conversion operators	
	operator float() const;
	operator double() const;
#ifdef TEMPLIFY
	template<class _TY> operator _TY() const;
#else
	operator long() const; 
    operator int() const;
    operator short() const;
    operator char() const;
    operator unsigned long() const;
    operator unsigned int() const;
    operator unsigned short() const;
    operator unsigned char() const;
	operator long long() const;
	operator unsigned long long() const;
#endif

    // Essential assignment operators
    int_precision& operator=( const int_precision& );
    int_precision& operator+=( const int_precision& );
    int_precision& operator-=( const int_precision& );
    int_precision& operator*=( const int_precision& );
    int_precision& operator/=( const int_precision& );
    int_precision& operator%=( const int_precision& );
    int_precision& operator>>=( const int_precision& );
    int_precision& operator<<=( const int_precision& );
	int_precision& operator&=( const int_precision& );
    int_precision& operator|=( const int_precision& );
	int_precision& operator^=( const int_precision& );

    // Specialization
	friend std::ostream& operator<<( std::ostream& strm, const int_precision& d );
	friend std::istream& operator>>( std::istream& strm, int_precision& d );

	
    // Exception class
    class bad_int_syntax {};
    class out_of_range   {};
    class divide_by_zero {};
	};

	// Arithmetic
	template <class _Ty> inline int_precision operator+(int_precision&, const _Ty&);
	template <class _Ty> inline int_precision operator+(const _Ty&, const int_precision&);
	inline int_precision operator+(const int_precision&);  // Unary
	inline int_precision operator++(int_precision&);       // Prefix Increment
	inline int_precision operator++(int_precision&, int);  // Postfix Increment

	template <class _Ty> inline int_precision operator-(int_precision&, const _Ty&);
	template <class _Ty> inline int_precision operator-(const _Ty&, const int_precision&);
	inline int_precision operator-(const int_precision&);  // Unary
	inline int_precision operator--(int_precision&);       // Prefix Decrement
	inline int_precision operator--(int_precision&, int);  // Postfix Decrement

	template <class _Ty> inline int_precision operator*(int_precision&, const _Ty&);
	template <class _Ty> inline int_precision operator*(const _Ty&, const int_precision&);
	template <class _Ty> inline int_precision operator/(int_precision&, const _Ty&);
	template <class _Ty> inline int_precision operator/(const _Ty&, const int_precision&);
	template <class _Ty> inline int_precision operator%(int_precision&, const _Ty&);
	template <class _Ty> inline int_precision operator%(const _Ty&, const int_precision&);
	template <class _Ty> inline int_precision operator<<(int_precision&, const _Ty&);
	inline int_precision operator<<(const int_precision&, const int_precision&);
	//template <class _Ty> inline int_precision operator<<( const _Ty&, const int_precision& );  // Dont allow to avoid overloading conflict in streams library??
	template <class _Ty> inline int_precision operator >> (int_precision&, const _Ty&);
	template <class _Ty> inline int_precision operator >> (const _Ty&, const int_precision&);

	template <class _Ty> inline int_precision operator&(int_precision&, const _Ty&);
	template <class _Ty> inline int_precision operator&(const _Ty&, const int_precision&);
	inline int_precision operator|(const int_precision&, const int_precision&);
	template <class _Ty> inline int_precision operator|(int_precision&, const _Ty&);
	//template <class _Ty> inline int_precision operator|(const _Ty&, const int_precision&);

	template <class _Ty> inline int_precision operator ^(int_precision&, const _Ty&);
	template <class _Ty> inline int_precision operator ^(const _Ty&, const int_precision&);
	//inline int_precision operator|(int_precision&, const int_precision&);
	int_precision operator~(const int_precision&);        // Unary Negate

	// Boolean Comparision Operators
	template <class _Ty> inline bool operator==(int_precision&, const _Ty&);
	template <class _Ty> inline bool operator==(const _Ty&, const int_precision&);
	template <class _Ty> inline bool operator!=(int_precision&, const _Ty&);
	template <class _Ty> inline bool operator!=(const _Ty&, const int_precision&);
	template <class _Ty> inline bool operator>(int_precision&, const _Ty&);
	template <class _Ty> inline bool operator>(const _Ty&, const int_precision&);
	template <class _Ty> inline bool operator>=(int_precision&, const _Ty&);
	template <class _Ty> inline bool operator>=(const _Ty&, const int_precision&);
	template <class _Ty> inline bool operator<=(int_precision&, const _Ty&);
	template <class _Ty> inline bool operator<=(const _Ty&, const int_precision&);
	template <class _Ty> inline bool operator<(int_precision&, const _Ty&);
	template <class _Ty> inline bool operator<(const _Ty&, const int_precision&);

	inline bool operator&&(const int_precision&, const int_precision&);	// Logical AND
	inline bool operator||(const int_precision&, const int_precision&);	// Logical OR
	inline bool operator!(const int_precision&);						// Logical NOT

	// Integer Precision functions
	extern int_precision abs(const int_precision&);							// return |a|
	extern int_precision ipow(const int_precision&, const int_precision&);  // return a^b
	extern int_precision ipow_modular(const int_precision&, const int_precision&, const int_precision&); // return a^b%c
	extern bool iprime(const int_precision&);
	template <class _TY> inline _TY gcd(const _TY lhs, const _TY rhs);
	extern int_precision gcd(const int_precision&, const int_precision&);		// return greatest comon divisor
	extern int_precision lcm(const int_precision&, const int_precision&);		// return least comon multiplier

	// Core Support functions 
	double _int_precision_iptod(const int_precision *);					// Explicit conversion to double
	std::vector<iptype> _int_precision_atoip(const char *, int *);		// char * string to int_precision
	std::vector<iptype> _int_precision_atoip(const std::string&, int *);// STL String to int_precision
	std::string itostring(const int, const unsigned);
	std::string _int_precision_itoa(int_precision&, const int base = BASE_10);
	std::string _int_precision_itoa(int_precision *, const int base = BASE_10);
	std::string _int_precision_itoa(const std::vector<iptype> *, const int base = BASE_10);
	std::string _int_precision_itoa(const std::vector<iptype>&, const int base = BASE_10);
	size_t _int_precision_ctz(const iptype);							// Count trailing zero bits in iptype
	size_t _int_precision_ctz(const std::vector<iptype> &);				// Cout trailing zero bits in a vector<iptype>
	size_t _int_precision_clz(const iptype);							// Count leading zero bits in iptype
	size_t _int_precision_clz(const std::vector<iptype> &);				// Count leading zero bits in a vector<iptype>

	// Core Binary functions that works directly on vector<iptype> class and unsigned arithmetic
	std::vector<iptype> _int_precision_uadd(const std::vector<iptype> *, const std::vector<iptype> *);
	std::vector<iptype> _int_precision_uadd_short(const std::vector<iptype> *, const iptype);
	std::vector<iptype> _int_precision_usub(int *, const std::vector<iptype> *, const std::vector<iptype> *);
	std::vector<iptype> _int_precision_usub_short(int *, const std::vector<iptype> *, const iptype);
	std::vector<iptype> _int_precision_umul(const std::vector<iptype> *, const std::vector<iptype> *);
	std::vector<iptype> _int_precision_umul_school(const std::vector<iptype> *, const std::vector<iptype> *);
	std::vector<iptype> _int_precision_umul_short(const std::vector<iptype> *, const iptype);
	std::vector<iptype> _int_precision_umul64(const std::vector<iptype> *, const std::vector<iptype> *);
	std::vector<iptype> _int_precision_umul_fourier(const std::vector<iptype> *, const std::vector<iptype> *, int = 8);
	std::vector<iptype> _int_precision_umul_karatsuba(const std::vector<iptype> *, const std::vector<iptype> *);
	std::vector<iptype> _int_precision_umul_linear(const std::vector<iptype> *, const std::vector<iptype> *);
	std::vector<iptype> _int_precision_udiv(const std::vector<iptype> *, const std::vector<iptype> *);
	std::vector<iptype> _int_precision_udiv_short(iptype *, const std::vector<iptype> *, const iptype);
	std::vector<iptype> _int_precision_urem(const std::vector<iptype> *, const std::vector<iptype> *);
	std::vector<iptype> _int_precision_urem_short(const std::vector<iptype> *, const iptype);
	std::vector<iptype> _int_precision_udivrem(std::vector<iptype> *, std::vector<iptype> *, std::vector<iptype> *);
	std::vector<iptype> _int_precision_ushiftright(const std::vector<iptype> *, const size_t);
	std::vector<iptype> _int_precision_ushiftleft(const std::vector<iptype> *, const size_t);
	std::vector<iptype> _int_precision_unegate(const std::vector<iptype>&);
	std::vector<iptype> _int_precision_uand(const std::vector<iptype>&, const std::vector<iptype>&);
	std::vector<iptype> _int_precision_uor(const std::vector<iptype>&, const std::vector<iptype>&);
	std::vector<iptype> _int_precision_uxor(const std::vector<iptype>&, const std::vector<iptype>&);

	int _int_precision_compare2(std::vector<iptype>&, std::vector<iptype>&);
	int _int_precision_compare(const std::vector<iptype> *, const std::vector<iptype> *);
	void _int_precision_strip_leading_zeros(std::vector<iptype> *);
	void _int_precision_strip_trailing_zeros(std::vector<iptype> *);


//////////////////////////////////////////////////////////////////////
//
//
//    Constructors
//
//
//////////////////////////////////////////////////////////////////////


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		int_precision::int_precision
//	@return		nothing
//	@param		"c"	-	the character integer to convert to multi precision number
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initilize with a signed character
//
  /*
   template<class _TY> inline int_precision::int_precision(_TY c)
   {
	   const unsigned int md = sizeof(c) / sizeof(iptype);
	   if (c < (_TY)0)
	   {
		   mSign = -1; c = -c;
	   }
	   else
		   mSign = 1;
	   if (md > 1)
	   {
		   for (int i = 0; i < md; ++i)
		   {
			   mNumber[i] = (iptype)c; c >>= (_TY)Bitsiptype;
		   }
	   }
	   else
		   mNumber.push_back(c);
   }
   */

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Mar/2022
//	@brief 		int_precision::int_precision
//	@return		nothing
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initilize with an int_precision
//
inline int_precision::int_precision()
	{
	mSign = +1;
	mLimit = SIZE_MAX;
	mNumber.assign(1, 0);
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Mar/2022
//	@brief 		int_precision::int_precision
//	@return		nothing
//	@param		"s"	-	the int_precision variable to assign to this
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initilize with an int_precision
//
inline int_precision::int_precision(const int_precision& s, const size_t limit)
	{
	mSign = s.mSign;
	mLimit = limit;
	mNumber = s.mNumber;
	if (mNumber.size() > mLimit)
		mNumber.resize(mLimit);
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Mar/2022
//	@brief 		int_precision::int_precision
//	@return		nothing
//	@param		"v"	-	the vector<iptype> to assign to mNumber
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initilize with a vector<iptype> 
//
inline int_precision::int_precision( const std::vector<iptype>& v, const size_t limit) 
	{
	mSign = 1;
	mLimit = limit;
	mNumber = v;
	if (mNumber.size() > mLimit)
		mNumber.resize(mLimit);
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		int_precision::int_precision
//	@return 	nothing
//	@param		"str"	-	Convert the character string number into a multi precision number
//
//	@todo
//
// Description:
//  Constructor
//  Validate input and convert to internal representation
//
inline int_precision::int_precision(const char *str, const  size_t limit)
	{
	std::string s(str);
	if (strlen(str) == 0)
		{ throw bad_int_syntax(); return; }
	mLimit = limit;
	mNumber = _int_precision_atoip(s, &mSign);
	if (mNumber.size() > mLimit)
		mNumber.resize(mLimit);
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		int_precision::int_precision
//	@return 	nothing
//	@param		"str"	-	Convert the std::string number into a multi precision number
//
//	@todo
//
// Description:
//  Constructor
//  Validate input and convert to internal representation
//
//
inline int_precision::int_precision(const std::string& str, const size_t limit)
	{
	if (str.empty())
		{ throw bad_int_syntax(); return; }
	
	mLimit = limit;
	mNumber= _int_precision_atoip(str, &mSign );
	if (mNumber.size() > mLimit)
		mNumber.resize(mLimit);
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		int_precision::int_precision
//	@return		nothing
//	@param		"c"	-	the character integer to convert to multi precision number
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initilize with a signed character
//
inline int_precision::int_precision( char c, const size_t limit)
	{
	mSign = 1;
	mLimit = limit;
	if (c < 0)
		{
		mSign = -1; c = -c;
		}
	mNumber.push_back(c);
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		int_precision::int_precision
//	@return 	nothing
//	@param		"uc"	-	the character integer to convert to multi precision number
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initilize with an unsigned character
//
inline int_precision::int_precision( unsigned char uc, const size_t limit )
	{
	mSign = 1;
	mLimit = limit;
	mNumber.push_back(uc);
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		int_precision::int_precision
//	@return 	nothing
//	@param		"us"	-	the binary integer to convert to multi precision number
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initialize with an integer
//
inline int_precision::int_precision( short us, const size_t limit )
	{
	mSign = 1;
	mLimit = limit;
	if (us < 0)
		{
		mSign = -1; us = -us;
		}
	mNumber.push_back(us);
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		int_precision::int_precision
//	@return 	nothing
//	@param		"s"	-	the binary integer to convert to multi precision number
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initialize with an integer
//
inline int_precision::int_precision( unsigned short s, const size_t limit )
	{
	mSign = 1;
	mLimit = limit;
	mNumber.push_back(s);
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		int_precision::int_precision
//	@return 	nothing
//	@param		"i"	-	the binary integer to convert to multi precision number
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initialize with an integer
//
inline int_precision::int_precision( int i, const size_t limit)
	{
	mSign = 1;
	mLimit = limit;
	if (i < 0)
		{
		mSign = -1; i = -i;
		}
	mNumber.push_back(i);
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		int_precision::int_precision
//	@return		nothing
//	@param		"ui"	-	the binary unsigned integer to convert to multi precision number
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initialize with an integer
//
inline int_precision::int_precision( unsigned int ui, const size_t limit)
	{
	mSign = 1;
	mLimit = limit;
	mNumber.push_back(ui);
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		int_precision::int_precision
//	@return 	nothing
//	@param		"l"	-	the binary long integer to convert to multi precision number
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initialize with an integer
//
inline int_precision::int_precision( long l, const size_t limit)
	{
	mSign = 1;
	mLimit = limit;
	if (l < 0)
		{
		mSign = -1; l = -l;
		}
	mNumber.push_back(l);
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		int_precision::int_precision
//	@return 	nothing
//	@param		"ul"	-	the binary unsigned long to convert to multi precision number
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initialize with an integer
//
inline int_precision::int_precision( unsigned long ul, const size_t limit)
	{
	mSign = 1;
	mLimit = limit;
	mNumber.push_back(ul);
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		int_precision::int_precision
//	@return 	nothing
//	@param		"ll"	-	the binary long long to convert to multi precision number
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initialize with an integer
//
inline int_precision::int_precision(long long ll, const size_t limit)
	{
	mSign = 1;
	mLimit = limit;
	if (ll < 0)
		{
		mSign = -1; ll = -ll;
		}
	mNumber.push_back(ll);
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		int_precision::int_precision
//	@return 	nothing
//	@param		"ull"	-	the binary unsigned long long to convert to multi precision number
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initialize with an integer
//
inline int_precision::int_precision(unsigned long long ull, const size_t limit)
	{
	mSign = 1;
	mLimit = limit;
	mNumber.push_back(ull);
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		26/Aug/2022
//	@brief 		int_precision::int_precision
//	@return 	nothing
//	@param		"d"	-	the floatto convert to int_precision number
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initialize with a float
//
inline int_precision::int_precision(float d, const size_t limit)
	{// Call the constructor for double
	*this = int_precision((double)d,limit);
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		26/Aug/2022
//	@brief 		int_precision::int_precision
//	@return 	nothing
//	@param		"d"	-	the double to convert to int_precision number
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initialize with a double
//
inline int_precision::int_precision(double d, const size_t limit)
	{
	int expo;
	uintmax_t fpb;

	mNumber.resize(1, 0);
	mSign = +1;
	mLimit = limit;
	if (d<1.0&&d>-1.0)
		return;			// return 0
						// d>=|1|
	if (d < 0)
		{
		mSign = -1; d = -d;
		}
	// d>=1 therefore expo>=0
	fpb = *(uintmax_t *)&d;
	expo = (fpb >> 52) & 0x7ff;	// Extract the exponent
	expo -= 1023;				// unbiased the double exponent
								// Put the imaginary 1 in front of the number
	fpb &= 0xfffffffffffff;		// Mask out exponent  to get the mantissa
	fpb |= 0x10000000000000ull;	// Add the implicit 1  (1ull << 52)
	if (expo <= 52)
		{
		fpb >>= 52 - expo;
		expo = 0;
		}
	else
		expo -= 52;
	mNumber[0] = fpb;
	if (expo > 0)	// Do the remaining left shift
		mNumber = _int_precision_ushiftleft(&mNumber, expo);
	if (mNumber.size() > mLimit)
		mNumber.resize(mLimit);
/*
	mSign = 1;
	*this = _int_precision_dtoip(d);
	mLimit = limit;
	*/
	}

//////////////////////////////////////////////////////////////////////
//
//
//    Implicit conversions to base types: int, short, long, char, float & double
//
//
//////////////////////////////////////////////////////////////////////

#ifdef TEMPLIFY
template<class _TY> inline int_precision::operator _TY() const
	{
	return (_TY)(mNumber[0] * mSign);
	}
#else
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		int_precision::operator long long
//	@return 	long long -
//
//	@todo
//
// Description:
//  This is the main operator from int_precision to regular long, int, short & char
//  Any explicit or implicit copnversion first convert to standard c long type and then to any other
//  inbuild type long long, long, int, short, char. As a type long long >= long >= int >= short >= char
//  As with regular C type conversion the conversion truncate to desired type and a possible
//  loss in precision is possible
//
inline int_precision::operator long long() const
	{// Conversion to long long
	return (long long)(mNumber[0] * mSign);
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/20221
//	@brief 		int_precision::operator long
//	@return 	long 
//
//	@todo
//
// Description:
//  This is the main operator from int_precision to regular long, int, short & char
//  Any explicit or implicit copnversion first convert to standard c long type and then to any other
//  inbuild type int, short, char. As a type long >= int >= short >= char
//  As with regular C type conversion the conversion truncate to desired type and a possible
//  loss in precision is possible
//
inline int_precision::operator long() const
	{// Conversion to long
	return (long)(mNumber[0] * mSign);
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		int_precision::operator int
//	@return 	int	-
//
//	@todo  Add to do things
//
// Description:
//  Any explicit or implicit copnversion first convert to standard c long type and then to int
//  As with regular C type conversion the conversion truncate to desired type and a possible
//  loss in precision is possible
//
inline int_precision::operator int() const
	{// Conversion to int
	return (int)(mNumber[0] * mSign);
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		int_precision::operator short
//	@return 	short	-
//
//	@todo  Add to do things
//
// Description:
//  Any explicit or implicit copnversion first convert to standard c long type and then to short
//  As with regular C type conversion the conversion truncate to desired type and a possible
//  loss in precision is possible
//
inline int_precision::operator short() const
	{// Conversion to short
	return (short)(mNumber[0] * mSign);
    }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		int_precision::operator char
//	@return 	char -
//
//	@todo  Add to do things
//
// Description:
//  Any explicit or implicit copnversion first convert to standard c long type and then to char
//  As with regular C type conversion the conversion truncate to desired type and a possible
//  loss in precision is possible
//
inline int_precision::operator char() const
	{// Conversion to char
	return (char)(mNumber[0] * mSign); 
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		int_precision::operator unsigned long long
//	@return 	unsigned long long -
//
//	@todo  Add to do things
//
// Description:
//  Any explicit or implicit copnversion first convert to standard c long type and then to int
//  As with regular C type conversion the conversion truncate to desired type and a possible
//  loss in precision is possible
//
inline int_precision::operator unsigned long long() const
	{// Conversion to unsigned long long
	return (unsigned long long)(mNumber[0]);
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		int_precision::operator unsigned long
//	@return 	unsigned long -
//
//	@todo  Add to do things
//
// Description:
//  Any explicit or implicit copnversion first convert to standard c long type and then to int
//  As with regular C type conversion the conversion truncate to desired type and a possible
//  loss in precision is possible
//
inline int_precision::operator unsigned long() const
	{// Conversion to unsigned long
	return (unsigned long)(mNumber[0] );
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8-Aug-2021
//	@brief 		int_precision::operator unsigned int
//	@return 	unsigned int -
//
//	@todo  Add to do things
//
// Description:
//  Any explicit or implicit copnversion first convert to standard c long type and then to int
//  As with regular C type conversion the conversion truncate to desired type and a possible
//  loss in precision is possible
//
inline int_precision::operator unsigned int() const
	{// Conversion to unsigned int
	return (unsigned int)(mNumber[0] );
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		int_precision::operator unsigned short
//	@return 	unsigned short -
//
//	@todo  Add to do things
//
// Description:
//  Any explicit or implicit copnversion first convert to standard c long type and then to short
//  As with regular C type conversion the conversion truncate to desired type and a possible
//  loss in precision is possible
//
inline int_precision::operator unsigned short() const
	{// Conversion to unsigned short
	return (unsigned short)(mNumber[0] );
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		int_precision::operator unsigned char
//	@return 	unsigned char -
//
//	@todo  Add to do things
//
// Description:
//  Any explicit or implicit copnversion first convert to standard c long type and then to char
//  As with regular C type conversion the conversion truncate to desired type and a possible
//  loss in precision is possible
//
inline int_precision::operator unsigned char() const
	{// Conversion to char
	return (unsigned char)(mNumber[0] );
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		9/Aug/2021
//	@brief 		int_preceision::operator double
//	@return 	return double -
//
//	@todo
//
// Description:
//  Conversion from int_precision to double
//
inline int_precision::operator double() const
	{// Conversion to double
	return _int_precision_iptod(this);
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		9/Aug/2021
//	@brief 		int_precision::operator float
//	@return 	return float -
//
//	@todo  Add to do things
//
// Description:
//  Conversion from int_precision to float
//  Using the double conversion frist and then trunk to float using standard c conversion
//
inline int_precision::operator float() const
	{// Conversion to float
	return (float)((double)*this);
	}
#endif

//////////////////////////////////////////////////////////////////////
//
//    Class Methods:
//			pointer
//			number
//			sign
//			change_sign
//			size
//			abs
//			setbit
//			resetbit
//			flipbit
//			testbit
//			ctz
//			clz
//			even
//			odd
//			iszero
//			toString
//
//////////////////////////////////////////////////////////////////////


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Nov/2021
//	@brief 		int_precision::pointer
//	@return 	return a pointer to mNumber
//
//	@todo  Add to do things
//
// Description:
//  Return a pointer to mNumber
//
inline std::vector<iptype> *int_precision::pointer() { return &mNumber; }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Nov/2021
//	@brief 		int_precision::number
//	@return 	return a copy of mNumber
//
//	@todo  Add to do things
//
// Description:
//  Return a copy of mNumber
//
inline std::vector<iptype> int_precision::number() const { return mNumber; }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Nov/2021
//	@brief 		int_precision::number
//	@param		"mNumber"	-	New mNumber vector<iptype>
//	@return 	Set and return a copy of mNumber
//
//	@todo  Add to do things
//
// Description:
//  Set and Return a copy of mNumber
//
inline std::vector<iptype> int_precision::number(std::vector<iptype> &mb) { return mNumber = mb; }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Nov/2021
//	@brief 		int_precision::index
//	@param		"inx"	-	Index into mNumber vector<iptype>
//	@return 	return the index of mNumber[inx]
//
//	@todo  Add to do things
//
// Description:
//  Return the indx of mNumber[inx]
//
inline iptype int_precision::index(const size_t inx) const { return mNumber[inx]; }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Nov/2021
//	@brief 		int_precision::sign
//	@return 	return a copy of the sign
//
//	@todo  Add to do things
//
// Description:
//  Return a copy of mSign
//
inline int int_precision::sign() const { return mSign; }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Nov/2021
//	@brief 		int_precision::number
//	@param		"mSign"	-	New mSign
//	@return 	Set and return a copy of the new sign
//
//	@todo  Add to do things
//
// Description:
//  Set and Return a copy of mSign
//
inline int int_precision::sign(int s) { return (mSign = s); }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Nov/2021
//	@brief 		int_precision::change_sign
//	@return 	Change and return a copy of the sign
//
//	@todo  Add to do things
//
// Description:
//  Change sign and Return a copy of mSign
//
inline int int_precision::change_sign() { mSign *= -1;  return mSign; }		// Toggle and return sign 

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Nov/2021
//	@brief 		int_precision::size
//	@return 	return the size of the mNumber vector<iptype>
//
//	@todo  Add to do things
//
// Description:
// Return the size of the mNumber vector<iptype>
//
inline size_t int_precision::size() const { return mNumber.size(); }		// Return the actual number of iptype digits

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Nov/2021
//	@brief 		int_precision::precision
//	@return 	return the maximum size the mNumber vector<iptype> cn hold
//
//	@todo  Add to do things
//
// Description:
// Return the maximum precision of the mNumber vector<iptype>
//
inline size_t int_precision::precision() const { return mLimit; }		// Return the maximum number of iptype digits

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Mar/2022
//	@brief 		int_precision::precision
//	@return 	return the size of the mNumber vector<iptype>
//  @param		"p"		-- The new fixed size integer preicsion or arbitrary precision
//
//	@todo  Add to do things
//
// Description:
// Set and Return the new maximum precision the mNumber vector<iptype> can hold
//
inline size_t int_precision::precision(const size_t p )  
	{ 
	mLimit = p==0 ? 1 : p;  //  Can only be set to a size >= 1
	if (mLimit < mNumber.size())
		mNumber.resize(mLimit);
	return mLimit; // Return the new number of maximum iptype digits mNumber can hold 
	}		

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Nov/2021
//	@brief 		int_precision::abs
//	@return 	return the absolute value of the int_precision object
//
//	@todo  Add to do things
//
// Description:
//  Return te absolue value of the int_precision object
//
inline int_precision& int_precision::abs() { mSign = 1; return *this; }		// Change sign to + and return number

// Start of Bit Methods
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Nov/2021
//	@brief 		int_precision::setbit
//	@param		"i"	-	Bit position
//	@return 	void
//
//	@todo  Add to do things
//
// Description:
//  Set the bit in mNumber at bit position i
//
inline void int_precision::setbit(size_t i) {
	size_t n = i / Bitsiptype;
	if (n >= mNumber.size())
		mNumber.insert(mNumber.end(), n + 1 - mNumber.size(), 0);
	mNumber[n] |= (iptype)(1) << (i % Bitsiptype);
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Nov/2021
//	@brief 		int_precision::resetbit
//	@param		"i"	-	Bit position
//	@return 	void
//
//	@todo  Add to do things
//
// Description:
// Resetset the bit in mNumber at bit position i
//
inline void int_precision::resetbit(size_t i) {
	size_t n = i / Bitsiptype;
	if (n >= mNumber.size())
		mNumber.insert(mNumber.end(), n + 1 - mNumber.size(), 0);
	mNumber[n] &= ~((iptype)(1) << (i % Bitsiptype));
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Nov/2021
//	@brief 		int_precision::flipbit
//	@param		"i"	-	Bit position
//	@return 	void
//
//	@todo  Add to do things
//
// Description:
//  Flip the bit in mNumber at bit position i
//
inline void int_precision::flipbit(size_t i) {
	size_t n = i / Bitsiptype;
	if (n >= mNumber.size())
		mNumber.insert(mNumber.end(), n + 1 - mNumber.size(), 0);
	mNumber[n] ^= ~((iptype)(1) << (i % Bitsiptype));
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Nov/2021
//	@brief 		int_precision::testbit
//	@param		"i"	-	Bit position
//	@return 	boolean true or false
//
//	@todo  Add to do things
//
// Description:
//  Test the bit in mNumber at bit position i and return true if bit is set otherwise false
//
inline bool int_precision::testbit(size_t i) {
	size_t n = i / Bitsiptype;
	if (n >= mNumber.size())
		mNumber.insert(mNumber.end(), n + 1 - mNumber.size(), 0);
	return mNumber[n] & (iptype)(1) << (i % Bitsiptype) ? true : false;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Nov/2021
//	@brief 		int_precision::ctz
//	@param		"i"	-	Bit position
//	@return 	number of trailing zero bit in mNumber 
//
//	@todo  Add to do things
//
// Description:
// Return the number of trailing zero bits in mNumber vector<iptype>
//
inline size_t int_precision::ctz() { return _int_precision_ctz(mNumber); }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Nov/2021
//	@brief 		int_precision::ctz
//	@param		"i"	-	Bit position
//	@return 	number of leadingzero  bit in mNumber 
//
//	@todo  Add to do things
//
// Description:
// Return the number of leading zero bits in mNumber vector<iptype>
//
inline size_t int_precision::clz() { return _int_precision_clz(mNumber); }

// End of bit methods

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Nov/2021
//	@brief 		int_precision::even
//	@return 	return true or false if mNumber is even or odd 
//
//	@todo  Add to do things
//
// Description:
// Return true if mNumber number is even otherwise false
//
inline bool int_precision::even() const { return (mNumber[0] & 0x1) ? false : true; }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Nov/2021
//	@brief 		int_precision::odfd
//	@return 	return true or false if mNumber is even or odd 
//
//	@todo  Add to do things
//
// Description:
// Return true if mNumber number is odd otherwise false
//
inline bool int_precision::odd() const { return (mNumber[0] & 0x1) ? true : false; }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Nov/2021
//	@brief 		int_precision::iszero
//	@return 	return true if mNumber number is zero
//
//	@todo  Add to do things
//
// Description:
// Return true if mNumber number is zero
//
inline bool int_precision::iszero()	const { return mNumber.size() == 1 && mNumber[0] == 0 ? true : false; }

// Conversion methods. Safer and less ambiguous than overloading implicit/explicit conversion operators
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Nov/2021
//	@brief 		int_precision::toString
//	@return 	return the decimal string number of the int_precision object
//
//	@todo  Add to do things
//
// Description:
// Return the decimal string number of the int_precision object
//
inline std::string int_precision::toString(const int base=BASE_10) { return _int_precision_itoa(this, base); }


//////////////////////////////////////////////////////////////////////
//
//
//    Essentialsoperators =, +=, -=, *=, /=, %=, <<=, >>=, &=, |=, ^=
//
//
//////////////////////////////////////////////////////////////////////

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/19/2005
//	@brief 		operator=
//	@return 	static int_precision	-	return a=b
//	@param		"a"	-	Assignment operand
//
//	@todo
//
// Description:
//  Assign operator
//
inline int_precision& int_precision::operator=( const int_precision& a )
	{
	mSign = a.mSign;
	mNumber = a.mNumber;
	if (mNumber.size() > mLimit)
		mNumber.resize(mLimit);
	return *this;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Aug/2021
//	@brief 		operator+=
//	@return		static int_precision	-	return a +=b
//	@param		"a"	-	Adding operand
//
//	@todo
//
// Description:
//  += operator. 
//
inline int_precision& int_precision::operator+=( const int_precision& a )
	{
	int wrap, cmp;

	if( a.mSign == mSign )
		mNumber = _int_precision_uadd( (std::vector<iptype> *)&a.mNumber, &mNumber );  // Add and no change of sign
	else
		{
		cmp = _int_precision_compare( (std::vector<iptype> *)&a.mNumber, &mNumber );
		//cmp = _int_precision_compare(const_cast<std::vector<iptype>& > (a.mNumber), mNumber);
		if (cmp > 0) // Since we subctract less the wrap indicater need not to be checked
			{
			mSign = a.mSign;
			mNumber = _int_precision_usub(&wrap, (std::vector<iptype> *)&a.mNumber, &mNumber);  // Subtract and change to sign1
			}
		else
			if( cmp < 0 )
				mNumber = _int_precision_usub( &wrap, &mNumber, (std::vector<iptype> *)&a.mNumber ); // Subtract and no change in sign
			else
				{// result is 0
				mSign = +1;  // Change to + sign, since -0 is not allowed for the internal representation
				mNumber = std::vector<iptype>(1,0);
				}
		}

	if (mNumber.size() > mLimit)
		mNumber.resize(mLimit);
	
	return *this;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/19/2005
//	@brief 		operator-=
//	@return		static int_precision	-	return a -=b
//	@param		"a"	-	Subtracting operand
//
//	@todo
//
// Description:
//  -= operator
//  The essential -= operator
//  n = n - a is the same as n = n + (-a);
//
inline int_precision& int_precision::operator-=( const int_precision& a )
	{
	int_precision b;

	b = a;
	b.change_sign();
	*this += b;

	return *this;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Aug/2021
//	@brief 		operator*=
//	@return 	static int_precision	-	return a *=b
//	@param		"a"	-	Multiplying operand
//
//	@todo
//
// Description:
//  *= operator
//
inline int_precision& int_precision::operator*=( const int_precision& a )
	{
	mSign *= a.mSign;  // Resulting sign
	mNumber = _int_precision_umul(&mNumber, (std::vector<iptype> *)&a.mNumber);
	if (mSign == -1 && mNumber.size() == 1 && mNumber[0] == (iptype)0)  // Avoid -0 as result +0 is right
		mSign = +1;

	if (mNumber.size() > mLimit)
		mNumber.resize(mLimit);
	
	return *this;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Aug/2021
//	@brief 		operator/=
//	@return 	static int_precision	-	return a /=b
//	@param		"a"	-	Dividing operand
//
//	@todo
//
// Description:
//  /= operator
//
inline int_precision& int_precision::operator/=( const int_precision& a )
	{
	iptype wrap;
#ifdef _INT_PRECISSION_FAST_DIV_REM
	// Not yet tested
	if (this->size()>a.size() + 8 && a.size() != 1 )  // Check that lhs is 8 digit larger and that rhs is not a single digit before do the fastremdiv operation
		{
		extern int_precision _int_precision_fastdiv( const int_precision&, const int_precision& );
		int_precision b=*this;
		*this = _int_precision_fastdiv( b, a );
		return *this;
		}
#endif

	mSign *= a.mSign;  // Resulting sign after division
	if (a.mNumber.size() == 1 && (a.mNumber[0]>>32)==0 ) // Make short div if denominator <= 32 bit integer.
		mNumber = _int_precision_udiv_short( &wrap, &mNumber, a.mNumber[0]);
	else
		{// Check for division of of number that can safely be done using 64bit binary division
		mNumber = _int_precision_udiv(&mNumber, (std::vector<iptype> *)&a.mNumber);
 		}
	
	if (mSign == -1 && mNumber.size() == 1 && mNumber[0] == (iptype)0)  // Avoid -0 as result +0 is right
		mSign = +1;

	if (mNumber.size() > mLimit)
		mNumber.resize(mLimit);
   
	return *this;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Aug/2021
//	@brief 		operator%=
//	@return 	static int_precision	-	return a %=b
//	@param		"a"	-	Modulus operand
//
//	@todo
//
// Description:
//  %= operator
//
inline int_precision& int_precision::operator%=( const int_precision& a )
	{
#ifdef _INT_PRECISSION_FAST_DIV_REM
	if(this->size()>a.size()+8 && a.size() != 2 )  // Check that lhs is 8 digit larger and that rhs is not a single digit before do the fastremdiv operation
		{
		extern int_precision _int_precision_fastrem( const int_precision&, const int_precision& );
		int_precision b=*this;
		*this =_int_precision_fastrem( b, a );
		return *this;
		}
#endif

	if (a.mNumber.size() == 1 && (a.mNumber.front() >> 32) == 0) // Make short rem 
		mNumber = _int_precision_urem_short( &mNumber, a.mNumber[0]);  // Short rem and sign stay the same
	else
		// Check for remainder of of number that can safely be done using 64bit binary remainder
		mNumber = _int_precision_urem(&mNumber, (std::vector<iptype> *)&a.mNumber);	// regular rem. sign stay the same
   
	if (mSign == -1 && mNumber.size() == 1 && mNumber[0] == (iptype)0)  // Avoid -0 as result +0 is right
	   mSign = +1;

	if (mNumber.size() > mLimit)
		mNumber.resize(mLimit);

	return *this;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Aug/2021
//	@brief 		operator<<=
//	@return		static int_precision	-	return shifting a<<= b
//	@param		"a"	-	Shifting number
//
//	@todo
//
// Description:
//  <<= operator
//
inline int_precision& int_precision::operator<<=( const int_precision& a )
	{
	unsigned int shift;
	shift = (unsigned int)(a.mNumber[0]);
	mNumber = _int_precision_ushiftleft(&mNumber, shift);
	if (mNumber.size() > mLimit)
		mNumber.resize(mLimit);

   return *this;
   }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Aug/2021
//	@brief 		operator>>=
//	@return		int_precision	-	return shifting a>>= b
//	@param		"a"	-	Shifting operand
//
//	@todo
//
// Description:
//  >>= operator
//
inline int_precision& int_precision::operator>>=( const int_precision& a )
	{
	unsigned int shift;
	shift = (unsigned int)(a.mNumber[0]);
	mNumber = _int_precision_ushiftright(&mNumber, shift);
	if (mNumber.size() > mLimit)
		mNumber.resize(mLimit);

	return *this;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Aug/2021
//	@brief 		operator&=
//	@return 	int_precision	-	return a &=b
//	@param		"a"	-	Anding operand
//
//	@todo
//
// Description:
//  &= operator
//
inline int_precision& int_precision::operator&=( const int_precision& a )
   {
   mNumber = _int_precision_uand( mNumber, a.mNumber);
   if (mNumber.size() > mLimit)
	   mNumber.resize(mLimit);

   return *this;
   }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Aug/2021
//	@brief 		operator|=
//	@return 	int_precision	-	return a |=b
//	@param		"a"	-	Oring operand
//
//	@todo
//
// Description:
//  |= operator
//
inline int_precision& int_precision::operator|=( const int_precision& a)
	{
	mNumber = _int_precision_uor( mNumber, a.mNumber);
	if (mNumber.size() > mLimit)
		mNumber.resize(mLimit);

	return *this;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Aug/2021
//	@brief 		operator^=
//	@return 	int_precision	-	return a ^=b
//	@param		"a"	-	Xoring operand
//
//	@todo
//
// Description:
//  ^= operator
//
inline int_precision& int_precision::operator^=(const int_precision& a)
	{
	mNumber = _int_precision_uxor( mNumber, a.mNumber);
	if (mNumber.size() > mLimit)
		mNumber.resize(mLimit);

	return *this;
	}

//////////////////////////////////////////////////////////////////////
//
//
//    Arithmetic
//
//
//////////////////////////////////////////////////////////////////////


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		2/19/2006
//	@brief 		operator+
//	@return 	int_precision	-	return addition of lhs + rhs
//	@param		"lhs"	-	First operand
//	@param		"rhs"	-	Second operand
//
//	@todo  Add to do things
//
// Description:
//  Add operator for int_precision + <any other type>
//  no const on the lhs parameter to prevent ambigous overload
//
template <class _Ty> inline int_precision operator+( int_precision& lhs, const _Ty& rhs )
	{
	return int_precision(lhs) += int_precision(rhs);
	}


// @author Henrik Vestermark (hve@hvks.com)
// @date		2/19/2006
//	@brief 		operator+
//	@return 	int_precision	-	return addition of lhs + rhs
//	@param		"lhs"	-	First operand
//	@param		"rhs"	-	Second operand
//
// Description:
//  Add operator for <any other type> + int_precision
//
template <class _Ty> inline int_precision operator+( const _Ty& lhs, const int_precision& rhs )
	{
	return int_precision(lhs) += rhs;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/19/2005
//	@brief 		operator unary +
//	@return 	int_precision	-	a
//	@param		"a"	-	operand
//
//	@todo
//
// Description:
//  Unary + operator
//  Do nothing
//
inline int_precision operator+( const int_precision& a )
	{
	return a;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/19/2005
//	@brief 		operator++ Prefix
//	@return 	int_precision	-	return the incremented a
//	@param		"a"	-	operand
//
//	@todo
//
// Description:
//  Increment operator
//
inline int_precision operator++( int_precision& a )
	{
	a += int_precision( 1 );
	return a;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/14/2005
//	@brief 		operator++ Postfix
//	@return 	int_precision	-	return the a before incrementation
//	@param		"a"	-	operand
//
//	@todo
//
// Description:
//  Postfix Increment operator
//
inline int_precision operator++( int_precision& a, int )
	{
	int_precision postfix_a(a);

	a += int_precision( 1 );
	return postfix_a;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		2/19/2006
//	@brief 		operator-
//	@return 	int_precision	-	return addition of lhs + rhs
//	@param		"lhs"	-	First operand
//	@param		"rhs"	-	Second operand
//
//	@todo  Add to do things
//
// Description:
//  Add operator for int_precision - <any other type>
//  no const on the lhs parameter to prevent ambigous overload
//
template <class _Ty> inline int_precision operator-( int_precision& lhs, const _Ty& rhs )
	{
	return int_precision(lhs) -= int_precision(rhs);
	}


// @author Henrik Vestermark (hve@hvks.com)
// @date		2/19/2006
//	@brief 		operator-
//	@return 	int_precision	-	return addition of lhs - rhs
//	@param		"lhs"	-	First operand
//	@param		"rhs"	-	Second operand
//
// Description:
//  Add operator for <any other type> - int_precision
//
template <class _Ty> inline int_precision operator-( const _Ty& lhs, const int_precision& rhs )
	{
	return int_precision(lhs) -= rhs;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/19/2005
//	@brief 		operator unary -
//	@return 	int_precision	-	-a
//	@param		"a"	-	operand for sign change
//
//	@todo
//
// Description:
//  Unary - operator
//  Change sign
//
inline int_precision operator-( const int_precision& a )
	{
	int_precision b(a);
	b.change_sign();
	return b;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/19/2005
//	@brief 		operator-- prefix
//	@return 	int_precision	-	return the decremented a
//	@param		"a"	-	operand
//
//	@todo
//
// Description:
//  Decrement operator
//
inline int_precision operator--( int_precision& a )
	{
	a -= int_precision( 1 );
	return a;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/14/2005
//	@brief 		operator-- postfix
//	@return 	int_precision	-	return the a before decrementation
//	@param		"a"	-	operand
//
//	@todo
//
// Description:
//  Postfix Decrement operator
//
int_precision operator--( int_precision& a, int )
	{
	int_precision postfix_a(a);
	a -= int_precision( 1 );
	return postfix_a;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		2/19/2006
//	@brief 		operator*
//	@return 	int_precision	-	return addition of lhs + rhs
//	@param		"lhs"	-	First operand
//	@param		"rhs"	-	Second operand
//
//	@todo  Add to do things
//
// Description:
//  Add operator for int_precision * <any other type>
//
template <class _Ty> inline int_precision operator*( int_precision& lhs, const _Ty& rhs )
	{
	return int_precision(lhs) *= int_precision(rhs);
	}


// @author Henrik Vestermark (hve@hvks.com)
// @date		2/19/2006
//	@brief 		operator*
//	@return 	int_precision	-	return addition of lhs - rhs
//	@param		"lhs"	-	First operand
//	@param		"rhs"	-	Second operand
//
// Description:
//  Add operator for <any other type> * int_precision
//  no const on the lhs parameter to prevent ambigous overload
//
template <class _Ty> inline int_precision operator*( const _Ty& lhs, const int_precision& rhs )
	{
	return int_precision(lhs) *= rhs;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		2/19/2006
//	@brief 		operator/
//	@return 	int_precision	-	return addition of lhs + rhs
//	@param		"lhs"	-	First operand
//	@param		"rhs"	-	Second operand
//
//	@todo  Add to do things
//
// Description:
//  Add operator for int_precision / <any other type>
//  no const on the lhs parameter to prevent ambigous overload
//
template <class _Ty> inline int_precision operator/( int_precision& lhs, const _Ty& rhs )
	{
	return int_precision(lhs) /= int_precision(rhs);
	}


// @author Henrik Vestermark (hve@hvks.com)
// @date		2/19/2006
//	@brief 		operator*
//	@return 	int_precision	-	return addition of lhs - rhs
//	@param		"lhs"	-	First operand
//	@param		"rhs"	-	Second operand
//
// Description:
//  Add operator for <any other type> / int_precision
//
template <class _Ty> inline int_precision operator/( const _Ty& lhs, const int_precision& rhs )
	{
	return int_precision(lhs) /= rhs;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		2/19/2006
//	@brief 		operator%
//	@return 	int_precision	-	return addition of lhs + rhs
//	@param		"lhs"	-	First operand
//	@param		"rhs"	-	Second operand
//
//	@todo  Add to do things
//
// Description:
//  Add operator for int_precision % <any other type>
//
template <class _Ty> inline int_precision operator%( int_precision& lhs, const _Ty& rhs )
	{
	return int_precision(lhs) %= int_precision(rhs);
	}


// @author Henrik Vestermark (hve@hvks.com)
// @date		2/19/2006
//	@brief 		operator%
//	@return 	int_precision	-	return addition of lhs - rhs
//	@param		"lhs"	-	First operand
//	@param		"rhs"	-	Second operand
//
// Description:
//  Add operator for <any other type> % int_precision
//
template <class _Ty> inline int_precision operator%( const _Ty& lhs, const int_precision& rhs )
	{
	return int_precision(lhs) %= rhs;
	}

// @author Henrik Vestermark (hve@hvks.com)
// @date		5/sep/2021
//	@brief 		operator<<
//	@return 	int_precision	-	return addition of lhs - rhs
//	@param		"lhs"	-	First operand
//	@param		"rhs"	-	Second operand
//
// Description:
//  Add operator for int_precision << int_precision
//
inline int_precision operator<<(const int_precision& lhs, const int_precision& rhs)
	{
	return int_precision(lhs) <<= rhs;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		2/20/2006
//	@brief 		operator<<
//	@return 	int_precision	-	return addition of lhs + rhs
//	@param		"lhs"	-	First operand
//	@param		"rhs"	-	Second operand
//
//	@todo  Add to do things
//
// Description:
//  Add operator for int_precision << <any other type>
//
template <class _Ty> inline int_precision operator<<( int_precision& lhs, const _Ty& rhs )
	{
	return int_precision(lhs) <<= int_precision(rhs);
	}


// @author Henrik Vestermark (hve@hvks.com)
// @date		2/20/2006
//	@brief 		operator<<
//	@return 	int_precision	-	return addition of lhs - rhs
//	@param		"lhs"	-	First operand
//	@param		"rhs"	-	Second operand
//
// Description:
//  Add operator for <any other type> << int_precision
//
//template <class _Ty> inline int_precision operator<<(  const _Ty& lhs, const int_precision& rhs )
//	{
//	return int_precision(lhs) <<= rhs;
//	}




//	@author Henrik Vestermark (hve@hvks.com)
//	@date		2/20/2006
//	@brief 		operator>>
//	@return 	int_precision	-	return addition of lhs + rhs
//	@param		"lhs"	-	First operand
//	@param		"rhs"	-	Second operand
//
//	@todo  Add to do things
//
// Description:
//  Add operator for int_precision >> <any other type>
//
template <class _Ty> inline int_precision operator>>( int_precision& lhs, const _Ty& rhs )
	{
	return int_precision(lhs) >>= int_precision(rhs);
	}


// @author Henrik Vestermark (hve@hvks.com)
// @date		2/20/2006
//	@brief 		operator>>
//	@return 	int_precision	-	return addition of lhs - rhs
//	@param		"lhs"	-	First operand
//	@param		"rhs"	-	Second operand
//
// Description:
//  Add operator for <any other type> >> int_precision
//
template <class _Ty> inline int_precision operator>>( const _Ty& lhs, const int_precision& rhs )
	{
	return int_precision(lhs) >>= rhs;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Aug/2021
//	@brief 		operator&
//	@return 	int_precision	-	return lhs & rhs
//	@param		"lhs"	-	First operand
//	@param		"rhs"	-	Second operand
//
//	@todo  Add to do things
//
// Description:
//  And operator for int_precision & <any other type>
//  no const on the lhs parameter to prevent ambigous overload
//
template <class _Ty> inline int_precision operator&( int_precision& lhs, const _Ty& rhs )
   {
   return int_precision(lhs) &= int_precision(rhs);
   }

// @author Henrik Vestermark (hve@hvks.com)
// @date		11/Aug/2021
//	@brief 		operator&
//	@return 	int_precision	-	return lhs & rhs
//	@param		"lhs"	-	First operand
//	@param		"rhs"	-	Second operand
//
// Description:
//  | operator for <any other type> & int_precision
//
template <class _Ty> inline int_precision operator&( const _Ty& lhs, const int_precision& rhs )
   {
   return int_precision(lhs) &= rhs;
   }

// @author Henrik Vestermark (hve@hvks.com)
// @date		2/Sep/2021
//	@brief 		operator|
//	@return 	int_precision	-	return lhs | rhs
//	@param		"lhs"	-	First operand
//	@param		"rhs"	-	Second operand
//
// Description:
//  | operator for <any other type> | int_precision
//
inline int_precision operator|(const int_precision& lhs, const int_precision& rhs)
	{
	return int_precision(lhs) |= rhs;
	}

// @author Henrik Vestermark (hve@hvks.com)
// @date		2/Sep/2021
//	@brief 		operator|
//	@return 	int_precision	-	return lhs | rhs
//	@param		"lhs"	-	First operand
//	@param		"rhs"	-	Second operand
//
// Description:
//  | operator for <any other type> | int_precision
//
//template <class _Ty> inline int_precision operator|(const _Ty& lhs, const int_precision& rhs)
//	{
//	return int_precision(lhs) |= rhs;
//	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		2/Sep/2021
//	@brief 		operator&
//	@return 	int_precision	-	return lhs | rhs
//	@param		"lhs"	-	First operand
//	@param		"rhs"	-	Second operand
//
//	@todo  Add to do things
//
// Description:
//  And operator for int_precision | <any other type>
//  no const on the lhs parameter to prevent ambigous overload
//
template <class _Ty> inline int_precision operator|(int_precision& lhs, const _Ty& rhs)
	{
	return int_precision(lhs) |= (int_precision)rhs;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/Sep/2021
//	@brief 		operator&
//	@return 	int_precision	-	return lhs ^ rhs
//	@param		"lhs"	-	First operand
//	@param		"rhs"	-	Second operand
//
//	@todo  Add to do things
//
// Description:
//  And operator for int_precision ^ <any other type>
//  no const on the lhs parameter to prevent ambigous overload
//
template <class _Ty> inline int_precision operator^(int_precision& lhs, const _Ty& rhs)
	{	
	return int_precision(lhs) ^= int_precision(rhs);
	}

// @author Henrik Vestermark (hve@hvks.com)
// @date		11/Aug/2021
//	@brief 		operator^
//	@return 	int_precision	-	return lhs ^ rhs
//	@param		"lhs"	-	First operand
//	@param		"rhs"	-	Second operand
//
// Description:
//  ^  operator for <any other type> ^ int_precision
//
template <class _Ty> inline int_precision operator^( const _Ty& lhs, const int_precision& rhs )
	{
	return int_precision(lhs) ^= rhs;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/Sep/2021
//	@brief 		operator unary ~(negate)
//	@return 	int_precision	-	-a
//	@param		"a"	-	operand for ~ (negate) operator
//
//	@todo
//
// Description:
//  Unary ~ operator
//  Negate Integer
//
inline int_precision operator~(const int_precision& a)
	{
	int_precision lhs;
	std::vector<iptype> des = a.number();
	std::vector<iptype>::iterator d_pos;

	for (d_pos = des.begin(); d_pos != des.end(); ++d_pos)
		{ // negating element of the number
		*d_pos = ~*d_pos;
		}
	lhs.number(des);
	return lhs; 
	}

//////////////////////////////////////////////////////////////////////
//
//
//    Comparison
//
//
//////////////////////////////////////////////////////////////////////

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Aug/2021
//	@brief 		operator==
//	@return 	bool	- 	the boolean value of the operator
//	@param		"a"	-	First operand number to compare
//	@param		"b"	-	Second operand number to
//
//	@todo
//
// Description:
//  Boolean equal of two precision numbers. Early out algorithm
//  1) if sign different then result is false. We actual don't do this test because of -0==+0 we should not ocuured but just in case
//  2) if length is different then the result is false
//  3) use core compare to determine boolean value
//
template <class _Ty> inline bool operator==( int_precision& a, const _Ty& b )
	{
	int_precision c(b);
	if( a.sign()==c.sign() &&  _int_precision_compare(const_cast<int_precision&>(a).pointer(), const_cast<int_precision&>(c).pointer())==0)// Same return true
		return true;
	return false;
	//_int_precision_compare(const_cast<int_precision&>(a).pointer(), const_cast<int_precision&>(c).pointer())
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Aug/2021
//	@brief 		operator==
//	@return 	bool	- 	the boolean value of the operator
//	@param		"a"	-	First operand number to compare
//	@param		"b"	-	Second operand number to
//
//	@todo
//
// Description:
//  Boolean equal of two precision numbers. Early out algorithm
//  1) if sign different then result is false. We actual don't do this test because of -0==+0 we should not ocuured but just in case
//  2) if length is different then the result is false
//  3) use core compare to determine boolean value
//
template <class _Ty> inline bool operator==( const _Ty& a, const int_precision& b )
	{
	int_precision c(a);
	if( c.sign()==b.sign() && _int_precision_compare(const_cast<int_precision&>(c).pointer(), const_cast<int_precision&>(b).pointer()) == 0 )    return true;
	return false;
	//_int_precision_compare( const_cast<int_precision&>(c).pointer(), const_cast<int_precision&>(b).pointer() )
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Aug/2021
//	@brief 		operator<
//	@return 	bool	- 	the boolean value of the operator
//	@param		"a"	-	First operand number to compare
//	@param		"b"	-	Second operand number to
//
//	@todo
//
// Description:
//  Boolean less of two precision numbers. Early out algorithm for higher performance
//  1) If sign different determine boolean result based on sign
//  2) Otherwise determine boolean result based length of number amnd the sign
//  3) Same sign and same length. Do a core comparison and return the result
//
template <class _Ty> inline bool operator<( int_precision& a, const _Ty& c )
	{
	int sign1, sign2, cmp;
	int_precision b(c);

	sign1 = a.sign();
	sign2 = b.sign();

	// Different signs
	if( sign1 > sign2 )
		return false;
	if( sign1 < sign2 )
		return true;

	// Same sign
	if( sign1 == 1 && a.size() < b.size() ) // Different therefore true
		return true;
	if( sign1 == 1 && a.size() > b.size() ) // Different therefore false
		return false;
	if( sign1 == -1 && a.size() > b.size() )
		return true;
	if( sign1 == -1 && a.size() < b.size() )
		return false;

	// Same sign and same length
	cmp = _int_precision_compare(const_cast<int_precision&>(a).pointer(), const_cast<int_precision&>(b).pointer());
	if( cmp < 0 && sign1 == 1 )
		return true;
	else
		if( cmp > 0 && sign1 == -1 )
			return true;

	return false;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Aug/2021
//	@brief 		operator<
//	@return 	bool	- 	the boolean value of the operator
//	@param		"a"	-	First operand number to compare
//	@param		"b"	-	Second operand number to
//
//	@todo
//
// Description:
//  Boolean less of two precision numbers. Early out algorithm for higher performance
//  1) If sign different determine boolean result based on sign
//  2) Otherwise determine boolean result based length of number amnd the sign
//  3) Same sign and same length. Do a core comparison and return the result
//
template <class _Ty> inline bool operator<( const _Ty& c, const int_precision& b )
	{
	int sign1, sign2, cmp;
	int_precision a(c);

	sign1 = a.sign();
	sign2 = b.sign();

	// Different signs
	if( sign1 > sign2 )
		return false;
	if( sign1 < sign2 )
		return true;

	// Same sign
	if( sign1 == 1 && a.size() < b.size() ) // Different therefore true
		return true;
	if( sign1 == 1 && a.size() > b.size() ) // Different therefore false
		return false;
	if( sign1 == -1 && a.size() > b.size() )
		return true;
	if( sign1 == -1 && a.size() < b.size() )
		return false;

	// Same sign and same length
	cmp = _int_precision_compare(const_cast<int_precision&>(a).pointer(), const_cast<int_precision&>(b).pointer());
	if( cmp < 0 && sign1 == 1 )
		return true;
	else
		if( cmp > 0 && sign1 == -1 )
			return true;

	return false;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/19/2005
//	@brief 		operator!=
//	@return 	bool	- 	the boolean value of the operator
//	@param		"a"	-	First operand number to compare
//	@param		"b"	-	Second operand number to
//
//	@todo
//
// Description:
//  Boolean not equal of two precision numbers
//
template <class _Ty> inline bool operator!=( int_precision& a, const _Ty& b )
	{
	return a == b ? false : true;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/19/2005
//	@brief 		operator!=
//	@return 	bool	- 	the boolean value of the operator
//	@param		"a"	-	First operand number to compare
//	@param		"b"	-	Second operand number to
//
//	@todo
//
// Description:
//  Boolean not equal of two precision numbers
//
template <class _Ty> inline bool operator!=( const _Ty& a, const int_precision& b )
	{
	return a == b ? false : true;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/19/2005
//	@brief 		operator>
//	@return 	bool	- 	the boolean value of the operator
//	@param		"a"	-	First operand number to compare
//	@param		"b"	-	Second operand number to
//
//	@todo
//
// Description:
//  Boolean greater of two precision numbers
//
template <class _Ty> inline bool operator>( int_precision& a, const _Ty& b )
	{
	return b < a ? true : false;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/19/2005
//	@brief 		operator>
//	@return 	bool	- 	the boolean value of the operator
//	@param		"a"	-	First operand number to compare
//	@param		"b"	-	Second operand number to
//
//	@todo
//
// Description:
//  Boolean greater of two precision numbers
//
template <class _Ty> inline bool operator>( const _Ty& a, const int_precision& b )
	{
	int_precision c(a);
	return b < c ? true : false;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/19/2005
//	@brief 		operator<=
//	@return 	bool	- 	the boolean value of the operator
//	@param		"a"	-	First operand number to compare
//	@param		"b"	-	Second operand number to
//
//	@todo
//
// Description:
//  Boolean less or equal of two precision numbers
//
template <class _Ty> inline bool operator<=( int_precision& a, const _Ty& b )
	{
	return b < a ? false : true;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/19/2005
//	@brief 		operator<=
//	@return 	bool	- 	the boolean value of the operator
//	@param		"a"	-	First operand number to compare
//	@param		"b"	-	Second operand number to
//
//	@todo
//
// Description:
//  Boolean less or equal of two precision numbers
//
template <class _Ty> inline bool operator<=( const _Ty& a, const int_precision& b )
	{
	int_precision c(a);
	return b < c ? false : true;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/19/2005
//	@brief 		operator>=
//	@return 	bool	- 	the boolean value of the operator
//	@param		"a"	-	First operand number to compare
//	@param		"b"	-	Second operand number to
//
//	@todo
//
// Description:
//  Boolean greater or equal of two precision numbers
//
template <class _Ty> inline bool operator>=( int_precision& a, const _Ty& b )
	{
	return a < b ? false: true;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/19/2005
//	@brief 		operator>=
//	@return 	bool	- 	the boolean value of the operator
//	@param		"a"	-	First operand number to compare
//	@param		"b"	-	Second operand number to
//
//	@todo
//
// Description:
//  Boolean less or equal of two precision numbers
//
template <class _Ty> inline bool operator>=( const _Ty& a, const int_precision& b )
	{
	return a < b ? false: true;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/Feb/2017, revised 20/JUL/2019
//	@brief 		gcd - Greatest Common Divisor
//	@return 	The greates common divisor or a & b
//	@param		"a"	-	First operand number 
//	@param		"b"	-	Second operand number
//
//	@todo
//
// Description:
//  gcd of two integer. Tis should work for both signed and unsigned operands
//  change the while loop while(b>0) to while(b!=0) to accomodate negative b
//
template<class _Ty> inline _Ty gcd(const _Ty lhs, const _Ty rhs)
	{
	_Ty tmp, a = lhs, b = rhs;
	// GCD(0,rhs)==rhs; GCD(lhs,0)==0; GCD(0,0)==0
	if (a == (_Ty)0) return b;
	if (b == (_Ty)0) return a;
	while (b !=(_Ty)0) { tmp = b; b = a%b; a = tmp; }
	return a;
	}

//////////////////////////////////////////////////////////////////////
//
//
//    Logical !,&&,||
//
//
//////////////////////////////////////////////////////////////////////

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/Nov/2021
//	@brief 		operator unary !(Logical NOT)
//	@return 	int_precision	-	-a
//	@param		"a"	-	operand for ! operator
//
//	@todo
//
// Description:
//  Unary ~ operator
//  Negate Integer
//
inline bool operator!(const int_precision& a)
	{
	if (a.iszero()) return true;
	return false;
 	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/Nov/2021
//	@brief 		operator unary && (Logical and)
//	@return 	int_precision	-	-a
//	@param		"a"	-	operand for && operator
//
//	@todo
//
// Description:
//  Unary ~ operator
//  Negate Integer
//
inline bool operator&&(const int_precision& a, const int_precision& b)
	{
	if (a.iszero() || b.iszero()) return false;
	return true;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/Nov/2021
//	@brief 		operator unary || (Logical or)
//	@return 	int_precision	-	-a
//	@param		"a"	-	operand for || operator
//
//	@todo
//
// Description:
//  Unary ~ operator
//  Negate Integer
//
inline bool operator||(const int_precision& a, const int_precision& b)
	{
	if (a.iszero() && b.iszero()) return false;
	return true;
	}

#endif