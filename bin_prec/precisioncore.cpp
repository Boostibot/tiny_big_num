/*
 *******************************************************************************
 *
 *
 *                       Copyright (c) 2007-2022
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
 * Module name     :precisioncore.cpp
 * Module ID Nbr   :   
 * Description     :Arbitrary precision core functions for integer and floating
 *					point precision class
 * -----------------------------------------------------------------------------
 * Change Record   :   
 *
 * Version	Author/Date		Description of changes
 * -------  ---------------	----------------------------------------------------
 * 03.01	HVE/14-Aug-2021	Switch to new internal binary format iptype & fptype
 * 03.02	HVE/19-Oct-2021	Passing all internal testing
 * 03.03	HVE/19-Nov-2021	Fixed compiler bugs reported by GNU version on a mac
 * 03.04	HVE/20-Nov-2021 A few bugs fixed and change to avoid compiler warnings
 * 03.05	HVE/21-Nov-2021 More cleaning up and improvements
 * 03.06	HVE/22-Nov-2021	Fix a bug in _float_precision_fptoa() where fraction was incorret truncated
 *							Also fix a bug in the modf().
 * 03.07	HVE/26-Nov-2021 Added 1/sqrt(2) and sqrt(2) to the "constant" _float_table
 * 03.08	HVE/5-Dec-2021	Improved the algorithm for float_precion_ftoa() that improved performanen for a 100K digits varianble
 *							from more than 21000sec -> 300msec.
 * 03.09	HVE/6-Dec-2021	Fix an float_precision issue where numbers with small precision but high negative exponent was converted to 
 *							0E0 instead of a very small number. Also added gigatrunk splitting in float_precision_fptoa 
 * 03.10	HVE/8-Dec-2021	Change eptype to an intmax_t insead of int to raise the exponent limit to more than 300M digits
 * 03.11	HVE/10-Dec-2021	Added two more trunking levels for faster handling of digits exceedig 10-100M decimal digits
 * 03.12	HVE/11-Dec-2021	Allowed separaors ' symbols to be part of a string based number for atoip() and atofp(). atofp() can be called with a std::string or a char *
 * 03.13	HVE/25-Dec-2021	remove of double *a,*b with new and replaced it with vector<double> va, vb. 
 * 03.14	HVE/25-Dec-2021 Added _float_precision_schonhagen_strassen_umul() to add better multiplication for medium size multiplication of digits<6,000 the function was modified 
 *							from the int_precision counterpart
 * 03.15	HVE/26-Dec-2021	Change the FFT _int_precision_umul_fourier() and _float_precision_umul_fourier() to be able to handle digits in excess of 25E9 digits 
 * 03.16	HVE/29-Dec-2021 Fixed a few smaller bugs in _flot_precision_fptoa()
 * 03.17	HVE/3-Jan-2022	Special test version plus a fix in _umul_fourier() where the 4bit version was kicked in aftter 18M digits instead of >150M digits
 * 03.18	HVE/3-Jan-2022	A bug was found in the float_precision_schonhage-straasen_umul() and was change to just a call the _umul_fourier() as a temporary fix
 * 03.19	HVE/6-Jan-2022	Fixed the overflow bug in _int_precision_umul_fourier() & _float_precision_umul_fourier() that happens only for very large multiplications > 10M digits
 * 03.20	HVE/6-Jan-2022	Added _INVSQRT3 and SQRT3 as a build in constant
 * 03.21	HVE/8-Jan-2022	Added Multi threading in _int_multiplication_umul_fourier(), _float_precision_umul_fourier() and float_table(_PI)
 * 03.22	HVE/10-Jan-2022	Minor optimazation when multiply by power of 2. It is faster to do x.exponent(x.exponent()+'power of two');
 * 03.23	HVE/19-Jan-2022	Change name of schoonhage-strassen to just _umul_linear to followed the name convention. change _umul to _umul_school nd added a new _umul as the entry point 
 *							for all multiplication algorithm
 * 03.24	HVE/21-Jan-2022	Use string2number in the _float_precision_atofp for higher performance, particular when digits exceed 100,000+ digits
 * 03.25	HVE/4-Feb-2022	Added method adjustExponent() for a faster way to multiply or divide by any power of 2
 * 03.26	HVE/8-Feb-2022	Fxied an bug in _int_precision_urem() that incorrectly called _int_precision_urem_short() when short argument was bigger than 32bit. Also
 *							fxied an issue in _int_precision_urem_short() where shortcut was not working
 * 03.27	HVE/9-Feb-2022	_int_precision_itoa() has improved performance for base 2, 8 & 16
 * 03.28	HVE/14-Mar-2022	Fixed a bug in _int_precision_umul_fourier for int multiplication exceeding 153M digits. fixed another bug in the same function for exceeding the accuracy of a double
 * 03.29	HVE/16-Mar-2022 Fixed a bug in _int_precision_itoa() for base==2, 8 or 16
 * 03.30	HVE/19-Mar-2022	Removed debug code that was left unattended
 * 03.31	HVE/25-Mar-2022	Fixed a issues in _int_precision_atoip() for base=2,8 and 16
 * 03.32	HVE/13-Apr-2022	Change _float_precision_inverse to add the Brent Improvement resulting in slightly better performance
 * 03.33	HVE/19-Apr-2022	Change _float_precision_inverse and sqrt() to now use Halley 3rd order method that is approx 20% faster than previous Newton (2nd order) method.
 * 03.34	HVE/10-May-2022	The sinh(), cosh() has been added automatic argument reduction (instead of fixed argument reduction) and coefficients scalling, exp() has also been updated
 * 03.35	HVE/20-May-2022 Added log2(), and support both an enhanced version of log() using Taylor series and a version using AGM. the mainfunction log() will detrmine which of these two to call
 * 03.36	HVE/21-May-2022 Added AGM() using float_precision parameters (Arithmetic-Geometric Mean)
 * 03.37	HVE/24-Jun-2022	The Trigonometric functions sin(x), cos(x), tan(x), asin(x), acos(x), atan(x) has been improved and adding more aggressive reductions factors and Taylor terms grouping
 * 03.38	HVE/14-Jul-2022	Improved sinh() and cosh() including increase coefficient scaling
 * 03.39	HVE/17-Jul-2022 Change the algorithm for calculating e to the binary splitting method resulting in a speedup of approx 10 times
 * 03.40	HVE/19-Aug-2022	Change the PI method to the Chudnovsky binary splitting, which is 8 tims faster than the Brent-Salamin method it replaced
 * 03.41	HVE/5-Sep-2022	Restore the _int_precision_compare() to the previous pointer arguments
 * 03.42	HVE/8-Sep-2022	Added multi threading for calculation of PI
 * 03.43	HVE/9-Sep-2022	Added four ways multi threading to calculation of e
 * 03.44	HVE/20-Sep-2022 Enhanced sqrt() and _float_precision_inverse() by using a hybrid algorithm that automatically determine if 2nd order or 3rd order is the most efficient algorithm
 *							Improved the termination criteria for sqrt(), _float_precision_inverse(), _INVSQRT2, _INVSQRT3
 *
 * End of Change Record
 * --------------------------------------------------------------------------
*/

/* define version string */
static char _VIP_[] = "@(#)precisioncore.cpp 03.44 -- Copyright (C) Henrik Vestermark";

#include <cstdint>
#include <ctime>
#include <cmath> 
#include <iostream>
#include <iomanip>
#include <cstring>
#include <vector>
#include <algorithm>
#include <thread>

#include "iprecision.h"
#include "fprecision.h"

// #define HVE_DEBUG
#define HVE_THREAD				// Add Multi threading for _int_precison_umul_fourier(), _float_precison_umul_fourier() and PI calculations


///////////////////////////////////////////////
//
//
//    Integer Precision Core
//
//
///////////////////////////////////////////////


///////////////////////////////////////////////
//
//
//    Integer Precision Input and output operator
//
//
///////////////////////////////////////////////

std::ostream& operator<<( std::ostream& strm, const int_precision& d ) 
	{
	return strm << _int_precision_itoa(const_cast<int_precision *>(&d)).c_str();
	}

std::istream& operator>>( std::istream& strm, int_precision& d )
         { 
         char ch; std::string s;
         strm.get(ch);// strm >> ch; 
         while( ch == ' ' ) strm.get(ch);  // Ignore leading white space.
         if( ch == '+' || ch == '-' ) { s += ch; strm.get(ch); } else s += '+';  // Parse sign

         if( ch == '0' ) // Octal, Binary or Hexadecimal number
            {
            strm.get( ch );
            if( ch == 'x' || ch == 'X' ) // Parse Hexadecimal
               for( s += "0x"; (ch >= '0' && ch <= '9') || (ch >='a' && ch <= 'f') || (ch >= 'A' && ch <= 'F'); strm.get( ch ) ) s += ch;
            else
               if( ch == 'b' || ch == 'B' )  // Parse Binary
                  for( s += "0b"; ch >= '0' && ch <= '1'; strm.get( ch ) ) s += ch;
               else // Parse Octal
                  for( s += "0"; ch >= '0' && ch <= '7'; strm.get( ch ) ) s += ch;
            }
         else // Parse Decimal number
            for( ; ch >= '0' && ch <= '9'; strm.get(ch) /*strm >> ch*/ ) s += ch;

         strm.putback( ch );  // ch contains the first character not part of the number, so put it back
         if(!strm.fail() && s.length() >= 2 )  // Valid number has at least a length of 2 or higher
            d = int_precision( const_cast<char *>( s.c_str() ) );

         return strm;
         }


///////////////////////////////////////////////
//
//
//    Miscellaneous
//
//
///////////////////////////////////////////////

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  1/19/2005
//	@brief 		std::string _int_precision_itoa Convert number to ascii string
//	@return		std::string -	the converted number in ascii string format
//	@param		"a"	-	Number to convert to ascii
//	@param		"base"	-	base for conversion to ascii
//
//	@todo 
//
// Description:
//   Convert int_precsion to ascii string
//   using base. Sign is only added if negative
//	  Base is default BASE_10 but can be anything from base 2..36
//
std::string _int_precision_itoa( int_precision *a, const int base )
   {
   return ( a->sign() < 0 ? "-": "" ) + _int_precision_itoa( a->pointer(), base );
   }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date  1/19/2005
//	@brief 		std::string _int_precision_itoa Convert number to ascii string
//	@return		std::string -	the converted number in ascii string format
//	@param		"a"	-	Number to convert to ascii
//	@param		"base"	-	base for conversion to ascii
//
//	@todo 
//
// Description:
//   Convert int_precsion to ascii string
//   using base. Sign is only added if negative
//	  Base is default BASE_10 but can be anything from base 2..36
//
std::string _int_precision_itoa( int_precision& a, const int base)
	{
	return (a.sign() < 0 ? "-" : "") + _int_precision_itoa(a.pointer(), base);
	}


///////////////////////////////////////////////
//
//
//    Core Support Functions. 
//
//
///////////////////////////////////////////////

//
// Core functions
// The core functions all perform unsigned arithmetic un elements of the string class!
//    _int_precision_strip_leading_zeros	-- Strips non significant leading zeros
//    _int_precision_uadd_short			-- add a short digit [0..RADIX] to the string
//    _vector_reverse_binary			-- Reverse bit in the data buffer
//    _vector_fourier					-- Fourier transformn the data
//    _vector_real _fourier				 -- Convert n discrete double data into a fourier transform data set
//

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  1/19/2005
//	@brief 	_int_reverse_binary
//	@return 	void	-	
//	@param   "data[]"	-	array of double complex number to permute
//	@param   "n"	-	number of element in data[]
//
//	@todo  
//
// Description:
//   Reverse binary permute
//   n must be a power of 2
//
static void _vector_reverse_binary(std::complex<double> data[], const size_t n)
	{
	size_t i, j, m;

	if (n <= 2) return;

	for (j = 1, i = 1; i < n; i++)
		{
		if (j > i)
			{
			//if (j - 1 >= n || i - 1 >= n)  // Debug
			//	j = j;// DEBUG Error
			std::swap(data[j - 1], data[i - 1]);
			}

		for (m = n >> 1; m >= 2 && j > m; m >>= 1)
			j -= m;

		j += m;
		}
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  1/19/2005
//	@brief 	_int_fourier do the fourier transformation
//	@return 	static void	-	
//	@param   "data[]"	-	complex<double> fourie data
//	@param   "n"	-	number of element in data (must be a power of 2)
//	@param   "isign"	-	transform in(1) or out(-1)
//
//	@todo
//
// Description:
//   Wk=exp(2* PI *i *j )  j=0..n/2-1
//   exp( k * i *j ) => exp( k * i * (j-1) + k * i ) => exp( t + o ) for each new step
//   exp( t + o ) => exp(t)-exp(t)*( 1 - cos(o) -isin(o) ) => exp(t)-exp(t)*(a-ib)
//   => exp(t)+exp(t)*(-a+ib) => exp(t)( 1 + (-a+b) )
//   sin(t+o)=sin(t)+[-a*sin(t)+b*cos(t)]
//   a=2sin^2(o/2), b=sin(o)
//   n must be a power of 2
//
static void _vector_fourier(std::complex<double> data[], const size_t n, const int isign)
	{
	double theta;
	std::complex<double> w, wp;
	size_t mh, m, r, j, i;

	_vector_reverse_binary(data, n);

	for (m = 2; n >= m; m <<= 1)
		{
		theta = isign * 2 * 3.14159265358979323846264 / m;
		wp = std::complex<double>(-2.0 * sin(0.5 * theta) * sin(0.5 * theta), sin(theta));
		w = std::complex<double>(1, 0); // exp(0) == exp( isign*2*PI*i/mmax * m-1 )
		mh = m >> 1;

		for (j = 0; j < mh; j++)      // m/2 iteration
			{
			for (r = 0; r <= n - m; r += m)
				{
				std::complex<double> tempc;
				i = r + j;
				//if (i + mh >= n) //DEBUG
				//	i = i;  // DEBUG
				tempc = w * data[i + mh];              // u=data[i]; v=data[j]*w; data[i]=u+v;data[j]=u-v;
				data[i + mh] = data[i] - tempc;
				data[i] += tempc;
				}

			w = w * wp + w;  // w = w(1+wp) ==> w *=1+wp;
			}
		}
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  1/19/2005
//	@brief 			_int_real_fourier
//	@return 			static void	-	
//	@param   "data[]"	-	
//	@param   "n"	-	number of data element in data. n must be a power of 2)
//	@param   "isign"	-	Converting in(1) or out(-1)
//
//	@todo  
//
// Description:
//   Convert n discrete double data into a fourier transform data set
//   n must be a power of 2
//
static void _vector_real_fourier(std::vector<double>& data, const size_t n, const int isign)
	{
	size_t i;
	double theta, c1 = 0.5, c2;
	std::complex<double> w, wp, h1, h2;
	double *ptr = data.data() ;

	theta = 3.14159265358979323846264 / (double)(n >> 1);
	if (isign == 1)
		{
		c2 = -c1;
		_vector_fourier((std::complex<double> *)ptr, n>>1, 1);  // n>> 1 == complex numbers
		}
	else
		{
		c2 = c1;
		theta = -theta;
		}
	wp = std::complex<double>(-2.0 * sin(0.5 * theta) * sin(0.5 * theta), sin(theta));
	w = std::complex<double>(1 + wp.real(), wp.imag());
	for (i = 1; i < (n >> 2); i++)
		{
		size_t i1, i2, i3, i4;
		std::complex<double> tc;

		i1 = i + i;
		i2 = i1 + 1;
		i3 = n + 1 - i2;
		i4 = i3 + 1;
		h1 = std::complex<double>(c1 * (data[i1] + data[i3]), c1 * (data[i2] - data[i4]));
		h2 = std::complex<double>(-c2 * (data[i2] + data[i4]), c2 * (data[i1] - data[i3]));
		tc = w * h2;
		data[i1] = h1.real() + tc.real();
		data[i2] = h1.imag() + tc.imag();
		data[i3] = h1.real() - tc.real();
		data[i4] = -h1.imag() + tc.imag();
		w *= (std::complex<double>(1) + wp);
		}
	if (isign == 1)
		{
		double t;
		data[0] = (t = data[0]) + data[1];
		data[1] = t - data[1];
		}
	else
		{
		double t;
		data[0] = c1*((t = data[0]) + data[1]);
		data[1] = c1*(t - data[1]);
		_vector_fourier((std::complex<double> *)ptr, n>>1, -1);
		}
	}


///////////////////////////////////////////////
//
//
//    Core Functions. BINARY
//
// The core functions all perform unsigned arithmetic of elements of the vector<iptype> class!
//
//	  build_i_number					--	Build integer number as vector<iptype>
//    _int_precision_strip_leading_zeros	-- Strips non significant leading zeros
//    _int_precision_strip_trailing_zeros	-- Strips non significant trailing zeros
//	  _int_precision_clz				-- Count leading zeros in an iptype
//	  _int_precision_clz				-- Count leading zeros in an vector<iptype>
//	  _int_precision_ctz				-- Count trailing zeros in an iptype
//	  _int_precision_ctz				-- Count trailing zeros in an vector<iptype>
//	  _int_precision_csb				-- Bit position of the most significant bit in vector<iptype>
//    _int_precision_compare			-- Compare two strings for numeric order
//    _int_precision_uadd_short			-- Add a short iptype digit e.g. [0..2^64] to the number
//    _int_precision_uadd				-- Add two unsigned binary numbers
//    _int_precision_usub_short			-- Subtract a short iptype digit [0..2^64] from the number
//    _int_precision_usub				-- Subtract two unsigned binary numbers
//    _int_precision_umul_short			-- Multiply a iptype digit [0..2^64] to the number
//	  _int_precision_umul				-- Multiply two unsigned binary numbers using the most optimal multiplication algorithm
//    _int_precision_umul_school		-- Multiply two unsigned binary numbers using old fashion school algorithm
//    _int_precision_udiv_short			-- Divide a iptype digit [0..^64] up in the number
//    _int_precision_udiv				-- Divide two unsinged strings
//    _int_precision_urem				-- Remainder of dividing two unsigned numbers
//	  _int_precision_urem_short			-- Get the remainder an iptype digit[0.. ^ 64] up in the number
//	  _int_precision_shiftright			-- Shift right the vector<iptype> number
//	  _int_precision_shiftleft			-- Shift left the vector<iptype> number
//    _int_precision_itoa				-- Convert internal precision to BASE_10 string
//	  _int_precision_uand				-- And the binary numbers together
//	  _int_precision_uor				-- Or the binary numbers together
//	  _int_precision_xor				-- Xor the beinary numbers together
//	  _int_precision_unegate			-- Negate the binary number
//    _build_i_number					-- Build the binary number from string representation
//	
//	  _vector_reverse_binary			-- Reverse bit in the data buffer
//    _vector_fourier					-- Fourier transformn the data
//    _vector_real_fourier				-- Convert n discrete double data into a fourier transform data set
//	  _precision_umul64					-- Generic multiplication of two vector<iptype> or vector<fptype> numbers
//	  _int_precision_umul_fourier		-- Multiply two binary numbers using FFT
//	  _int_precision_umul_karatsuba		-- Multiply two unsigned binary numbers using karatsuba algorithm
//	  _int_precision_umul_linear		-- Multiply two unsigned binary numbers using schonhagen-strassen algorithm with linearconvolution
//
//		decimal2number
//		trunk2number
//		string2number
//		stringbase2number
//		stringbase8_2number
//		_int_precision_iptod			- Convert int_precision number to double
//		
//
///////////////////////////////////////////////

// Use for various conversions to and from strings
static size_t _powerof10Table[20] = { 1,10,100,1000,10 * 1000,100 * 1000,1000 * 1000,10 * 1000000,100 * 1000000,1000 * 1000000,
									10000ull * 1000000,100000ull * 1000000,1000000ull * 1000000,10ull * 1000000 * 1000000,100ull * 1000000 * 1000000,
									1000ull * 1000000 * 1000000, 10000ull * 1000000 * 1000000, 100000ull * 1000000 * 1000000,1000000ull * 1000000 * 1000000,
									1000'000ull * 1000'000 * 1000'000 * 10 };


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief		Build a binary repesentation of signle digit number
//	@return		std::vector<iptype>	- The integer precision string	
//	@param		"digit"		- The next digit to be added to the integer point number to convert. Is always positive
// @param		"base"		- The base of the digit being added
//
//	@todo 	
//
// Description:
//   Add a digit to the number being build for the integer precision number
//   The function dosnt create any leading significant zeros
//   To run it efficiently is is better to take advantages of iptype (64bit) instead of just one decimal,
//		binary,octal or hexdecimal digit at a time
//    
static inline std::vector<iptype> build_i_number(std::vector<iptype> &number, iptype digit, iptype base)
	{
	number = _int_precision_umul_short(&number, base);
	number = _int_precision_uadd_short(&number, digit);
	return number;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/Dec/2021
//	@brief		Build a binary repesentation of a string of single digit number
//	@return		std::vector<iptype>	- The integer precision string
//  @param		"s"			- The decimal string 
//	@param		"start"		- The start index of the string 
//  @param		"end"		- The End index of the string
//
//	@todo 	
//
// Description:
//	Build and collect th binary number that correspond to the decimal representation of the number
//
static std::vector<iptype> decimal2number(std::vector<iptype>& number, std::string s, size_t start, size_t end)
	{
	size_t i=start, length=end-start;
	const size_t max_digits = 19;
	iptype pwr;

	for(i=start;length>=max_digits;length-=max_digits, i+=max_digits)
		{
		std::string s2 = s.substr(i, max_digits);
		uint64_t n = strtoull(s2.c_str(), NULL, BASE_10);
		pwr = _powerof10Table[max_digits];
		build_i_number(number, n, pwr);
		}
	
	if (length!= 0)
		{
		std::string s2 = s.substr(i, length);
		uint64_t n = strtoull(s2.c_str(), NULL, BASE_10);
		pwr = _powerof10Table[length];
		build_i_number(number, n, pwr);
		}

	return number;
	}

// Do  a single trunk and return it
// It is guarantee that there is trunk size of data
static std::vector<iptype> trunk2number(std::vector<iptype>& number, std::string s, size_t start)
	{
	size_t i;
	std::vector<iptype> sum(1,0);
	size_t thr = MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS;
	static std::vector<iptype> _trunkPowerof10;

	if ( _trunkPowerof10.size()==0) // is _trunkPowerof10 build or created
		{
		std::vector<iptype> p(1,_powerof10Table[MAX_DECIMAL_DIGITS]);
		for (i = MAX_TRUNK_SIZE, _trunkPowerof10.assign(1,1); i > 0; i >>= 1)
			{// Build multiply factor for trunk size
			if (i & 1) _trunkPowerof10 =_int_precision_umul(&_trunkPowerof10, &p);	// Odd
			if (i > 1) p = _int_precision_umul(&p, &p);					// square it
			}
		}

	decimal2number(sum, s, start, start + thr);
	number =_int_precision_umul(&number, &_trunkPowerof10);
	number =_int_precision_uadd(&number, &sum );

	return number;
	}

#if false
// Do a single kilo trunk and return it
// It is guarantee that there is trunk size of data
static std::vector<iptype> kilo2number(std::vector<iptype>& number, std::string s, size_t start)
	{
	size_t i;
	std::vector<iptype> sum(1, 0);
	size_t thr = MAX_KILOTRUNK_SIZE;
	static std::vector<iptype> _trunkPowerof10;

	if (_trunkPowerof10.size() == 0) // is _trunkPowerof10 build or created
		{
		std::vector<iptype> p(1, _powerof10Table[MAX_DECIMAL_DIGITS]);
		for (i = MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE, _trunkPowerof10.assign(1, 1); i > 0; i >>= 1)
			{// Build multiply factor for trunk size
			if (i & 1) _trunkPowerof10 = _int_precision_umul_fourier(&_trunkPowerof10, &p);	// Odd
			if (i > 1) p = _int_precision_umul_fourier(&p, &p);					// square it
			}
		}

	for (i=start;thr>0; --thr, i+=MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS)
		{
		trunk2number(sum, s, i);
		number = _int_precision_umul_fourier(&number, &_trunkPowerof10);
		number = _int_precision_uadd(&number, &sum);
		sum.assign(1, 0);
		}

	return number;
	}
#endif


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Jan/2022
//	@brief		Build a binary repesentation of a string of single digit decimal numbers
//	@return		std::vector<iptype>	- The integer precision number
//  @param		"s"			- The decimal string 
//	@param		"start"		- The start index of the string 
//  @param		"end"		- The End index of the string
//
//	@todo 	
//
// Description:
//	Build and collect the binary number that correspond to the decimal representation of the number
//
static std::vector<iptype> string2number(const std::string s, const size_t start, size_t end)
	{
	size_t i = end - start, j=0, radix_inx=0;
	std::string s2;
	std::vector<std::vector<iptype> > vn(0);
	std::vector<iptype> num(1);
	static std::vector<std::vector<iptype> > radix;
	
	vn.reserve(i/MAX_DECIMAL_DIGITS+16);
	// Step 1 partition the string into a binary vector with 1 binary digit in order of least to most significant 
	for (; i > MAX_DECIMAL_DIGITS; i -= MAX_DECIMAL_DIGITS, end-= MAX_DECIMAL_DIGITS)
		{
		s2 = s.substr(end - MAX_DECIMAL_DIGITS, MAX_DECIMAL_DIGITS);
		num[0] = strtoull(s2.c_str(), NULL, BASE_10);
		vn.push_back(num);
		}
	s2 = s.substr(start, i);
	num[0] = strtoull(s2.c_str(), NULL, BASE_10);
	vn.push_back(num);
	
	// Step2 collected into higher binary values by reducing the vector with 2,3,...,n MAX_DECIMAL_DIGITS
	if (radix.size() == 0)
		{
		num[0] = _powerof10Table[MAX_DECIMAL_DIGITS];
		radix.push_back(num);
		}
	for (;vn.size() > 1; ++radix_inx )
	{
		if (radix_inx >= radix.size())
		{
			num = _int_precision_umul(&radix[radix_inx - 1], &radix[radix_inx - 1]); // replace by _int_precision_square_fourier() when ready
			radix.push_back(num);  
		}
		for (i = 0, j = 0; j < vn.size(); ++i, j += 2)
		{
			if (j + 1 < vn.size())
			{
				num = _int_precision_umul(&vn[j + 1], &radix[radix_inx]);
				vn[i] = _int_precision_uadd(&vn[j], &num);
			}
			else
				vn[i] = vn[j];
		}
		vn.resize(i);
//		std::cout << "\tResize vn =" << vn.size() << std::endl;  // DEBUG HVE
	}

	return vn[0];
 	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Jan/2022
//	@brief		Build a binary repesentation of a string of single digit numbers in base  or base 16
//	@return		std::vector<iptype>	- The integer precision number
//  @param		"s"			- The decimal string 
//	@param		"start"		- The start index of the string 
//  @param		"end"		- The End index of the string
//	@param		"base"		- The base of the string number (either base 2 or base 16)
//
//	@todo 	
//
// Description:
//	Build and collect the binary number that correspond to the decimal representation of the number
//	Base is either base 2 or base 16
//
static std::vector<iptype> stringbase2number(const std::string s, const size_t start, size_t end, int base)
	{
	const uintmax_t baselength = base == BASE_2 ? 64 : 16;
	size_t i = end - start;
	std::string s2;
	std::vector<iptype> vn;
	iptype n;

	vn.reserve(i / baselength + 16);
	// Step 1 partition the string into a binary vector with 1 binary digit in order of least to most significant 
	for (; i > baselength; i -= baselength, end -= baselength)
		{
		s2 = s.substr(end - baselength, baselength);
		n = strtoull(s2.c_str(), NULL, base);
		vn.push_back(n);
		}
	s2 = s.substr(start, i);
	n= strtoull(s2.c_str(), NULL, base);
	vn.push_back(n);

	return vn;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		22/Jan/2022
//	@brief		Build a binary repesentation of a string of single digit numbers in base 8 (octal)
//	@return		std::vector<iptype>	- The integer precision number
//  @param		"s"			- The decimal string 
//	@param		"start"		- The start index of the string 
//  @param		"end"		- The End index of the string
//
//	@todo 	
//
// Description:
//	Build and collect the binary number that correspond to the decimal representation of the number
//	Currently Only base 8 is supproted
//
static std::vector<iptype> stringbase8_2number(const std::string s, const size_t start, size_t end)
	{
	size_t i = end - start, j, radix_inx = 0;
	std::string s2;
	std::vector<std::vector<iptype> > vn(0);
	std::vector<iptype> num(1);
	static std::vector<std::vector<iptype> > radix;

	vn.reserve(i / MAX_OCTAL_DIGITS + 16);
	// Step 1 partition the string into a binary vector with 1 binary digit in order of least to most significant 
	for (; i > MAX_OCTAL_DIGITS; i -= MAX_OCTAL_DIGITS, end -= MAX_OCTAL_DIGITS)
		{
		s2 = s.substr(end - MAX_OCTAL_DIGITS, MAX_OCTAL_DIGITS);
		num[0] = strtoull(s2.c_str(), NULL, BASE_8);
		vn.push_back(num);
		}
	s2 = s.substr(start, i);
	num[0] = strtoull(s2.c_str(), NULL, BASE_8);
	vn.push_back(num);

	// Step2 collected into higher binary values by reducing the vector with 2,3,...,n MAX_OCTAL_DIGITS
	if (radix.size() == 0)
		{
		num[0] = 01'000'000'000'000'000'000'000ull;  // MAX_OCTAL DIGITS cant be changed without changing th radix constant
		radix.push_back(num);
		}
	for (; vn.size() > 1; ++radix_inx)
		{
		if (radix_inx >= radix.size())
			{
			num = _int_precision_umul(&radix[radix_inx - 1], &radix[radix_inx - 1]); // replace by _int_precision_square_fourier() when ready
			radix.push_back(num);
			}
		for (i = 0, j = 0; j < vn.size(); ++i, j += 2)
			{
			if (j + 1 < vn.size())
				{
				num = _int_precision_umul(&vn[j + 1], &radix[radix_inx]);
				vn[i] = _int_precision_uadd(&vn[j], &num);
				}
			else
				vn[i] = vn[j];
			}
		vn.resize(i);
		}

	return vn[0];
	}

// @author Henrik Vestermark(hve@hvks.com)
//	@date		26/Aug/2022
//	@brief 		Convert inp_precision to double (IEE754) 
//	@return 	double		- The converted int_precision number
//
//	@todo 	
//
// Description:
//   Convert int_precision into a double number
//	ip can be huge but double can only contain 52bit mantissa with an implicit 1 bit in front.
//  First convert iprecision -> float_precision
//	Second Convert float_precision -> double
//
double _int_precision_iptod(const int_precision* ip)
	{
	uintmax_t t = 0, expo;
	double d;
	std::vector<fptype> fnum(0);
	const std::vector<iptype> inum(ip->number());

	if (ip->iszero())
		return 0.0;
	for (size_t i = inum.size(); i > 0; --i)
		fnum.push_back((fptype)inum[i - 1]);
	expo = _float_precision_normalize(&fnum);
	expo += (eptype)((inum.size() - 1) * Bitsfptype);
	expo += _float_precision_rounding(&fnum, ip->sign(), 16, ROUND_NEAR);// IEEE754 double can hold 15.955 decimal digit
	expo += 1023;  // Biased double exponent format 
	if (fnum.size()>1)
		t = (uintmax_t)fnum[1];
	t >>= 12;				// Make room for exponent
	t |= (expo << 52) & 0x7fffffffffffffff; // Add exponent
	d = *(double *)&t;		// Make it a double
	if (ip->sign() < 0)		// Check sign
		d = -d;				// Make it negative

	return d;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		4/Aug/2021
//	@brief 		_int_precision_strip_leading_zeros
//	@return		void	-	
//	@param		"s"	-	pointer to source operand
//
//	@todo
//
// Description:
//   Remove leading nosignificant zeros of the binary number
//	This is from the start of the vector<iptype>
//
void _int_precision_strip_leading_zeros(std::vector<iptype> *s)
	{
	std::vector<iptype>::iterator pos;

	// Strip leading zeros
	for (pos = s->begin(); pos != s->end() && *pos == (iptype)0; ++pos);	// Find first not zero digit

	if (s->begin() != pos)
		s->erase(s->begin(), pos);
	if (s->empty())
		s->assign(1,(iptype)0);

	return;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		13/Aug/2021
//	@brief 		_int_precision_strip_trailing_zeros
//	@return		void	-	
//	@param		"s"	-	pointer to source operand
//
//	@todo
//
// Description:
//   Remove trailing nosignificant zeros of the binary number
//		this is from the top of the vector<iptype> 
//
void _int_precision_strip_trailing_zeros(std::vector<iptype> *s)
	{
	size_t i;
	std::vector<iptype>::reverse_iterator pos;

	// Strip leading zeros
	for (i = s->size() - 1, pos=s->rbegin(); i > 0 && *pos == (iptype)0; ++pos, --i);

	s->resize(i + 1);  // Keep at least one digit by default

	return;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		13/Aug/2021
//	@brief 		_int_precision_clz
//	@return		unsigned int	-	the count of leading zero bits in "a"
//	@param		"a"	-	iptype operand
//
//	@todo
//
// Description:
//   Count leading nosignificant zeros of the binary iptype 
//
size_t _int_precision_clz(const iptype a)
	{
	iptype x = a;
	size_t offset = 0;
	static const unsigned char lookup[16] = { 4,3,2,2,1,1,1,1,0,0,0,0,0,0 };

	if (sizeof(iptype) > 4 && x & 0xffffffff00000000u)
		x >>= 32; else offset += 32;

	if (sizeof(iptype) > 2 && x & 0xffff0000u)
		x >>= 16; else offset += 16;

	if (sizeof(iptype) > 1 && x & 0xff00u)
		x >>= 8; else offset += 8;

	if ( x & 0xf0u)
		x >>= 4; else offset += 4;
	offset += lookup[(unsigned char)x];
	return offset;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Sep/2021
//	@brief 		_int_precision_clz
//	@return		unsigned int	-	the count of leading zero bits in "a"
//	@param		"mb"	-	vector<iptype> operand
//
//	@todo
//
// Description:
//   Count leading nosignificant zeros of the binary iptype 
//
size_t _int_precision_clz(const std::vector<iptype> &mb ) 
	{
	size_t tot_cnt = 0, cnt;
	for (size_t i = mb.size(); i > 0; --i)
		{
		cnt = _int_precision_clz(mb[i - 1]);
		tot_cnt += cnt;
		if (cnt != 64) break;
		}
	return tot_cnt;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		13/Aug/2021
//	@brief 		_int_precision_ctz
//	@return		unsigned int	-	the count of trailing zero bits in "a"
//	@param		"a"	-	iptype operand
//
//	@todo
//
// Description:
//   Count trailing nosignificant zeros of the binary iptype 
//	  iptype bit word input to count zero bits on right
//   cnt will be the number of zero bits on the right,
//   so if a is 1101000 (base 2), then c will be 3
// NOTE: if 0 == a, then c = 64.
//
size_t _int_precision_ctz(const iptype a)
	{
	iptype x = a;  // sizeof(iptype) can be 8, 4, 2, or 1 
	size_t cnt;
	if (x & 0x1) 
		cnt = 0;
	else
		{
		cnt = 1;
		if (sizeof(iptype) >= 8 && (x & 0xffffffffu) == 0)
			{
			x >>= 32; cnt += 32;
			}
		if (sizeof(iptype) >= 4 && (x & 0xffffu) == 0)
			{
			x >>= 16; cnt += 16;
			}
		if (sizeof(iptype) >= 2 && (x & 0xff) == 0)
			{
			x >>= 8; cnt += 8;
			}
		if ((x & 0xf) == 0)
			{
			x >>= 4; cnt += 4;
			}
		if ((x & 0x3) == 0)
			{
			x >>= 2; cnt += 2;
			}
		cnt -= x & 0x1;
		}
	
	return cnt;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Sep/2021
//	@brief 		_int_precision_ctz
//	@return		unsigned int	-	the count of trailing zero bits in "a"
//	@param		"mb"	-	vector<iptype> operand
//
//	@todo
//
// Description:
//   Count trailing nosignificant zeros of the binary iptype 
//	  iptype bit word input to count zero bits on right
//   cnt will be the number of zero bits on the right,
//   so if a is 1101000 (base 2), then c will be 3
// NOTE: if 0 == a, then c = 64.
//
size_t _int_precision_ctz( const std::vector<iptype> &mb )
	{
	size_t tot_cnt = 0, cnt;
	for (size_t i = 0; i < mb.size(); ++i)
		{
		cnt = _int_precision_ctz(mb[i]);
		tot_cnt += cnt;
		if (cnt != 64) break;
		}
	return tot_cnt;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		7/Sep/2021
//	@brief 		_int_precision_csb
//	@return		size_t	-	the number of significant bits in vector<iptype>
//	@param		"a"	-	vector<iptype> operand
//
//	@todo
//
// Description:
//   Count number of significant bits in a vector<iptype> number
//
size_t _int_precision_csb(const std::vector<iptype> &a)
	{
	size_t bit_pos = 0, cnt;
	for (size_t i = a.size(); i > 0; --i)
		{
		if (a[i - 1] == 0) continue;
		cnt = _int_precision_clz(a[i - 1]);
		bit_pos = Bitsiptype - cnt;
		if (i - 1 > 0) 
			bit_pos += Bitsiptype * (i - 1);
		break;
		}
	return bit_pos;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		7/Aug/2021
//	@brief 		_int_precision_compare
//	@return 	int	-	The compare result. 0==equal, 1==s1>s2, -1==s1<s2
//	@param		"s1"	-	First operand to compare
//	@param		"s2"	-	Second operand to compare
//
//	@todo  
//
// Description:
//   Compare two unsigned vector<iptype> binary numbers 
//   and return 0 is equal, 1 if s1 > s2 otherwise -1
//   Optimized check length first and determine 1 or -1 if equal
//   compare the digits until a determination can be made.
//
int _int_precision_compare(const std::vector<iptype> *s1, const std::vector<iptype> *s2)
	{
	std::vector<iptype>::reverse_iterator s1_pos, s2_pos;
	int cmp; size_t i;

	if (s1->size() > s2->size())
		cmp = 1;
	else
		if (s1->size() < s2->size())
			cmp = -1;
		else
			{// Same size.
			s1_pos = const_cast<std::vector<iptype> *>(s1)->rbegin(); s2_pos = const_cast<std::vector<iptype> *>(s2)->rbegin();
			for (cmp = 0, i = s1->size(); i > 0; --i, ++s1_pos, ++s2_pos)
				{
				if (*s1_pos == *s2_pos) continue;
				if (*s1_pos > *s2_pos)
					{
					cmp = 1; break;
					}
				else
					{
					cmp = -1; break;
					}
				//if (i == 0) break;
				}
			}

	return cmp;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		7/Aug/2021, 31/Aug/2022
//	@brief 		_int_precision_compare
//	@return 	int	-	The compare result. 0==equal, 1==s1>s2, -1==s1<s2
//	@param		"s1"	-	First operand to compare
//	@param		"s2"	-	Second operand to compare
//
//	@todo  
//
// Description:
//   Compare two unsigned vector<iptype> binary numbers 
//   and return 0 is equal, 1 if s1 > s2 otherwise -1
//   Optimized check length first and determine 1 or -1 if equal
//   compare the digits until a determination can be made.
// 	modified from pointers to reference
//
int _int_precision_compare2( std::vector<iptype>& s1,  std::vector<iptype>& s2)
	{
	std::vector<iptype>::reverse_iterator s1_pos, s2_pos;
	int cmp; size_t i;

	if (s1.size() > s2.size())
		cmp = 1;
	else
		if (s1.size() < s2.size())
			cmp = -1;
		else
			{// Same size.
			s1_pos = s1.rbegin(); s2_pos = s2.rbegin();
			for (cmp = 0, i = s1.size(); i > 0; --i, ++s1_pos, ++s2_pos)
				{
				if (*s1_pos == *s2_pos) continue;
				if (*s1_pos > *s2_pos)
					{
					cmp = 1; break;
					}
				else
					{
					cmp = -1; break;
					}
				//if (i == 0) break;
				}
			}

	return cmp;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		20/Sep/2021
//	@brief 		std::vector<_TY> _precision_uadd64
//	@return 	std::vector<_TY> - 	the result of the addition. 2 dimensional vector
//	@param      "a"		-	_TY operand a
//	@param      "b"	   -	_TY operand b   
//
//	@todo
//
// Description:
//	 Generic Addition function for two vector<iptype> or vector<fptype> numbers
//  Add two unsigned iptype numbers togeher and return the result as a vector<iptype> [2] or vctor<fptype> [2]:
//	 where [0] is the lower operand and [1] is the upper operand
//
template<class _TY> std::vector<_TY> _precision_uadd64(const _TY a, const _TY b)
	{
	std::vector<_TY> res(2);
	res[0] = a + b;
	res[1] = res[0] < a ? 1 : 0;  // Carry
	return res;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		4/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_uadd_short
//	@return 	std::vector<iptype> - 	the result of the add
//	@param      "src1"	-	Source binary to add short number
//	@param      "d"	   -	Number to add.   
//
//	@todo
//
// Description:
//   Short Add: The digit d of iptype [0..2^64] for iptype=uint64_t is added to the unsigned binary vector 
//   Optimized 0 add or early out add is implemented
//
std::vector<iptype> _int_precision_uadd_short(const std::vector<iptype> *src1, const iptype d)
	{
	iptype carry;
	std::vector<iptype>::const_iterator s_pos;
	std::vector<iptype>::iterator d_pos;
	std::vector<iptype> des;

	if (d == 0)   // Zero add
		return *src1;

	carry = d;
	des = *const_cast<std::vector<iptype> *> (src1);		// Copy source to des1
	d_pos = des.begin();
	s_pos = src1->begin(); 

	for (; carry != 0 && d_pos != des.end(); ++s_pos, ++d_pos)
		{
		*d_pos += carry;
		if (*d_pos < *s_pos ) 
			carry = 1;  // Set Carry
		else carry = 0;
		}

	// Exhaust the smalles of the number, so only the carry can changes the uppper radix digits
	for (; carry != 0 && d_pos != des.end(); )
		{
		iptype tmp = *d_pos;
		*d_pos = tmp + carry;
		if (*d_pos < tmp) 
			carry = 1;  // Set Carry
		else carry = 0;
		++d_pos;
		}

	// No more carry or end of upper radix number. 
	if (carry != 0) // If carry add the carry as a extra digit to the front of the number
		des.push_back(1); 

	_int_precision_strip_trailing_zeros(&des);

	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		4/Aug/2021, 7-Aug-2022
//	@brief 		std::vector<iptype> _int_precision_uadd
//	@return		std::vector<iptype>	-	the result of adding src1 and src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   Add two unsigned decimal strings
//   Optimized: Used early out add
//	Now can handle a double carry correctly when adding the elements together
//
std::vector<iptype> _int_precision_uadd(const std::vector<iptype> *src1, const std::vector<iptype> *src2)
	{
	iptype carry = 0, tmp;
	std::vector<iptype> des;
	std::vector<iptype>::const_iterator pos, end;
	std::vector<iptype>::iterator d_pos;

	if (src1->size() >= src2->size())
		{
		des = *src1;
		pos = src2->begin();  
		end = src2->end(); 
		}
	else
		{
		des = *src2;
		pos = src1->begin(); 
		end = src1->end(); 
		}
	d_pos = des.begin();

	for (; pos != end; )
		{ // Adding element by element for the two numbers and correctly handle the carry along the way
		tmp = *pos + carry;
		carry = tmp < *pos ? 1 : 0;  // Next carry
		*d_pos += tmp;
		carry = *d_pos < tmp ? 1 : carry;
	//		*d_pos = *pos + *d_pos + carry;  // Wrong does not take into account a double carry
	//	carry = *d_pos < *pos ? 1 : 0;
		++pos;
		++d_pos;
		}

	// Exhaust the smallest of the number, so only the carry can changes the uppper radix digits
	for (; carry != 0 && d_pos != des.end();  )
		{
		tmp = *d_pos;
		*d_pos = tmp + carry;
		carry = *d_pos < tmp ? 1 : 0; 
		++d_pos;
		}

	// No more carry or end of upper radix number. 
	if (carry != 0) // If carry add the carry as a extra digit to the front of the number
		des.push_back(1);

	_int_precision_strip_trailing_zeros(&des);

	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		5/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_usub_short
//	@return 	std::vector<iptype> - 	the result of the add
//	@param      "src1"	-	Source string to add short number
//	@param      "d"	   -	iptype Number to add.   
// @param		"result" - Indicated wrap around (1) or not (0)
//
//	@todo
//
// Description:
//   Short subtract: The iptype digit d [0..2^64] is subtracted from the unsigned vector 
//   if src1 < d result is set to -1 (wrap around) otherwise result is set to  0 (no wrap around)
//   Optimized for 0 subtract
std::vector<iptype> _int_precision_usub_short(int *result, const std::vector<iptype> *src1, const iptype d)
	{
	iptype r, borrow=0;
	std::vector<iptype>::const_iterator pos;
	std::vector<iptype> des;

	if (d == 0) // Nothing to subtract
		{
		*result = 0;
		return *src1;
		}

	des.reserve(src1->capacity());  // Reserver space to avoid time consuming reallocation
	pos = src1->begin();
	r = *pos - (d + borrow);
	borrow = *pos < (d + borrow) ? 1 : 0;
	des.push_back(r);
	for (++pos; borrow>0 && pos != src1->end(); ++pos)
		{
		r = *pos - borrow;
		borrow = *pos < borrow ? 1 : 0;
		des.push_back(r);
		}
	_int_precision_strip_trailing_zeros(&des);

	*result = borrow > 0 ? -1 : 0;
	return des;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		4/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_usub
//	@return 	std::vector<iptype>	-	the result of subtracting src2 from src1
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
// @param		"result" - Return indicate wrap around (-1) otherwise 0
//
//	@todo
//
// Description:
//   Subtract two unsigned decimal strings
//   if src1 < src2 return -1 (wrap around) otherwise return 0 (no wrap around)
//
std::vector<iptype> _int_precision_usub(int *result, const std::vector<iptype> *src1, const std::vector<iptype> *src2)
	{
	iptype r, borrow=0;
	std::vector<iptype>::const_iterator pos1, pos2;
	std::vector<iptype> des1;

	if (src1->size() > src2->size())
		des1.reserve(src1->capacity());  // Reserver space to avoid time consuming reallocation
	else
		des1.reserve(src2->capacity());  // Reserver space to avoid time consuming reallocation
	pos1 = src1->begin();
	pos2 = src2->begin();

	for (; pos1 != src1->end() || pos2 != src2->end();)
		{
		if (pos1 != src1->end() && pos2 != src2->end())
			{
			r = *pos1 - (*pos2 + borrow);
			borrow = *pos1 < (*pos2 + borrow) ? 1 : 
					 *pos1==0? borrow : 0;      // if borrow was not paid then propagate it to next iptype subtraction
			++pos1; ++pos2;
			}
		else
			if ( pos1 != src1->end())
				{
				r = *pos1 - (borrow);
				borrow = *pos1 < borrow ? 1 : 0; 
				++pos1;
				}
			else
				{
				r = 0-(*pos2 + borrow);
				//borrow = 0 < (*pos2 + borrow) ? 1 : 0;
				borrow = 1;
				++pos2;
				}
		des1.push_back(r);
		}
	_int_precision_strip_trailing_zeros(&des1);

	*result = borrow>0 ? -1 : 0;
	return des1;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		5/Aug/2021
//	@brief 		std::vector<_TY> _precision_umul64
//	@return 	std::vector<_TY> - 	the result of the multiplication. 2 dimensional vector
//	@param      "a"		-	_TY operand a
//	@param      "b"	   -	_TY operand b   
//
//	@todo
//
// Description:
//	 Generic multiplication function for two vector<iptype> or vector<fptype> numbers
//  Multiply two unsigned iptype numbers togeher and return the result as a vector<iptype> [2]: where [0] is the lower operand and [1] is the upper operand
//
template<class _TY> inline std::vector<_TY> _precision_umul64(const _TY a, const _TY b)
	{
	const _TY mask = 0xffffffff;
	const unsigned int shift = 32;
	const _TY a0 = a & mask, a1 = a >> shift;
	const _TY b0 = b & mask, b1 = b >> shift;
	const _TY a0b0 = a0*b0, a0b1 = a0*b1, a1b0 = a1*b0, a1b1 = a1*b1;
	_TY carry0, carry1, mid;
	std::vector<_TY> res(2);
	mid = a0b1 + a1b0; carry1 = mid < a0b1 ? 1 : 0;
	res[0] = a0b0 + ( ( mid&mask ) << shift ); carry0 = res[0] < a0b0 ? 1 : 0;
	res[1] = (mid >> shift) + a1b1 + (carry1 << shift) + carry0; // no overflow
	return res;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		5/Aug/2021
//	@brief 		std::string _int_precision_umul_short
//	@return 	std::string - 	the result of the short multiplication
//	@param      "src1"	-	Source string to multiply short number
//	@param      "d"	   -	Number to multiply   
//
//	@todo
//
// Description:
//   Short Add: The digit d [0..RADIX] is multiplied to the unsigned decimal string
//   Optimized Multiply with zero yields zero, Multiply with one return the original 
//   
//
std::vector<iptype> _int_precision_umul_short(const std::vector<iptype> *src1, const iptype d)
	{
	iptype carry=0;
	std::vector<iptype>::const_iterator pos, end;
	std::vector<iptype> des;
	std::vector<iptype> tmp(2);

	if (d == 0)  // Multiply by zero is zero.
		{
		des.push_back(0);
		return des;
		}

	if (d == 1)  // Multiply by one dont change the src1.
		{
		des = *const_cast<std::vector<iptype> *> (src1);
		_int_precision_strip_trailing_zeros(&des);
		return des;
		}

	des.reserve(src1->capacity());  // Reserver space to avoid time consuming reallocation   
	pos = src1->begin();
	end = src1->end(); 

	for (; pos != end; ++pos)
		{
		tmp=_precision_umul64(*pos, d);
		tmp[0] += carry;
		carry = (tmp[0] < carry) ? tmp[1]+1 : tmp[1];
		des.push_back(tmp[0]);
		}

	if (carry != 0)
		des.push_back(carry);
	_int_precision_strip_trailing_zeros(&des);

	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		19/Jan/2022
//	@brief 		std::vector<iptype>  _int_precision_umul
//	@return		std::vector<iptype> -	the result of multiplying src1 and src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
// Multiply two unsigned vector<iptype> using the most optimal multiplication algorithm based on operand sizes
//
std::vector<iptype> _int_precision_umul(const std::vector<iptype> *src1, const std::vector<iptype> *src2)
	{
	std::vector<iptype> des;
	// Check for multiplication with 1 digit and use umul_short().
	if (src1->size() == 1)
		des = _int_precision_umul_short(src2, *src1->begin());
	else
		if (src2->size() == 1)
			des = _int_precision_umul_short(src1, *src2->begin());
		else
			if (src1->size() + src2->size() < 4000) // Use Schonhage-Strassen for multiplication
				des = _int_precision_umul_linear(src1, src2);
			else // Use FFT for multiplication
				des = _int_precision_umul_fourier(src1, src2);
	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/Aug/2021 
//	@brief 		std::vector<iptype>  _int_precision_umul_school
//	@return		std::vector<iptype> -	the result of multiplying src1 and src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
// Multiply two unsigned vector<iptype>.
// Not used anymore since the complexity is o(n^2)
//
std::vector<iptype> _int_precision_umul_school(const std::vector<iptype> *src1, const std::vector<iptype> *src2)
	{
	int disp;
	std::vector<iptype> des, tmp, offset;
	std::vector<iptype>::const_iterator pos2, pos2_end;

	pos2 = src2->begin();
	pos2_end = src2->end();

	des = _int_precision_umul_short(src1, *pos2);
	for (pos2++, disp = 1; pos2 != pos2_end; disp++, pos2++)
		{
		if (*pos2 != 0)
			{
			offset.push_back(0);
			tmp = _int_precision_umul_short(src1, *pos2);

			tmp.insert(tmp.begin(), disp, 0); // = offset + tmp;
			des = _int_precision_uadd(&des, &tmp);
			}
		}

	_int_precision_strip_trailing_zeros(&des);

	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		16/Aug/2021
//	@brief 		void convertbinary2double
//	@return		void	-	
//	@param		"dp"	-	pointer to array of doubles
//	@param		"d"		-	the binary iptype or fptype
//	@param	`   "bits"  -	Splitting bits (8 or 4)
//
//	@todo
//
// Description:
// convert an iptype or fptype into an array of doubles[]
//
template<class _TY> static size_t convertbinary2double(double *dp,  _TY d, const bool first, const int bits=8 )
	{
	const unsigned int mask = bits==8 ?0xff : 0xf;
	size_t k = 0;
	for (int i = sizeof(d) * 8 - bits; i >= 0; i -= bits)
		{
		_TY val = (d >> i) & mask;
		if (first == true && val == 0 && i != 0 && k==0) continue;
		k++;
		*dp++ = (double)(val);
		}
	return k;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		16/Aug/2021
//	@brief 		void convertbinary2double
//	@return		void	-	
//	@param		"dp"	-	pointer to array of doubles
//	@param		"d"		-	the binary iptype or fptype
//	@param	`   "bits"  -	Splitting bits (8 or 4)
//
//	@todo
//
// Description:
// convert an iptype or fptype into an array of doubles[]
// used by both _int_precision_umul_fourier and _float_precision_umul_fourier
//
/*
template<class _TY> static size_t convertbinary2double(std::vector<double>::iterator dp, _TY d, const bool first, const int bits = 8)
	{
	const unsigned int mask = bits == 8 ? 0xff : 0xf;
	size_t k = 0;
	for (int i = sizeof(d) * 8 - bits; i >= 0; i -= bits)
		{
		if (first == true && ((d >> i) & mask) == 0 && i != 0 && k == 0) continue;
		k++;
		*dp++ = (double)((d >> i) & mask);
		}
	return k;
	}
*/
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		16/Aug/2021
//	@brief 		void convertbinary2uint8
//	@return		void	-	
//	@param		"dp"	-	pointer to array of unsigned chars
//	@param		"d"		-	the binary iptype
//
//	@todo
//
// Description:
// convert an iptype into an array of unsigned bytes. Used in Schonhagen-Strassen
//
static inline size_t convertbinary2uint8(unsigned char *dp, iptype d, const bool first)
	{
	int k = 0;
	for (int i = Bitsiptype - 8; i >= 0; i -= 8)
		{
		if (first == true && ((d >> i) & 0xff) == 0 && i != 0 && k==0) continue;
		k++;
		*dp++ = (unsigned char)((d >> i) & 0xff);
		}
	return k;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		16/Oct/2021
//	@brief 		void convertbinary2uint16
//	@return		void	-	
//	@param		"dp"	-	pointer to array of unsigned shorts (16bits)
//	@param		"d"		-	the binary iptype
//
//	@todo
//
// Description:
// convert an iptype into an array of unsigned shorts (16bit). Used in Schonhagen-Strassen
//
static inline size_t convertbinary2uint16(unsigned short *dp, iptype d, const bool first)
	{
	int k = 0;
	for (int i = Bitsiptype - 16; i >= 0; i -= 16)
		{
		if (first == true && ((d >> i) & 0xffff) == 0 && i != 0 && k == 0) continue;
		k++;
		*dp++ = (unsigned short)((d >> i) & 0xffff);
		}
	return k;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		17/Oct/2021
//	@brief 		void convertbinary2Halfiptype
//	@return		void	-	
//	@param		"dp"	-	pointer to array of unsigned shorts (16bits)
//	@param		"d"		-	the binary iptype
//
//	@todo
//
// Description:
// convert an iptype into an array of unsigned shorts (16bit). Used in Schonhagen-Strassen
//
static inline size_t convertbinary2Halfiptype(iptype *dp, iptype d, const bool first)
	{
	const unsigned int HalfBitsiptype = Bitsiptype / 2;
	const iptype mask = (~(iptype)0) >> HalfBitsiptype;
	int k = 0;
	if (first == false || ((d >> HalfBitsiptype) & mask) != 0)
		{*dp++ = d >> HalfBitsiptype; ++k; }
	*dp++ = d & mask; k++;
	return k;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		16/Aug/2021
//	@brief 		void convertdoublebinary
//	@return		iptype	-	return the binary constructed number
//	@param		"dp"	-	pointer to array of doubles
//	@param		"bits"	-	Number of bits in *dp (either 8 or 4 bits)
//
//	@todo
//
// Description:
// convert an an array of doubles[] into a binary iptype
// used by both _int_precision_umul_fourier and _float_precision_umul_fourier
//
static iptype convertdouble2binary(double *dp, size_t maxinx, const double cy=0, const int bits=8)
	{
	const unsigned int mask = bits==8? 0xff: 0xf;
	iptype d=(unsigned char)cy;
	for (int i = 0; i < maxinx; ++i)
		{
		d <<= bits;
		d |= (unsigned char)(*dp++)&mask;
		if (bits == 4 && ++i < maxinx)
			{
			d <<= bits;
			d |= (unsigned char)(*dp++) & mask;
			}
		}
	return d;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		17/Aug/2021
//	@brief 		void convertHalfiptype2binary
//	@return		iptype	-	return the binary constructed number
//	@param		"dp"	-	pointer to array of half iptypes 
//
//	@todo
//
// Description:
// convert an an array of half iptypes [] into a binary iptype. Used in Schonhagen-Strassen
//
static inline iptype convertHalfiptype2binary( uintmax_t *dp, const size_t maxinx, const unsigned int cy = 0)
	{
	iptype d = (iptype)cy;
	if( maxinx > 1 )
		d |= dp[ 1 ];
	d <<= ( Bitsiptype / 2);
	d |= dp[0];
	return d;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		10/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_umul_fourier
//	@return 	std::vector<iptyp> -	the result of multplying src1 and src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//		Multiply two unsigned binary numbers
//		Optimized: Used FFT algorithm to performed the multiplication
//		Since we convert the binary numbers into float we have to ensure proper accuracy in the calculation.
//		In numerical recipies in C (2nd edition, chaper 20.6) they state that using double the equations that need to be fulfilled for accuracy
//		1byte binary:
//			log2(256^2)+log2(N)+safeguard*log2(log2(N))< digits in double which is 53. using safeguard as 1 we get for N=10^9 digits
//			16+30+4.9=50.9  which should be just Ok for 1 byte binary digits.
//		2byte binary:
//			log2(256^2^2)+log2(N)+safeguard*log2(log2(N))< digits in double which is 53. using safeguard as 1 we get for N=10^5 digits
//			32+16.6+4.05=52.65  10^5 digits is not enough for arbitrary precsion so we are only using 1byte.	
//
std::vector<iptype> _int_precision_umul_fourier(const std::vector<iptype> *src1, const std::vector<iptype> *src2, int nbits)
	{
	std::vector<iptype> des;
	size_t n, l, l1, l2, j;
	double cy;
	int bits = nbits == 0 ? 8 : nbits;
	int radix = bits == 8 ? 256 : 16;
	size_t sz = sizeof(iptype);
	std::vector<double> va, vb;

	l1 = src1->size();
	l2 = src2->size();
	des.reserve(l1 + l2 + 16);  // Ensure enough space to hold the Multiplication result to avoid reallocation of des
	l = l1 < l2 ? l2 : l1;  
	// Since we split the 64bit numbers into chunk of 8bit to ensure we have enough accuray when using double 
	l *= sizeof(iptype);  // Convert to byte 
	if (l > 6'000'000*sizeof(iptype) || bits == 4)
		{
		bits = 4; l <<= 1; radix = 16; sz *= 2;  // use 2^4 instead of 2^8
	//	std::cout << "int_umul_fourier" << " do 4bit" << std::endl;  // DEBUG
		}
	for (n = 1; n < l; n <<= 1) ;
	n <<= 1;
	
#ifdef HVE_THREAD
	// Using parallel sections below speeds up the performance of the two calls to _int_real_Fourier() with a factor of 1.8 
	if (nbits == 0 || l1 + l2>10'000)
		{// Starting thread using lambda expressions
		// L1, l2, va, vb by reference since it is used after the thread has terminated
		std::thread first([&, n, bits]()
			{std::vector<iptype>::const_reverse_iterator pos, end;
			size_t i;
			va.resize(n);
			for (i = 0, pos = src1->rbegin(), end = src1->rend(); pos != end; ++pos)
				i += convertbinary2double(&va[i], *pos, i == 0, bits );
			l1 = i; // L1 now Number of bytes or nibbles
			_vector_real_fourier(va, n, 1); // FFT va
			});

		std::thread second([&, n, bits]()
			{std::vector<iptype>::const_reverse_iterator pos, end;
			size_t i;
			vb.resize(n);
			for (i = 0, pos = src2->rbegin(), end = src2->rend(); pos != end; ++pos)
				i += convertbinary2double(&vb[i], *pos, i == 0, bits);
			l2 = i; // L2 now Number of bytes or nibbles
			_vector_real_fourier(vb, n, 1); // FFT vb
			});

		first.join();
		second.join();
	}
	else
#endif
		{
		std::vector<iptype>::const_reverse_iterator pos, end;
		va.resize(n);
		vb.resize(n);
		// Start with most significant fptype e.g. src1[0]
		for (l1 = 0, pos = src1->rbegin(), end = src1->rend(); pos != end; ++pos)
			l1 += convertbinary2double(&va[l1], *pos, l1 == 0, bits);
		// L1 now Number of bytes or nibbles
		// Start with most significant fptype e.g. src2[0]
		for (l2 = 0, pos = src2->rbegin(), end = src2->rend(); pos != end; ++pos)
			l2 += convertbinary2double(&vb[l2], *pos, l2 == 0, bits);
		// L2 now number of bytes or nibbles
		_vector_real_fourier(va, n, 1); // FFT va
		_vector_real_fourier(vb, n, 1); // FFT vb
		}

	vb[0] *= va[0];
	vb[1] *= va[1];
	for (j = 2; j < n; j += 2)
		{
		double t;
		vb[j] = (t = vb[j])*va[j] - vb[j + 1] * va[j + 1];
		vb[j + 1] = t*va[j + 1] + vb[j + 1] * va[j];
		}
	_vector_real_fourier(vb, n, -1);
	for (cy = 0, j = 0; j <= n - 1; ++j)
		{
		double t;
		t = vb[n - 1 - j] / (n >> 1) + cy + 0.5;
		cy = (unsigned long)(t / radix);  // Byte Radix 2^8 or 2^4
		vb[n - 1 - j] = t - cy * radix;
		}

	// Now collect then back into a vector<fptype> format
	l1 += l2 - 1;				// max number of bytes or nibbles plus the carry
	l2 = l1 / sz;   // Number of full 64bit digits
	for (l = l2; l > 0; --l)	// do the full 64bit integers first starting backwards from b
		{
		iptype num;
		size_t inx = l1 - sz*(l2 - l + 1);
		num = convertdouble2binary(&vb[inx], sz, 0, bits);
		des.push_back(num);
		}
	l2 = l1 % sz;   // Number of remaing 8bits or 4bits digits
	if (l2>0 || cy != 0)		// do the the last 64bit integers from b
		{
		iptype num;
		num = convertdouble2binary(&vb[0], l2, cy, bits);
		des.push_back(num);
		}
	_int_precision_strip_trailing_zeros(&des);

	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		17-Aug-2021
//	@brief 		std::vector<iptype> _int_precision_umul_karatsuba_
//	@return 	std::vector<iptype>	-	the result of multiplying src1 and src2
//	@param		"lhs"	-	First unsigned source argument
//	@param		"rhs"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   Multiply two unsigned binary numbers of vector<iptype>, using the karatsuba method
//	 Karatsuba is faster than umul_fourier with operands up to	 xx decimal digits whereafter _umul_fourier is faster
//   Notice when operands can fit into a 64bit integer we switch to native multiplications.
//
std::vector<iptype> _int_precision_umul_karatsuba(const std::vector<iptype> *lhs, const std::vector<iptype> *rhs)
	{
	std::vector<iptype> result, z0, z1, z2, z3;
	std::vector<iptype> lhshigh, lhslow, rhshigh, rhslow;
	int wrap;
	size_t half_length, length, l_length = lhs->size(), r_length = rhs->size(), tot_len;
	length = l_length < r_length ? r_length : l_length;
	tot_len = l_length + r_length;

	// Short cuts
	if ((l_length == 1 && lhs->front() == 0) || (r_length==1 && rhs->front() == 0) )
		{//lhs * 0 or rhs*0 is zero
		result=std::vector<iptype> (1,0);
		return result;
		}
	if (l_length == 1 && lhs->front() == 1 )
		{
		//rhs * 1 is rhs
		result = *rhs;
		return result;
		}
	if( r_length == 1 && rhs->front() == 1)
		{
		//lhs * 1 is lhs
		result = *lhs;
		return result;
		}
	if (l_length==1 && r_length==1)  // If max digits in lhs & rhs less than to fit into a 64 bit integer then do it the binary way
		{
		result = _precision_umul64(lhs->front(), rhs->front());
		if (result[0] == 0 || result[1] == 0)
			wrap = 0;
		_int_precision_strip_trailing_zeros(&result);
		return result;
		}

	// Splitting
	half_length = length >> 1;
	if (l_length <= half_length)
		{
		lhshigh.insert(lhshigh.begin(), 1, 0); lhslow.insert(lhslow.begin(),lhs->begin(), lhs->end() );
		}
	else if (l_length < length)
		{
		if (half_length >= l_length )
			lhshigh.insert(lhshigh.begin(), 1, 0);
		else
			lhshigh.insert(lhshigh.begin(),lhs->begin()+half_length,lhs->end());
		lhslow.insert(lhslow.begin(), lhs->begin(), lhs->begin() + half_length );
		}
	else
		{
		lhslow.insert(lhslow.begin(), lhs->begin(), lhs->begin() + half_length);
		lhshigh.insert(lhshigh.begin(), lhs->begin() + half_length, lhs->end());
		}

	if (r_length <= half_length)
		{
		rhshigh.insert(rhshigh.begin(), 1, 0); 
		rhslow.insert(rhslow.begin(),rhs->begin(), rhs->end() );
		}
	else if (r_length < length)
		{
		if (half_length >= r_length ) 
			rhshigh.insert(rhshigh.begin(), 1, 0);
		else 
			rhshigh.insert(rhshigh.begin(), rhs->begin()+half_length, rhs->end());
		rhslow.insert(rhslow.begin(), rhs->begin(), rhs->begin()+half_length);
		}
	else
		{
		rhslow.insert(rhslow.begin(), rhs->begin(), rhs->begin() + half_length);
		rhshigh.insert(rhshigh.begin(), rhs->begin() + half_length, rhs->end());
		}

	// Evaluation
	z0 = _int_precision_uadd(&lhshigh, &lhslow);
	z1 = _int_precision_uadd(&rhshigh, &rhslow);
	z2 = _int_precision_umul_karatsuba(&z0, &z1);
	z1 = _int_precision_umul_karatsuba(&lhslow, &rhslow);
	z0 = _int_precision_umul_karatsuba(&lhshigh, &rhshigh);
	z3 = _int_precision_uadd(&z0, &z1);
	z3 = _int_precision_usub(&wrap, &z2, &z3);

	// Recomposition
	z0 = _int_precision_ushiftleft( &z0, Bitsiptype*2*half_length);
	z3 = _int_precision_ushiftleft( &z3, Bitsiptype*half_length);
	z0 = _int_precision_uadd(&z0, &z1);
	result = _int_precision_uadd(&z0, &z3);
	return result;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		16/Oct/2021
//	@brief 		std::vector<iptype> _int_precision_schonhage_strassen_linear_umul
//	@return		std::vector<iptype>	-	the result of multiplying src1 and src2
//	@param		"lhs"	-	First unsigned source argument
//	@param		"rhs"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   Multiply two unsigned decimal strings, using the Schonhage-Strassen (linear convolution) method
//   Using the full advantages of Half bit size of iptype in the linear convolution
//
std::vector<iptype> _int_precision_umul_linear(const std::vector<iptype> *lhs, const std::vector<iptype> *rhs)
	{
	size_t i, j;
	size_t l_length = lhs->size(), r_length = rhs->size();
	const unsigned int HalfBitsiptype = Bitsiptype >> 1;  // same as / 2
	const iptype mask = (~(iptype)0) >> HalfBitsiptype;
	const uintmax_t radix = (uintmax_t)1 << HalfBitsiptype;
	std::vector<iptype> res;
	std::vector<iptype>::const_reverse_iterator pos, end;
	std::vector<uintmax_t> linearconvolution( 2 * (l_length + r_length), 0);  // initialize it with zero
	std::vector<iptype> ua( l_length * 2), ub( r_length * 2 );
	uintmax_t nextCarry = 0;

	// Convert to half iptype from vector<iptype> and notice we stored in reverse order 
	// by first converting lhs onto ua and then rhs into ub
	// e.g. lhs=a0+a1*R+a2*R^2,...an-1*R^n-1, a0 can be subdivied into half iptype  from iptype by mapping each mNumber number into 2 half iptype numbers 
	// the function convertbinary2Halfiptype() does this job per mNumber iptype number
	for (i = 0, pos = lhs->rbegin(), end=lhs->rend(); pos != end; ++pos)
		i += convertbinary2Halfiptype(&ua[i], *pos, i == 0);
	l_length = i;  // l_length now in half iptype's  instead of iptype's
	for (j = 0, pos = rhs->rbegin(), end=rhs->rend(); pos != end; ++pos)
		j += convertbinary2Halfiptype(&ub[j], *pos, j == 0);
	r_length = j;  // l_length now in half iptype's instead of iptype's

	// do the linear convolution
	for (i = 0; i < r_length; ++i)
		for (j = 0; j < l_length; ++j)
			{
			uintmax_t m = (uintmax_t)ub[r_length - 1 - i] * (uintmax_t)ua[l_length - 1 - j];
			linearconvolution[i + j] += m;
			if (linearconvolution[i + j] < m) // carry
				{
				// Propagate carry
				for (size_t k = i + j + 1; ; ++k)
					{
					linearconvolution[k] += radix;	// Add carry
					if (linearconvolution[k] >= radix)	// Continue adding carry?
						{
#ifdef DEBUG_HVE
						if(k-(i+j) > 1)
							std::cout << "Propagate end in linear convolution with " << k-(i+j) << " Steps" << std::endl; // DEBUG HVE
#endif
						break;
						}
					}
				}
			}

	res.reserve( r_length + l_length + 2 );
	for (i = 0; i < l_length + r_length - 1; ++i)
		{
#ifdef DEBUG_HVE
		if (linearconvolution[i] +nextCarry < nextCarry)  // Carry
			std::cout << "Overflow in linear convolution before" << " Nextcarry=" << nextCarry << " LC[i]=" << linearconvolution[i] << " i=" << i << std::endl; // DEBUG HVE
#endif
		linearconvolution[i] += nextCarry;
		if (linearconvolution[i] < nextCarry)  // Carry
			{
			size_t k;
#ifdef DEBUG_HVE
			std::cout << "Overflow in linear convolution after" << " Nextcarry=" << nextCarry << " LC[i]=" << linearconvolution[i] << std::endl; // DEBUG HVE
#endif
			// Propagate carry
			for (k = i + 1; ; ++k)
				{
				linearconvolution[k] += radix;	// Add carry
				if (linearconvolution[k] >= radix)	// Continue adding carry?
					break;
				}
#ifdef DEBUG_HVE
			std::cout << "Overflow Propagate Overflow end in linear convolution with " << k - (i) << " Steps" << std::endl; // DEBUG HVE
#endif
			}
		nextCarry = linearconvolution[i] >> HalfBitsiptype;  //  same as / radix;
		linearconvolution[i] &= mask;  // same as %= radix;
		}
	if (nextCarry != 0)
		{
	//	std::cout << "End Carry=" << nextCarry << std::endl; // DEBUG HVE
		linearconvolution[i++] = nextCarry & mask; // same as % radix;
		}

	//linearconvolution now holds the result with [0] as the most significant byte number and [i-1] as the least sinificant byte as Halfbitsiptype numbers
	// Now convert then back into a vector<iptype> format. i is the number of HalfBitsiptype's in the result
	// do the full 64bit integers first starting from least significant HalfBitsiptype
	for ( j = 0; j < i; j+=2 ) 
		{
		iptype num;
		num = convertHalfiptype2binary(&linearconvolution[j], 2);
		res.push_back(num);
		}

	_int_precision_strip_trailing_zeros(&res);
	return res;
	}

// Short Division: The ptype digit d  is divide up into the unsigned ptype vector
//
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/Aug/2021
//	@brief 		std::vector<iptype>	_int_precision_udiv_short
//	@return 	std::vector<iptype>	- The result of the short division
//	@param      "src1"				- Source string to divide with the short number
//	@param      "d"					- Number to divide
//	@param		"remaind"			- The remaind of the short division
//
//	@todo
//
// Description:
//   Short divide: The ptype digit d [0..2^32] is divided up in the unsigned vector<iptype> 
//	  Notice only up to int 32bit can be handle as short div.
//   Divide with zero throw an exception
//
std::vector<iptype> _int_precision_udiv_short(iptype *remaind, const std::vector<iptype> *src1, const iptype d)
	{
	const unsigned int shifts = 4 * sizeof(iptype);
	const iptype mask = (~((iptype)0) ) >> shifts;
	iptype ir;
	std::vector<iptype>::const_reverse_iterator s_pos, s_end;
	std::vector<iptype> des;

	if (d == 0)
		{
		throw int_precision::divide_by_zero();
		}

	if (d == 1)  // Divide by one dont change the src1.
		{
		des = *const_cast<std::vector<iptype> *>(src1);
		_int_precision_strip_trailing_zeros(&des);
		*remaind = 0;
		return des;
		}

	des.reserve(src1->capacity());  // Reserver space to avoid time consuming reallocation
	s_pos = src1->rbegin();
	s_end = src1->rend();
	for ( ir=0; s_pos != s_end; ++s_pos)
		{
		iptype n, qh, ql;
		/*if (ir == 0)
			{// Shortcut when carry is zero
			des.push_back(*s_pos / d );
			ir = *s_pos % d;
			}
		else*/
			{
			n = *s_pos >> shifts;
			n |= ir << shifts;
			qh = n / d;	ir = n % d;
			n = *s_pos & mask;
			n |= ir << shifts; 
			ql = n / d;	ir = n % d;
			n = (qh << shifts) | ql;
			des.push_back(n);
			}
		}

	reverse(des.begin(), des.end());
	_int_precision_strip_trailing_zeros(&des);
	*remaind = ir;
	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_udiv
//	@return		std::vector<iptype>-	the result of disivison
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   Divide two unsigned binary numbers
//   Optimized: Used early out add and multiplication w. zero
//
std::vector<iptype> _int_precision_udiv(const std::vector<iptype> *src1, const std::vector<iptype> *src2)
	{
	int plusdigit, plusbits, wrap, i;
	iptype underflow;
	std::vector<iptype> des, tmp, quotient, divisor;

	des.push_back(0);
	divisor = *const_cast<std::vector<iptype> *> (src1);
	if (src2->size() == 1 && ( src2->front() >> 32 ) == 0 ) // Make short div if denominator <= 32 bit integer.
		return _int_precision_udiv_short( &underflow, &divisor, src2->front());

	plusdigit = (int)divisor.size() - (int)src2->size();
	if (plusdigit < 0)  //src1 / src2 == 0
		return des;

	plusbits = (int)_int_precision_clz(src2->back())- (int)_int_precision_clz(divisor.back()) ;
	plusbits=plusdigit * Bitsiptype + plusbits;
	for(i=0; plusbits >= 1 ; ++i) 
		{
		tmp = _int_precision_ushiftleft(src2, plusbits);
		if (_int_precision_compare(&divisor, &tmp) < 0)
			{ // Too much reduce with one power of radix
			--plusbits; continue;
			}
		divisor = _int_precision_usub(&wrap, &divisor, &tmp);
		quotient.clear();
		quotient.insert(quotient.begin(), (plusbits / Bitsiptype ) + 1, 0);
		quotient[quotient.size() - 1] = (iptype)(1) << (plusbits % Bitsiptype );
		des = _int_precision_uadd(&des, &quotient);
		}

	for (wrap = 0; wrap == 0; )
		{
		divisor = _int_precision_usub(&wrap, &divisor, src2);
		if (wrap == 0) // src1 was indeed > src2
			des = _int_precision_uadd_short(&des, 1);
		}

	_int_precision_strip_trailing_zeros(&des);
	return des;
	}

// Short Remainder: The iptype digit d [1..2^64] is divide up into the unsigned vector<iptype> and the remaing is returned
//
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6-Aug-2021
//	@brief 		std::vector<iptype> _int_precision_urem_short
//	@return 	std::vector<iptype> - 	the result of the short remainder
//	@param      "src1"	-	Source string to divide with the short number
//	@param      "d"	   -	Number to divide
//
//	@todo
//
// Description:
//   Short remainder: The iptype digit d [0..2^64] is divided up in the unsigned vector<iptype> and the remaing is retuened
//   Divide with zero throw an exception
//   if d==1 then result == 0, for d==2,4,5,8,10 we only test the last few digits to get the result. This speed up rem for large integers with small rem value
//   since we dont have to run through every digits in src1
//
std::vector<iptype> _int_precision_urem_short(const std::vector<iptype> *src1, const iptype d)
	{
	const unsigned int shifts = 4 * sizeof(iptype);
	const iptype mask = (~((iptype)0)) >> shifts;
	iptype ir;
	std::vector<iptype>::const_reverse_iterator s_pos, s_end;
	std::vector<iptype> des;

	if (d == 0)
		{
		throw int_precision::divide_by_zero();
		}

	if (d == 1)  // Remainer is always 0 for d==1
		{
		des.push_back(0);
		return des;
		}

	// Short cut
	ir = *const_cast<std::vector<iptype> *>(src1)->begin();
	switch (d)
		{
		case 2: des.push_back(ir % 2); return des; break;
		case 4: des.push_back(ir % 4); return des; break;
	//	case 5: des.push_back(ir % 5); return des; break;
		case 8: des.push_back(ir % 8); return des; break;
	//	case 10: des.push_back(ir % 10); return des; break;
		case 16: des.push_back(ir % 16); return des; break;
		default:;   // No Short cut
		}
	

	s_pos = src1->rbegin();
	s_end = src1->rend();
	for (ir = 0; s_pos != s_end; ++s_pos)
		{
		iptype n;
		if (ir == 0)
			{// Shortcut when carry is zero
			ir = *s_pos % d;
			}
		else
			{
			n = *s_pos >> shifts;
			n += ir << shifts;
			ir = n % d;
			n = *s_pos & mask;
			n += ir << shifts;
			ir = n % d;
			}
		}

	des.push_back(ir);
	return des;
	 }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		7/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_urem
//	@return		std::vector<iptype>	-	the remaing result of divide src1 with src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   Find the remainder when divide two unsigned vector<iptype> numbers
//   Optimized: Used early out add and multiplication w. zero
//

std::vector<iptype> _int_precision_urem(const std::vector<iptype> *src1, const std::vector<iptype> *src2)
	{
	int wrap;
	std::vector<iptype> des, tmp;

	des.push_back(0);
	if (src2->size() == 1 &&  (src2->front() >> 32) == 0 ) // Make short rem 
		{
		iptype rem;
		//des = _int_precision_urem_short( src1, src2->front());
		//return _int_precision_udiv_short(&underflow, &divisor, src2->front());
		_int_precision_udiv_short(&rem, src1, src2->front());
		des[0] = rem;
		return des;
		}

	tmp = _int_precision_udiv( src1, src2);
	tmp = _int_precision_umul( &tmp, src2);
	des = _int_precision_usub(&wrap, src1, &tmp);

	_int_precision_strip_trailing_zeros(&des);

	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		23/Oct/2021
//	@brief 		std::vector<iptype> _int_precision_udivrem
//	@return		std::vector<iptype>	-	the quotient.  Result of divide src1 with src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
// @param		"r"		-	Pointer to remainder of the division
//
//	@todo
//
// Description:
//   Find both the remainder and the qoutient when divide two unsigned vector<iptype> numbers
//   
//
std::vector<iptype> _int_precision_udivrem(std::vector<iptype> *src1, std::vector<iptype> *src2, std::vector<iptype> *r )
	{
	int wrap;
	std::vector<iptype> des, tmp;

	des.push_back(0);
	if (src2->size() == 1 && (src2->front() >> 32) == 0) // Make short rem 
		{
		iptype rem;
		//des = _int_precision_urem_short( src1, src2->front());	
		//return _int_precision_udiv_short(&underflow, &divisor, src2->front());
		des = _int_precision_udiv_short(&rem, src1, src2->front());
		tmp.assign(1, rem);
		*r = tmp;
		return des;
		}

	des = _int_precision_udiv(src1, src2);
	tmp = _int_precision_umul(&des, src2);
	tmp = _int_precision_usub(&wrap, src1, &tmp);

	_int_precision_strip_trailing_zeros(&des);
	*r = tmp;
	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_ushiftright
//	@return 	std::vector<iptype> - The negated number	
//	@param      "src"	-	The number to be shifted
// @param		"shift"	-	The shift count
//
//	@todo
//
// Description:
// Implement the shift left operation src >> shift
//
std::vector<iptype> _int_precision_ushiftright(const std::vector<iptype> *src, const size_t shift )
	{
	size_t shiftwidth = Bitsiptype;
	std::vector<iptype>::const_reverse_iterator pos, end;
	std::vector<iptype> des; 
	size_t discard, within;
	size_t i;
	iptype carry, mask;
	
	// Determine how many full digit shift and the last shift (last shift = shift count % Bitsiptype
	if (src->size()==1 && src[0]==0)  // Short cut: a zero number zero shifting right is still zero.
		return *src;
	if (shift == 0)  // Short cut: shift zero right does not change the number.
		return *src;

	discard = shift / Bitsiptype;
	within = shift % Bitsiptype;
	shiftwidth -= within;
	mask = (~((iptype)0)) >> shiftwidth;  // check is resukltis still ok for shiftwidth == 64
	i = src->size();
	for (carry=0, pos = src->rbegin(), end=src->rend(); pos != end && i>discard; --i, ++pos)
		{
		iptype n, nextcarry;
		n = *pos;
		nextcarry = n & mask;
		n >>= within;
		carry <<= shiftwidth;  // check shiftwidth==64 abd the resukt is stilkl ok
		n |= carry;
		carry = nextcarry;
		des.push_back(n);
		}
	reverse(des.begin(), des.end());// check for des.size()==0
	if (des.size() == 0)
		des.push_back(0);
	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_ushiftleft
//	@return 	std::vector<iptype> - The negated number	
//	@param      "src"	-	The number to be shifted
// @param		"shift"	-	The shift count
//
//	@todo
//
// Description:
// Implement the shift left operation src1 << shift
//
std::vector<iptype> _int_precision_ushiftleft(const std::vector<iptype> *src, const size_t shift)
	{
	const unsigned int shiftwidth = 8 * sizeof(iptype);
	std::vector<iptype>::const_iterator pos, end;
	std::vector<iptype> des;
	size_t adding, within;
	iptype carry, mask;

	// Determine how many full digit shift and the last shift (last shift = shift count % sizeof(iptype)*8
	if (src->size() == 1 && src[0] == 0)  // Short cut: a zero number zero shifting left is still zero.
		return *src; 
	if (shift == 0)  // Short cut: shift zero left does not change the number.
		return *src; 

	adding = shift / Bitsiptype;
	within = shift % Bitsiptype;
	//for (; adding > 0; --adding)
	//	des.push_back(0);
	if( adding > 0 )
		des.insert(des.begin(), adding, 0);
	mask = (~((iptype)0)) >> within;  mask = ~mask;
	for (carry = 0, pos = src->begin(), end=src->end(); pos != end; ++pos)
		{
		iptype n, nextcarry;
		n = *pos;
		nextcarry = n & mask;
		n <<= within;
		n |= carry >> (shiftwidth - within);    // check shiftwidth 64 and witin == 0
		carry = nextcarry;
		des.push_back(n);
		}
	if (carry != 0)
		des.push_back(carry >> (shiftwidth - within)); // check shiftwidth 64 and witin == 0

	return des;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		7/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_uand
//	@return		std::vector<iptype>	-	the result of anding src1 and src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   And two unsigned vector<iptype> numbers
//   Optimized: Used early out and
//
std::vector<iptype> _int_precision_uand(const std::vector<iptype>& src1, const std::vector<iptype>& src2)
	{
	std::vector<iptype> des;
	std::vector<iptype>::const_iterator pos;
	std::vector<iptype>::iterator d_end, d_pos;

	// Making the shortest operand the result operand since that will be the maximum number of digits
	if (src1.size() >= src2.size()) 
		{
		des = src2; // *const_cast<std::vector<iptype> *>(src2);
		pos = src1.begin();
		d_end = des.end();
		}
	else
		{
		des = src1; // *const_cast<std::vector<iptype> *>(src1);
		pos = src2.begin();
		d_end = des.end();
		}
	d_pos = des.begin();

	for (; d_pos != d_end; ++pos, ++d_pos )
		{ // Anding element by element for the two numbers
		*d_pos &= *pos;
		}

	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		7/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_uor
//	@return		std::vector<iptype>	-	the result of oring src1 and src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   Or two unsigned vector<iptype> numbers
//   Optimized: Used early out and
//
std::vector<iptype> _int_precision_uor(const std::vector<iptype>& src1, const std::vector<iptype>& src2)
	{
	std::vector<iptype> des;
	std::vector<iptype>::const_iterator pos, end;
	std::vector<iptype>::iterator d_pos;

	// Making the shortest operand the result operand since that will be the maximum number of digits
	if (src1.size() >= src2.size()) 
		{
		des = src1; //  *const_cast<std::vector<iptype> *>(src1);
		pos = src2.begin();
		end = src2.end();
		}
	else
		{
		des = src2; // *const_cast<std::vector<iptype> *>(src2);
		pos = src1.begin();
		end = src1.end();
		}
	d_pos = des.begin();

	for (; pos != end; ++pos, ++d_pos)
		{ // oring element by element for the two numbers
		*d_pos |= *pos;
		}

	return des;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		7/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_uxor
//	@return		std::vector<iptype>	-	the result of xoring src1 and src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   Xor two unsigned vector<iptype> numbers
//   Optimized: Used early out and
//
std::vector<iptype> _int_precision_uxor(const std::vector<iptype>& src1, const std::vector<iptype>& src2)
	{
	std::vector<iptype> des;
	std::vector<iptype>::const_iterator pos, end;
	std::vector<iptype>::iterator d_pos;

	// Making the shortest operand the result operand since that will be the maximum number of digits
	if (src1.size() >= src2.size()) 
		{
		des = src1; // *const_cast<std::vector<iptype> *>(src1);
		pos = src2.begin();
		end = src2.end();
		}
	else
		{
		des = src2; //  *const_cast<std::vector<iptype> *>(src2);
		pos = src1.begin();
		end = src1.end();
		}
	d_pos = des.begin();

	for (; pos != end; ++pos, ++d_pos)
		{ // xoring element by element for the two numbers
			*d_pos ^= *pos;
		}

	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/Sep/2021
//	@brief 		std::vector<iptype> _int_precision_unegate
//	@return		std::vector<iptype>	-	the result of negating src1 
//	@param		"src1"	-	First unsigned source argument
//
//	@todo
//
// Description:
//   Negate unsigned vector<iptype> number
//
std::vector<iptype> _int_precision_unegate(const std::vector<iptype>& src1)
	{
	std::vector<iptype> des = src1;
	std::vector<iptype>::iterator d_pos;

	for (d_pos = des.begin(); d_pos != des.end(); ++d_pos)
		{ // negating element of the number
		*d_pos = ~*d_pos;
		}

	return des;
	}

///////////////////////////////////////////////
//
//
//    To and from string conversion including power of 10 tables for up to 64bit unsigned integers
//		or higher power of tables for splitting numbers into trunks, kilotrunks, megatrunks or
//		gigatrunks.
//
//
///////////////////////////////////////////////

// Same table as above but converted to float_precisions. Need to be sure that defaul precision is >=20
static float_precision _fpPowerof10Table[20] = {float_precision(_powerof10Table[0]),float_precision(_powerof10Table[1]),float_precision(_powerof10Table[2]),
												float_precision(_powerof10Table[3]),float_precision(_powerof10Table[4]),float_precision(_powerof10Table[5]),
												float_precision(_powerof10Table[6]),float_precision(_powerof10Table[7]),float_precision(_powerof10Table[8]),
												float_precision(_powerof10Table[9]),float_precision(_powerof10Table[10]),float_precision(_powerof10Table[11]),
												float_precision(_powerof10Table[12]),float_precision(_powerof10Table[13]),float_precision(_powerof10Table[14]), 
												float_precision(_powerof10Table[15]),float_precision(_powerof10Table[16]),float_precision(_powerof10Table[17]),
												float_precision(_powerof10Table[18]),float_precision(_powerof10Table[19]) };


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		25/Oct/2021
//	@brief 		std::string check_digits
//	@return 	bool	-	Return true if digits is valid other throw an exception and return false
//	@param		"s"	-		String of digits to check
//	@param		"start"	-	Value to convert to ascii string based on RADIX
//	@param		"end"	-	RADIX value of conversion 
//	@param		"base"	-	Check for base (default BASE_10)
//
//	@todo  Add to do things	
//
// Description:
//		Check the ascii string for valid digits according to base.
//
static bool check_digits(const std::string s, const size_t start, const size_t end, const int base = BASE_10)
	{
	if (base <= BASE_10)
		{
		for (size_t i = start; i < end; i++)
			{
			if (s[i] < '0' || s[i] >= base + '0')
				{
				throw int_precision::bad_int_syntax();
				return false;
				}
			}
		}
	else
		for (size_t i = start; i < end; i++)
			{ // base 11..36
			if ((s[i] < '0' || s[i] > '9') && (tolower(s[i]) < 'a' || tolower(s[i]) > base -BASE_10 +'a' ) )
				{
				throw int_precision::bad_int_syntax();
				return false;
				}
			}
	return true;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Dec/2021
//	@brief 		std::string remove_separators from string
//	@return 	std::string	-	Return the string without seperators '\''
//	@param		"s"	-		String of digits to remove separators from
//
//	@todo  Add to do things	
//
// Description:
//	Remove thousand seperators from string
//	A thousand separator is either a ' or a space  ' '
//
static std::string remove_separators(const std::string& s)
	{
	std::string r;
	std::string::const_iterator pos = s.begin();
	r.reserve(s.length() + 16);
	for (; pos != s.end(); ++pos)
		if (*pos != '\'' && *pos != ' ')
			r.push_back(*pos);

	return r;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		31/Oct/2021
//	@brief 		std::string uitostring10
//	@return 	static std::string	-	Return the ascii representation of number
//	@param		"value"	-	Value to convert to ascii string based on RADIX
//	@param		"minlength"	-	minlength of string. default to 0
//
//	@todo  Add to do things	
//
// Description:
//   This convert an uintmax_t (64bit) integer into a decimal string in base 10
//
static inline std::string uitostring10(const uintmax_t value, const unsigned minlength=0 )
	{
	std::string s;
	unsigned digit;
	uintmax_t uvalue=value;

	do
		{
		digit = (unsigned)(uvalue % BASE_10 );
		uvalue /= BASE_10;
		s.push_back((char)ICHARACTER10((unsigned char)digit));
	} while (uvalue > 0);

	while (s.length() < minlength)
		s.push_back('0');
	reverse(s.begin(), s.end());
	return s;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		31/Oct/2021
//	@brief 		std::string uitostringbase
//	@return 	static std::string	-	Return the ascii representation of number
//	@param		"value"	-	Value to convert to ascii string based on RADIX
//	@param		"radix"		-	base of number
//	@param		"minlength"	-	minlength of string. default to 0
//
//	@todo  Add to do things	
//
// Description:
//   This convert an uintmax_t (64bit) integer into a string in base radix with a minimum lenght of minlength
//		padded with 0
//	The string is returned in reverse order
//
static inline std::string uitostringbase(const uintmax_t value, const unsigned radix, const unsigned minlength = 0)
	{
	std::string s;
	unsigned digit;
	uintmax_t uvalue = value;

	if (radix < BASE_2 || radix > 36)
		return s;  // Conversion not supported

	do
		{
		digit = (unsigned)(uvalue % radix);
		uvalue /= radix;
		if (digit < 10)
			s.push_back((unsigned char)( digit + '0'));
		else
			s.push_back((unsigned char)(digit - 10 + 'a'));
	} while (uvalue > 0);

	while (s.length() < minlength)
		s.push_back('0');
	//reverse(s.begin(), s.end());
	return s;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		2/3/2006
//	@brief 		std::string itostring
//	@return 	static std::string	-	Return the ascii representation of number
//	@param		"value"	-	Value to convert to ascii string based on RADIX
//	@param		"radix"	-	RADIX value of conversion 
//
//	@todo  Add to do things	
//
// Description:
//   This function replace Microsoft _itoa() to a generic function that return the
//   string representation of the number in Base Radix. 
//   Radix can be in the range from 2..256 (only 2..36 deliveres a readable string)
//   only if value is < 0 will it return with a leading sign
std::string itostring( const int value, const unsigned radix )
   {
   std::string s;
   unsigned digit;
   unsigned uvalue;
   
   if (radix < BASE_2 || radix > 36 )
      return s;  // Conversion not supported

   if( radix == BASE_10 && value < 0 )
      uvalue = -value;
   else
      uvalue = (unsigned)value;

   do 
      {
      digit = (unsigned) (uvalue % radix);
      uvalue /= radix;

      if( radix <= 36 )
         {
         // Convert to ascii and store
         if( digit < 10 )
            s.push_back( (char)ICHARACTER10( (unsigned char)digit ) );      
         else
            s.push_back( (char)( digit - 10 + 'a' ) );      
         }
      else
         { // Keep it 'binary' not readable string
         s.push_back( (unsigned char)digit );      
         }
   } while (uvalue > 0);

   if( radix == BASE_10 && value < 0 )
      s.push_back( '-' );
   reverse(s.begin(), s.end());
   return s;
   }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		std::string _int_precision_atoip
//	@return 	vector<iptype>	-	The integer precision number
//	@param		"str"	-	The arbitrary precision string as a std::string
// @param		"sign"	-	Returned the sign as either +1 or -1
//
//	@todo  
//
// Description:
// Convert ascii string to vector<iptype> number
// A leading 0 is intepreted as a octal number
// a leading 0x is interpreted as a hexadecimal number 
// a leading 0b is interpreted as a binary number
// otherwise it's a decimal number.
//
std::vector<iptype> _int_precision_atoip(const std::string& str, int *sign)
	{
	std::vector<iptype> number;

	number = _int_precision_atoip(str.c_str(), sign);
	return number;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_atoi
//	@return 	string	-	The integer precision string
//	@param		"s"	-		The arbitrary precision string as char *
//	@param		"sign"	-	Returned the sign as either +1 or -1
//
//	@todo  
//
// Description:
// Convert ascii string to a binary number
// A leading 0 is intepreted as a octal number (BASE_8)
// a leading 0x is interpreted as a hexadecimal number  (BASE_16)
// a leading 0b is interpreted as a binary number (BASE_2)
// otherwise it's a decimal number.(BASE_10)
// The resulting number is a vector<iptype> type
//
std::vector<iptype> _int_precision_atoip(const char *str, int *sign)
	{
	std::string s(str);
	std::string::const_iterator pos;
	std::vector<iptype> number(1, 0);
	int base = BASE_10;

	s = remove_separators(s);
	number.reserve(s.size() + 16);
	*sign = +1;
	pos = s.begin();
	if (*pos == '+' || *pos == '-')
		{
		*sign = CHAR_SIGN(*pos);
		++pos;
		if (pos == s.end())
			throw int_precision::bad_int_syntax();
		}

	if (*pos == '0') // Octal, binary or hex representation
		{
		if (pos + 1 != s.end() && tolower(pos[1]) == 'x')
			{
			base = BASE_16;
			pos += 2;
			}
		else
			if (pos + 1 != s.end() && tolower(pos[1]) == 'b')
				{// Collec Binary representation
				base = BASE_2;
				pos += 2;
				}
			else
				{ // Collect octal representation
				base = BASE_8;
				pos += 1;
				}
		}

	check_digits(s, pos - s.begin(), s.size(), base);
	if(base == BASE_10)
		number = string2number(s, pos - s.begin(), s.size());
	else
		if(base==BASE_2 || base==BASE_16)
			number = stringbase2number(s, pos - s.begin(), s.size(),base);
		else
			{// BASE_8
			number = stringbase8_2number(s, pos - s.begin(), s.size() );
			}

	if (number.size() == 1 && number[0] == 0 && *sign == -1)
		*sign = +1;

	return number;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		9/Aug/2021
//	@brief 		std::vector<iptype>  _int_precision_itoa
//	@return 	std::string -	Return the inter precision as a string (no leading sign)
//	@param      "a"	-	the internal integer precision string
// @param		"base" -	the base to convert to. default is BASE_10
//
//	@todo 
//
// Description:
//   Convert int_precsion to ascii string using base
//   The string returned has no leading sign
//
std::string _int_precision_itoa(const std::vector<iptype>  *a, const int base )
	{
	std::vector<iptype> src, tmp_rem;
	std::string s;
	size_t i;

	src = *a;
	s.reserve(src.capacity());
	if (src.size() == 1 && base==BASE_10)  // Only one iptype binary number in vector. Convert it directly to string and return
		return s = uitostring10((uintmax_t)src[0]);

	switch (base)
		{
		case BASE_2:
			for (i = 0; i < src.size()-1; ++i)
				s+=uitostringbase(src[i], BASE_2, 64);
			s+=uitostringbase(src[i], BASE_2);
			s += "b0";
			break;
		case BASE_16:
			for (i = 0; i < src.size() - 1; ++i)
				s += uitostringbase(src[i], BASE_16, 16);
			s += uitostringbase(src[i], BASE_16);
			s += "x0";
			break;
		case BASE_10: // Base 10 only. Loop until we have a max of one iptype number left
			{std::vector<iptype> cbasepower10_div(1, _powerof10Table[9]);// 1E9
			for (; src.size() > 1;)
				{ // Take 9 decimal digits at a time
				src = _int_precision_udivrem(&src, &cbasepower10_div, &tmp_rem);
				s += uitostringbase((uintmax_t)tmp_rem[0], BASE_10, 9);
				}
			if (src[0] != 0)
				{// Do the last iptype directly via "native" functions. yielding up to 18 digits, instead of a single digit
				s += uitostringbase((uintmax_t)src[0], BASE_10);
				}
			}
			break;
		case BASE_8:
			{std::vector<iptype> cbase(1, 1 << 30);
			for (; src.size() > 1;)
				{// Take 10 octal digit at  time. 
				src = _int_precision_udivrem(&src, &cbase, &tmp_rem);
				s += uitostringbase((uintmax_t)tmp_rem[0], BASE_8, 10);
				}	
			if (src[0] != 0)
				{// Do the last iptype directly via "native" functions. yielding up to 21 octal digits
				s += uitostringbase((uintmax_t)src[0], BASE_8);
				}
			s += "0";
			}
			break;
		default:
		  // All other bases
		  std::vector<iptype> cbase(1, base);
		  for (; src.size() > 1 || src[0] != 0;)
				{// Take one digit at  time. 
				tmp_rem = _int_precision_urem(&src, &cbase);
				src = _int_precision_udiv(&src, &cbase);
				if (base == BASE_16 && tmp_rem[0] >= BASE_10)
					s.push_back((char)((unsigned char)(tmp_rem[0] - BASE_10 + 'a')));
				else
					s.push_back((char)ICHARACTER10((unsigned char)tmp_rem[0]));
				}
			break;
		}

	if (s.size() == 0)
		s.push_back((char)ICHARACTER10((unsigned char)0));
	reverse(s.begin(), s.end());

	return s;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		9/Aug/2021, 2/Sep/2021
//	@brief 		std::vector<iptype>  _int_precision_itoa
//	@return 	std::string -	Return the inter precision as a string (no leading sign)
//	@param      "a"	-	the internal integer precision string
// @param		"base" -	the base to convert to. default is BASE_10
//
//	@todo 
//
// Description:
//   Convert int_precsion to ascii string using base
//   The string returned has no leading sign
//	Change from a pointer parameter to a reference parameter
//
std::string _int_precision_itoa(const std::vector<iptype>& a, const int base)
	{
	std::vector<iptype> src, tmp_rem;
	std::string s;
	size_t i;

	src = a;
	s.reserve(src.capacity());
	if (src.size() == 1 && base == BASE_10)  // Only one iptype binary number in vector. Convert it directly to string and return
		return s = uitostring10((uintmax_t)src[0]);

	switch (base)
		{
		case BASE_2:
			for (i = 0; i < src.size() - 1; ++i)
				s += uitostringbase(src[i], BASE_2, 64);
			s += uitostringbase(src[i], BASE_2);
			s += "b0";
			break;
		case BASE_16:
			for (i = 0; i < src.size() - 1; ++i)
				s += uitostringbase(src[i], BASE_16, 16);
			s += uitostringbase(src[i], BASE_16);
			s += "x0";
			break;
		case BASE_10: // Base 10 only. Loop until we have a max of one iptype number left
			{std::vector<iptype> cbasepower10_div(1, _powerof10Table[9]);// 1E9
			for (; src.size() > 1;)
				{ // Take 9 decimal digits at a time
				src = _int_precision_udivrem(&src, &cbasepower10_div, &tmp_rem);
				s += uitostringbase((uintmax_t)tmp_rem[0], BASE_10, 9);
				}
			if (src[0] != 0)
				{// Do the last iptype directly via "native" functions. yielding up to 18 digits, instead of a single digit
				s += uitostringbase((uintmax_t)src[0], BASE_10);
				}
			}
			break;
		case BASE_8:
			{std::vector<iptype> cbase(1, 1 << 30);
			for (; src.size() > 1;)
				{// Take 10 octal digit at  time. 
				src = _int_precision_udivrem(&src, &cbase, &tmp_rem);
				s += uitostringbase((uintmax_t)tmp_rem[0], BASE_8, 10);
				}
			if (src[0] != 0)
				{// Do the last iptype directly via "native" functions. yielding up to 21 octal digits
				s += uitostringbase((uintmax_t)src[0], BASE_8);
				}
			s += "0";
			}
			break;
		default:
			// All other bases
			std::vector<iptype> cbase(1, base);
			for (; src.size() > 1 || src[0] != 0;)
				{// Take one digit at  time. 
				tmp_rem = _int_precision_urem(&src, &cbase);
				src = _int_precision_udiv(&src, &cbase);
				if (base == BASE_16 && tmp_rem[0] >= BASE_10)
					s.push_back((char)((unsigned char)(tmp_rem[0] - BASE_10 + 'a')));
				else
					s.push_back((char)ICHARACTER10((unsigned char)tmp_rem[0]));
				}
			break;
		}

	if (s.size() == 0)
		s.push_back((char)ICHARACTER10((unsigned char)0));
	reverse(s.begin(), s.end());
	return s;
	}



///////////////////////////////////////////////
//
//
//    Miscellaneous function
//		abs()
//		ipow()
//		ipow_modulo()
//		iprime()
//		gcd()
//		lcm()
//
//
///////////////////////////////////////////////

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/25/2012
//	@brief 		Calculate abs(x)
//	@return 	int_precision -	Return absolute value of x
//	@param      "x"	- The argument
//
//	@todo  
//
// Description:
//   int precision abs()
//    
//
int_precision abs(const int_precision& x)
	{
	int_precision i;

	if (x.sign() < 0)
		{
		i = x; i.change_sign();
		return i;
		}

	return x;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		26/Aug/2007
//	@brief 		return the integer power of x^y
//	@return 	int_precision	-	The integer precision power of x^y
//	@param		"x"	-	The int precision x
//	@param		"y"	-	The int precision y. Max y is 2^32-1
//
//	@todo  
//
// Description:
// Return the integer power of x^y. For any pratical purpose the power y is restricted to 2^32-1
//
int_precision ipow( const int_precision& x, const int_precision& y )
   {
   int_precision p(x), r(1);

   for(uintmax_t n = (uintmax_t)y; n > 0; n >>= 1) 
        {
        if( ( n & 0x1 ) != 0 ) r *= p;  // Odd
        if( n > 1 ) p *= p;				// Square it				 
        }
   return r;
   }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		24/Aug/2012
//	@brief 		return the integer power of x^y%z
//	@return 	int_precision	-	The integer precision power of x^y%z
//	@param		"x"	-	The int precision x
//	@param		"y"	-	The int precision y. Max y is 2^32-1
// @param		"z"	-	The int precision z.
//	@todo  
//
// Description:
// Return the integer power of x^y%z. For any pratical purose the power y is restricted to 2^32-1
//
int_precision ipow_modulo( const int_precision& x, const int_precision& y, const int_precision& z )
   {
   int_precision p(x), r(1);

   p%=z;
   for(int n = y; n > 0; n >>= 1) 
        {
        if( ( n & 0x1 ) != 0 ) { r *= p; r %= z; } // Odd
        p *= p;	p %= z;					 
        }
   return r;
   }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		2/Sep/2012
//	@brief 		Check a number for a prime
//	@return 	bool-	true is the integer is a prime number false otherwise
//	@param		"prime"	-	The int precision prime
//	@todo  
//
// Description:
// Return true if the integer prime is a prime number.
// All integers are of the form 30k + i for i = 0, 1, 2,...,29 and k an integer from 0..  However, 2 divides 0, 2, 4,...,28 and 3 divides 0, 3, 6,...,27 and 5 divides 0, 5, 10,...,25. 
// So all prime numbers are of the form 30k + i for i = 1, 7, 11, 13, 17, 19, 23, 29 (i.e. for i < 30 such that gcd(i,30) = 1). 
// Note that if i and 30 are not coprime, then 30k + i is divisible by a prime divisor of 30, namely 2, 3 or 5, and is therefore not a prime.
// 
//
bool iprime(const int_precision& prime)
	{
	int precheck[11] = { 10, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 };
	int primes[9] = { 8, 1, 7, 11, 13, 17, 19, 23, 29 };
	int_precision count, kp(30), mod;
	int i;

	for (i = 1; i <= precheck[0]; i++)
	if ((int)(prime % (int_precision)precheck[i]) == 0) return prime==int_precision(precheck[i]);

	for (; kp * kp < prime; kp += 30)   //Loop to divide the number by every number 6*count-1 and 6*count+1 and count < sqrt(i)
		{
		for (i = 1; i <= primes[0]; i++)
			{
			count = kp + primes[i];
			mod = prime % count;
			if (mod == 0)				// Statement to change the variable 'check' to 1 if the number gets divided
				return false;			//meaning its not prime
			}
		}

	return true;						// It is a prime
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		30/Aug/2021
//	@brief 		Greatest Common Divisor (gcd)
//	@return 	int_precision - gcd(a,b)
//	@param		"a"	-	first number - positive
// @param		"b"	-	second number - positive
//
//	@todo  
//
// Description:
// Return the greatest common divisor of the two numbers a & b.
// It used the Binary gcd method only using shifting and subtraction
// Changed to also handle negative arguments a,b;
//	Algorithm: while (b>0) { tmp = b; b = a%b; a = tmp; } return a;
//
int_precision gcd( const int_precision& a, const int_precision& b )
	{
	size_t shift;
	int_precision u, v, tmp, i0(0), i1(1);

	// GCD(0,v)==v; GCD(u,0)==0; GCD(0,0)==0
	if (a.iszero() ) return b;
	if (b.iszero() ) return a;
	u = a; v = b; 
	if(u < i0) u = -u; 
	if(v < i0) v = -v;
#if false
	while (b > i0)
		{
		tmp = v; v = u%v; u = tmp;
		}
	return u;
#else
		tmp = u | v;
		shift = tmp.ctz();
		u >>=u.ctz();
		do {
			v >>= v.ctz();
			if (u > v) {
				tmp = v; v = u; u = tmp;
				}
			v -= u;
		} while (v != i0);
		return u << shift;
	}
#endif


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/Feb/2017
//	@brief 		Least Common Multiplier
//	@return 	int_precision - lcm(a,b)
//	@param		"a"	-	first number
// @param		"b"	-	second number
//
//	@todo  
//
// Description:
// Return the least common multiplier of the two numbers a & b.
// It used the Binary lcm method only using shifting, subtraction and one multiplication and one division
// 
//
int_precision lcm(const int_precision& a, const int_precision& b)
	{
	int_precision r(a), gcd_ab;

	gcd_ab = gcd(a, b);
	r /= gcd_ab;
	r *= b;
	return r;
	}



///////////////////////////////////////////////
//
//
//    End of Integer Precision Core
//
//
///////////////////////////////////////////////

///////////////////////////////////////////////
//
//
//    Floating point Precision Core
//
//
///////////////////////////////////////////////

class float_precision_ctrl float_precision_ctrl(PRECISION,ROUND_NEAR);

///////////////////////////////////////////////
//
//
//    Floating point Precision Input, Output operator
//
//
///////////////////////////////////////////////


std::ostream& operator<<( std::ostream& strm, const float_precision& d )
    { return strm << _float_precision_fptoa( const_cast<float_precision *>(&d) ).c_str();}

std::istream& operator>>( std::istream& strm, float_precision& d )
         { char ch; std::string s; int cnt, exp_cnt=0;
         strm >> ch;  while( ch == ' ' ) strm.get(ch);  // Ignore leading white space.
         if( ch == '+' || ch == '-' ) { s += ch; strm >> ch; } else s += '+';  // Parse sign
         for( cnt = 0; ch >= '0' && ch <= '9'; cnt++, strm >> ch ) s += ch;  // Parse integer part
         if( ch == '.' )  // Parse fraction
            for( s += '.', strm >> ch; ch >= '0' && ch <= '9'; cnt++, strm >> ch ) s += ch;   // Parse fraction part
         if( ch == 'e' || ch == 'E' )
            {
            s += 'e'; strm >> ch; if( ch == '+' || ch == '-' ) { s += ch; strm >> ch; } else s += '+';  // Parse Expo sign 
            for( exp_cnt =0; ch >= '0' && ch <= '9'; exp_cnt++, strm >> ch ) s += ch;  // Parse expo number
            }

         std::cin.putback( ch );  // ch contains the first character not part of the number, so put it back
         if( !strm.fail() && ( cnt > 0 || exp_cnt > 0 ) )  // Valid number 
            d = float_precision( const_cast<char *>( s.c_str() ), float_precision_ctrl.precision(), float_precision_ctrl.mode() );
         return strm;
         }


///////////////////////////////////////
//
// CONVERT BINARY FLOAT PRECISION to and from ascii representation
//
//   _float_precision_fptoa()		-- Convert from float_precision to string format
//	  _float_precision_fptoainteger()-- Convert from float_precision to integer string format. Obsolete?
//   _float_precision_atofp()		-- Convert from string format to float_precision format
//	  _float_precision_fptoip()		-- Convert from float_precision to int_precision
//	  static check_float_digits()	-- Check string for illegal number characers
//
///////////////////////////////////////

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		25/Oct/2021
//	@brief 		std::string check_'float_digits
//	@return 	bool	-	Return true if digits is valid other throw an exception and return false
//	@param		"s"	-		String of digits to check
//	@param		"start"	-	Value to convert to ascii string based on RADIX
//	@param		"end"	-	RADIX value of conversion 
//
//	@todo  Add to do things	
//
// Description:
//		Check the ascii string for valid digits according to base.
//
static bool check_float_digits(const std::string& s, const size_t start, const size_t end, const int base=BASE_10 )
	{
	for (size_t i = start; i < end; i++)
		{
		if (base <= BASE_10 && (s[i] < '0' || s[i] > (base - 1 + '0')))
			{
			throw float_precision::bad_float_syntax();
			//return false;
			}
		else
			if (!((s[i] >= '0' && s[i] <= '9') || (s[i] >= 'a' && s[i] <= 'a' + base - 1)))
				{
				throw float_precision::bad_float_syntax();
				//return false;
				}
		}
	return true;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		4/Dec/2021
//	@brief 		Convert part of float_precision numbers into string (decimal representation)
//	@return		std::string -	The partly converted decimal floating point string	
//	@param		"fp"		-	float_precision number to partly convert
//	@param		"digits"	-	The remainf digis in the float_precision number 
//
//	@todo 	
//
// Description:
//   Convert partly a float_precision numbers into string (decimal representation)
//		A call to the function convert the next up to max_digits decimal number into a string
//		and return the string
//		Digits to be converted is inthe range 1..max_digits
//
static std::string number2Decimal(float_precision& fp, const int digits )
	{
	uintmax_t di;
	int min_width = digits;						// So many wanted digits. always >= 1

	if (min_width <= 0) min_width = 1;					// at least 1
	if (min_width > MAX_DECIMAL_DIGITS) min_width = MAX_DECIMAL_DIGITS;	// Take max 18 digits at a time
	fp *= _fpPowerof10Table[min_width];			// float_precision(_powerof10table[min_width]);
	di = (uintmax_t)fp.toFraction();  			// Get the next max digit decimal number										
	fp.precision(fp.precision() - min_width);	// Reduce precision with min_width decimal digits
	return uitostring10(di, min_width);			// Return min_width length string
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		4/Dec/2021
//	@brief 		Convert a trunk of float_precision numbers into string (decimal representation)
//	@return		std::string -	The trunk converted to a decimal floating point string	
//	@param		"fp"		-	float_precision number to partly convert
//
//	@todo 	
//
// Description:
//   Convert a trunk of a float_precision numbers into string (decimal representation)
//		A call to the function convert the next trunk size decimal number into a string
//		and return the trunk string
//		A trunk size is defined as threshold*max_digits number of decimal digits
//

static std::string trunk2Decimal(float_precision& fp)
	{
	int i;
	std::string str;
	float_precision trunk(fp); 
	static float_precision _trunkPowerof10(0, MAX_TRUNK_SIZE * 20, ROUND_DOWN);

	if (_trunkPowerof10.iszero())// is _trunkPowerof10 build or created
		{
		float_precision p(_powerof10Table[MAX_DECIMAL_DIGITS], fp.precision(), fp.mode());
		for (i = MAX_TRUNK_SIZE, _trunkPowerof10 = float_precision(1); i > 0; i >>= 1 )	
			{// Build multiply factor for trunk size
			if (i & 1) _trunkPowerof10 *= p;	// Odd
			if (i > 1) p *= p;					// square it
			}
		}

	// Set trunk precision to reduce precision to at least hold the needed number of digits we can ge out of the trunk size
	trunk.precision(MAX_TRUNK_SIZE * 20);	// We overallocate the precision 20 versus 18 digits in the trunk
	str.reserve(MAX_TRUNK_SIZE * 20);
	// Do it trunk size of thr (eshold) uintmax_t digits
	for (i = 0; i < MAX_TRUNK_SIZE; ++i)
		{
		str += number2Decimal(trunk, MAX_DECIMAL_DIGITS);
		}
	fp *= _trunkPowerof10;
	fp.toFraction();
	fp.precision(fp.precision() - MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS);
	return str;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		4/Dec/2021
//	@brief 		Convert a Megatrunk of float_precision numbers into string (decimal representation)
//	@return		std::string -	The trunk converted to a decimal floating point string	
//	@param		"fp"		-	float_precision number to partly convert
//
//	@todo 	
//
// Description:
//   Convert a trunk of a float_precision numbers into string (decimal representation)
//		A call to the function convert the next trunk size decimal number into a string
//		and return the trunk string
//		A trunk size is defined as threshold*max_digits number of decimal digits
//

static std::string kilotrunk2Decimal(float_precision& fp)
	{
	int i;
	std::string str;
	float_precision trunk(fp);
	static float_precision _kilotrunkPowerof10(0, MAX_KILOTRUNK_SIZE * MAX_TRUNK_SIZE * 20, ROUND_DOWN);

	if (_kilotrunkPowerof10.iszero())// is _megatrunkPowerof10 build or created
		{
		float_precision p(_powerof10Table[MAX_DECIMAL_DIGITS], fp.precision(), fp.mode());
		for (i = MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE, _kilotrunkPowerof10 = float_precision(1); i > 0; i >>= 1)
			{// Build multiply factor for trunk size
			if (i & 1) _kilotrunkPowerof10 *= p;	// Odd
			if (i > 1) p *= p;						// square it
			}
		}

	// Set trunk precision to reduce precision to at least hold the needed number of digits we can ge out of the trunk size
	trunk.precision(MAX_KILOTRUNK_SIZE * MAX_TRUNK_SIZE * 20);	// We overallocate the precision 20 versus 18 digits in the trunk
	str.reserve(MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE * 20);
	// Do it trunk size of thr (eshold) uintmax_t digits
	for (i = 0; i < MAX_KILOTRUNK_SIZE; ++i)
		{
		str += trunk2Decimal(trunk);
		}
	fp *= _kilotrunkPowerof10;
	fp.toFraction();
	fp.precision(fp.precision() - MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS);
	return str;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		4/Dec/2021
//	@brief 		Convert a Megatrunk of float_precision numbers into string (decimal representation)
//	@return		std::string -	The trunk converted to a decimal floating point string	
//	@param		"fp"		-	float_precision number to partly convert
//
//	@todo 	
//
// Description:
//   Convert a trunk of a float_precision numbers into string (decimal representation)
//		A call to the function convert the next trunk size decimal number into a string
//		and return the trunk string
//		A trunk size is defined as threshold*max_digits number of decimal digits
//
static std::string megatrunk2Decimal(float_precision& fp)
	{
	int i;
	std::string str;
	float_precision trunk(fp);
	static float_precision _megatrunkPowerof10(0, MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE * MAX_TRUNK_SIZE * 20, ROUND_DOWN);

	if (_megatrunkPowerof10.iszero())// is _megatrunkPowerof10 build or created
		{
		float_precision p(_powerof10Table[MAX_DECIMAL_DIGITS], fp.precision(), fp.mode());
		for (i = MAX_MEGATRUNK_SIZE*MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE, _megatrunkPowerof10 = float_precision(1); i > 0; i >>= 1)
			{// Build multiply factor for trunk size
				if (i & 1) _megatrunkPowerof10 *= p;	// Odd
				if (i > 1) p *= p;						// square it
			}
		}

	// Set trunk precision to reduce precision to at least hold the needed number of digits we can ge out of the trunk size
	trunk.precision(MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE * 20);	// We overallocate the precision 20 versus 18 digits in the trunk
	str.reserve(MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE * 20);
	// Do it trunk size of thr (eshold) uintmax_t digits
	for (i = 0; i < MAX_MEGATRUNK_SIZE; ++i)
		{
		str += kilotrunk2Decimal(trunk);
		}
	fp *= _megatrunkPowerof10;
	fp.toFraction();
	fp.precision(fp.precision() - MAX_MEGATRUNK_SIZE*MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS);
	return str;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Dec/2021
//	@brief 		Convert a Mega10trunk of float_precision numbers into string (decimal representation)
//	@return		std::string -	The trunk converted to a decimal floating point string	
//	@param		"fp"		-	float_precision number to partly convert
//
//	@todo 	
//
// Description:
//   Convert a trunk of a float_precision numbers into string (decimal representation)
//		A call to the function convert the next trunk size decimal number into a string
//		and return the trunk string
//		A trunk size is defined as threshold*max_digits number of decimal digits
//
static std::string mega10trunk2Decimal(float_precision& fp)
	{
	int i;
	std::string str;
	float_precision trunk(fp);
	static float_precision _mega10trunkPowerof10(0, MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE * MAX_TRUNK_SIZE * 20, ROUND_DOWN);

	if (_mega10trunkPowerof10.iszero())// is _gigatrunkPowerof10 build or created
		{
		float_precision p(_powerof10Table[MAX_DECIMAL_DIGITS], fp.precision(), fp.mode());
		_mega10trunkPowerof10 = float_precision(1);
		for (i = MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE*MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE; i > 0; i >>= 1)
			{// Build multiply factor for trunk size
			if (i & 1) _mega10trunkPowerof10 *= p;	// Odd
			if (i > 1) p *= p;						// square it
			}
		}

	// Set trunk precision to reduce precision to at least hold the needed number of digits we can ge out of the trunk size
	trunk.precision(MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE * 20);	// We overallocate the precision 20 versus 18 digits in the trunk
	str.reserve(MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE * 20);
	// Do it trunk size of thr (eshold) uintmax_t digits
	for (i = 0; i < MAX_MEGA10TRUNK_SIZE; ++i)
		{
		str += megatrunk2Decimal(trunk);
		}
	fp *= _mega10trunkPowerof10;
	fp.toFraction();
	fp.precision(fp.precision() - MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE*MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS);
	return str;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Dec/2021
//	@brief 		Convert a Mega100trunk of float_precision numbers into string (decimal representation)
//	@return		std::string -	The trunk converted to a decimal floating point string	
//	@param		"fp"		-	float_precision number to partly convert
//
//	@todo 	
//
// Description:
//   Convert a trunk of a float_precision numbers into string (decimal representation)
//		A call to the function convert the next trunk size decimal number into a string
//		and return the trunk string
//		A trunk size is defined as threshold*max_digits number of decimal digits
//
static std::string mega100trunk2Decimal(float_precision& fp)
	{
	int i;
	std::string str;
	float_precision trunk(fp);
	static float_precision _mega100trunkPowerof10(0, MAX_MEGA100TRUNK_SIZE*MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE * MAX_TRUNK_SIZE * 20, ROUND_DOWN);

	if (_mega100trunkPowerof10.iszero())// is _gigatrunkPowerof10 build or created
		{
		float_precision p(_powerof10Table[MAX_DECIMAL_DIGITS], fp.precision(), fp.mode());
		_mega100trunkPowerof10 = float_precision(1);
		for (i = MAX_MEGA100TRUNK_SIZE*MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE*MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE; i > 0; i >>= 1)
			{// Build multiply factor for trunk size
			if (i & 1) _mega100trunkPowerof10 *= p;	// Odd
			if (i > 1) p *= p;						// square it
			}
		}

	// Set trunk precision to reduce precision to at least hold the needed number of digits we can ge out of the trunk size
	trunk.precision(MAX_MEGA100TRUNK_SIZE*MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE * 20);	// We overallocate the precision 20 versus 18 digits in the trunk
	str.reserve(MAX_MEGA100TRUNK_SIZE*MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE * 20);
	// Do it trunk size of thr (eshold) uintmax_t digits
	for (i = 0; i < MAX_MEGA100TRUNK_SIZE; ++i)
		{
		str += mega10trunk2Decimal(trunk);
		}
	fp *= _mega100trunkPowerof10;
	fp.toFraction();
	fp.precision(fp.precision() - MAX_MEGA100TRUNK_SIZE*MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE*MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS);
	return str;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Dec/2021
//	@brief 		Convert a Gigatrunk of float_precision numbers into string (decimal representation)
//	@return		std::string -	The trunk converted to a decimal floating point string	
//	@param		"fp"		-	float_precision number to partly convert
//
//	@todo 	
//
// Description:
//   Convert a trunk of a float_precision numbers into string (decimal representation)
//		A call to the function convert the next trunk size decimal number into a string
//		and return the trunk string
//		A trunk size is defined as threshold*max_digits number of decimal digits
//
static std::string gigatrunk2Decimal(float_precision& fp)
	{
	int i;
	std::string str;
	float_precision trunk(fp);
	static float_precision _gigatrunkPowerof10(0, MAX_GIGATRUNK_SIZE*MAX_MEGA100TRUNK_SIZE*MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE * MAX_TRUNK_SIZE * 20, ROUND_DOWN);

	if (_gigatrunkPowerof10.iszero())// is _gigatrunkPowerof10 build or created
		{
		float_precision p(_powerof10Table[MAX_DECIMAL_DIGITS], fp.precision(), fp.mode());
		_gigatrunkPowerof10 = float_precision(1);
		for (i = MAX_GIGATRUNK_SIZE*MAX_MEGA100TRUNK_SIZE*MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE*MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE; i > 0; i >>= 1)
			{// Build multiply factor for trunk size
			if (i & 1) _gigatrunkPowerof10 *= p;	// Odd
			if (i > 1) p *= p;						// square it
			}
		}

	// Set trunk precision to reduce precision to at least hold the needed number of digits we can ge out of the trunk size
	trunk.precision(MAX_GIGATRUNK_SIZE*MAX_MEGA100TRUNK_SIZE*MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE * 20);	// We overallocate the precision 20 versus 18 digits in the trunk
	str.reserve(MAX_GIGATRUNK_SIZE*MAX_MEGA100TRUNK_SIZE*MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE * 20);
	// Do it trunk size of thr (eshold) uintmax_t digits
	for (i = 0; i < MAX_GIGATRUNK_SIZE; ++i)
		{
		str += mega100trunk2Decimal(trunk);
		}
	fp *= _gigatrunkPowerof10;
	fp.toFraction();
	fp.precision(fp.precision() - MAX_GIGATRUNK_SIZE*MAX_MEGA100TRUNK_SIZE*MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE*MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS);
	return str;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		14/Dec/2021
//	@brief 		Calculate the expoenent reduction for very small numbers.
//	@return		int -	The decimal 10 exponent reductions
//	@param		"f"		-	float_precision number to partly convert
//
//	@todo 	
//
// Description:
//   The number f much be less than 1. e.g. exponent() < 0 
//	 This claulate the decimal exponent reducion from he base 2 internal exponent
//
static int exponent2reduction(float_precision& f)
	{
	int expo10reduction = 0;
	eptype expo = f.exponent();  // start with the power of 2 exponent
	static float_precision e1kfactor(0, 1'000 + 2);
	static float_precision e10kfactor(0, 10'000 + 2);
	static float_precision e100kfactor(0, 100'000 + 2);
	float_precision p(10, e1kfactor.precision());

	if (expo < 0)
		{
		expo10reduction = -(int)(expo / log2(10));
		if (expo10reduction != 0)
			expo10reduction -= 1;  // Calculate the exponent decimal 10 reduction that ensure that 10^expo10reduction< 2^exponent
		expo = expo10reduction; // Now it is converted to a power of 10 exponent
		f.precision(f.precision() /*+ 5*/);  // Add extra guard bits
											 //Need to biuld any factors?
		if (expo >= 1'000 && e1kfactor.iszero())
			{// Build the e10factor first time it is needed.
			e1kfactor = float_precision(1); p = float_precision(10);
			for (eptype i = 1'000; i > 0; i >>= 1)
				{// Build multiply factor for trunk size
				if (i & 1) e1kfactor *= p; // Odd
				if (i > 1) p *= p;			// square it
				}
			}
		if (expo >= 10'000 && e10kfactor.iszero())
			{// Build the e10factor first time it is needed.
			p.precision(e10kfactor.precision());
			e10kfactor = float_precision(1); p = e1kfactor;
			for (eptype i = 10; i > 0; i >>= 1)
				{// Build multiply factor for trunk size
				if (i & 1) e10kfactor *= p; // Odd
				if (i > 1) p *= p;			// square it
				}
			}

		if (expo >= 100'000 && e100kfactor.iszero())
			{// Build the e10factor first time it is needed.
			p.precision(e100kfactor.precision());
			e100kfactor = float_precision(1); p = e10kfactor;
			for (eptype i = 10; i > 0; i >>= 1)
				{// Build multiply factor for trunk size
				if (i & 1) e100kfactor *= p; // Odd
				if (i > 1) p *= p;			// square it
				}
			}

		for (; expo >= 100'000; expo -= 100'000)
			f *= e100kfactor;

		for (; expo >= 10'000; expo -= 10'000)
			f *= e10kfactor;

		for (; expo >= 1'000; expo -= 1'000)
			f *= e1kfactor;

		for (; expo > 0; expo -= std::min((int)MAX_DECIMAL_DIGITS, (int)expo))
			f *= _fpPowerof10Table[std::min((int)MAX_DECIMAL_DIGITS, (int)expo)];

		f.precision(f.precision() /*- 5*/);  // reverse the extra guard bits
		}

	return expo10reduction;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		2/Dec/2021
//	@brief 		Convert float_precision numbers into string (decimal representation)
//	@return		std::string -	The decimal floating point string	
//	@param		"a"			-	float_precision number to convert
//
//	@todo 	
//
// Description:
//   Convert float_precision numbers into string (decimal representation)
//   We do it in multiple steps.
//	 1)	Repeat doing a megatrunk size by extrating a megatrunk size from the fptoa number
//		in the same manner as step 2
//   2) Repeat doing a trunk size by extracting a trunk size from fptoa number 
//		Then repeat multiple group within the trunk size of max_digis number at the time
//		The most efficient trunk size has be found to be 50*max_digits decimal in a trunk size
//		For every loop we reduce the original fptoa precision  with the trunk size precision 
//	 3)	Take the remaning fpto number in the range 0..max_digits
//
std::string _float_precision_fptoa(const float_precision *a)
	{
	const size_t thr = MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS;
	const size_t kthr = thr * MAX_KILOTRUNK_SIZE;
	const size_t mthr = kthr * MAX_MEGATRUNK_SIZE;
	const size_t m10thr = mthr * MAX_MEGA10TRUNK_SIZE;
	const size_t m100thr = m10thr * MAX_MEGA100TRUNK_SIZE;
	const size_t gigathr = m100thr * MAX_GIGATRUNK_SIZE;
	float_precision fracp, intp;
	std::string str;
	int expo10, sign = 1;
	size_t found, len;
	
	if (a->iszero() ) 
		return std::string("0E0");
	sign = a->sign();
	fracp.precision(a->precision());		// ensure fracp and intp has the same precision as a
	intp.precision(a->precision());
	fracp = modf(*a, &intp);				// Separate the integer and fraction 
	intp = fabs(intp);
	str += _float_precision_fptoainteger(&intp);  // Convert integer part to string 
	expo10 = (int)(str.size() - 1);			// Extract Expo as the size of string - 1
	if (!fracp.iszero() || str.size() > 1)	// If fraction part then add "." to string
		str.insert(str.begin() + 1, '.');
	fracp = fabs(fracp);					// Remove any fraction sign if any
	//fracp.precision(fracp.precision()- (int)(fracp.exponent()/log2(10)));  // Adjust for large negative exponent
	fracp.precision(fracp.precision() + 2 );  // Add extra guard bits
	fracp.mode(ROUND_DOWN);

	// Check for large negative exponent to avoid generating leading zero that will be cut of anyway at the end
	if(intp.iszero())
		expo10 -= exponent2reduction(fracp);	// remove large negative exponent from fracp and adjust fracp accordingly
	
	str.reserve(fracp.precision() + str.length() + 32);  // Ensure enough room,for the string to avoid reallocating
	len = fracp.precision() -2 + str.length();// The expected len of the fraction before we have enough
	 
	// Do it in Giga trunks size
	for (; (str.length() + gigathr < len) && !fracp.iszero(); )
		str += gigatrunk2Decimal(fracp);

	// Do it in 100M trunks size
	for (; (str.length() + m100thr < len) && !fracp.iszero(); )
		str += mega100trunk2Decimal(fracp);

	// Do it in 10M trunks size
	for (; (str.length() + m10thr < len) && !fracp.iszero(); )
 		str += mega10trunk2Decimal(fracp);
//-------
	// Do it in mega trunks size
	for (; (str.length() + mthr < len) && !fracp.iszero(); )
		str += megatrunk2Decimal(fracp);
			
	// Do it in kilo trunks size
	for (; (str.length() + kthr < len) && !fracp.iszero(); )
		str += kilotrunk2Decimal(fracp);

	// Do it in trunks of treshold times (64bit numbers). threshold*max_digits
	for (; (str.length() + thr < len) && !fracp.iszero(); )
		str += trunk2Decimal(fracp);

	// Do it for the remaining less than a trunk size of threshold
	for (; (str.length() < len || abs(fracp.exponent()) > 3) && !fracp.iszero(); )
		str += number2Decimal(fracp, (int)(len - str.length()));

	// Remove trailing zeros if any
	/*found = str.find_last_not_of('0');
	if (found != std::string::npos && str.size() > std::max(2, (int)found + 1))
		{
		str.erase(std::max(2, (int)found + 1));
		if (str[str.size() - 1] == '.')
			str.erase(1);
		}*/

	// check and find negative exponent.
	if (str[0] == '0')
		{
		found = str.find_first_not_of('0', 2);
		if (found == std::string::npos)
			{// Everything is zero
			return std::string("0E0");
			}
		else
			{// Reduce exponent with the number of leading zeros found
			expo10 += -(int)(found - 1);
			str.erase(0, found);
			if (str.size()>1)
				str.insert(str.begin() + 1, '.');
			}
		}
	// str is on the form x.yyyyyyyyyyyy
	// Remove any digits exceeding the precision
	if (str.size() > a->precision() + 2)
		str.erase(a->precision() + 2);

	// Remove trailing zeros if any
	found = str.find_last_not_of('0');
	if (found != std::string::npos && str.size() > std::max(2, (int)found + 1))
		{
		str.erase(std::max(2, (int)found + 1));
		if (str[str.size() - 1] == '.')
			str.erase(1);
		}

	// Add E[+-]exponent to back of str and sign if negative to the front
	str = (sign < 0 ? "-" : "") + str+"E"+ itostring(expo10, BASE_10);
	return str;
	}
   
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		7/Sep/2021
//	@brief 		Convert a string decimal number into a float_precision number
//	@return 	std::string - The decimal floating point string	
//	@param		"str"		-	ascii string of floating point number to convert
// @param		"p"			- The precision of the number
// @param		"m"			- The round mode of the number
//
//	@todo 	
//
// Description:
//   Convert ascii string into a float_precision numbers 
//    The ascii float format is based on standard C notation
//	
//
float_precision _float_precision_atofp(const std::string& str, size_t p, enum round_mode m)
	{
	return _float_precision_atofp(str.c_str(), p, m);
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		7/Sep/2021
//	@brief 		Convert a string decimal number into a float_precision number
//	@return 	std::string - The decimal floating point string	
//	@param		"str"		-	ascii string of floating point number to convert
// @param		"p"			- The precision of the number
// @param		"m"			- The round mode of the number
//
//	@todo 	
//
// Description:
//  Convert ascii string into a float_precision numbers 
//  The ascii float format is based on standard C notation
//	Using the string2number function for better performance
//	However hexadecimal float is aso supported
float_precision _float_precision_atofp(const char *str, size_t p, enum round_mode m)
	{
	int sign, sign_expo, base = 10;
	eptype expo_e10 = 0;
	size_t f_digits=0, i_digits=0;
	std::string::size_type i, nidx, idx;
	std::string s(str);
	std::string::iterator pos;
	std::vector<iptype> number(1,0), sum(1,0);
	bool ipart = false, fpart = false, epart = false;

	number.reserve(s.size() + 16);
	s=remove_separators(s);
	idx = 0;
	sign = +1;
	// Parse leading sign if any
	pos = s.begin();
	if (*pos == '+' || *pos == '-')  // 
		{
		sign = CHAR_SIGN(*pos);
		pos++;
		idx = 1;
		if (pos == s.end())
			throw float_precision::bad_int_syntax();
		}

	// Determine any significant, fraction sign or exponent sign
	nidx = s.find_first_of(".eE", idx);
	if (nidx == std::string::npos) // Only digits (INTEGER) if any
		{
		int_precision ip(str);			// Construct integer to int_precision
		float_precision fp(ip, p, m );	// Construct float_precision from int_precision 
		return fp;
		} // End of Integer parsing

	// Check for hexadecimal or decimal constant
	if (pos[0] == '0' && pos + 1 != s.end() && tolower(pos[1]) == 'x')
		{
		base = 16; 
		idx += 2; 
		}

	// Floating point number starts here
	// Pick up significant beteen idx and nidx 
	if (nidx > idx) // Number of digits before the . sign or exponent Ee
		{
		ipart = true;
		// Strip leading zeros
		for (i = idx; i != nidx; i++) if (s[i] != '0') break;
		if (check_float_digits(s, i, nidx, base ))
			{
			i_digits += (int)( nidx - i );
			if (s[nidx] != '.')  // No fraction part so collect the number
				{
				if (base == BASE_16)
					number = stringbase2number(s, i, nidx, base);
				else
					number = string2number(s, i, nidx); // Faster than decimal2number for large integer portion of the number
				}
			}
		}

	// Floating point representation
	if (s[nidx] == '.') // Any fraction ?
		{
		idx = nidx + 1;                    // Find start of fraction
		nidx = s.find_first_of("eE", idx); // Find end of fraction
		if (nidx == std::string::npos)
			nidx = s.length();

		if (idx < nidx)
			fpart = true;
		// Remove trailing zero digits
		for (i = nidx - 1; i >= idx; i--, nidx--) if (s[i] != '0') break;
		if (check_float_digits(s, idx, nidx, base))
			{
			f_digits += (int)(nidx - idx);
			// collect both the integer and fraction portion of the number
			if (i_digits != 0)
				{
				if (i_digits < f_digits)
					{
					size_t j;
					for (i = idx - 1, j=i_digits; j>0; --i, --j )
						s[i] = s[i - 1];
					}
				else
					{
					s.erase(idx - 1, 1); --nidx; --idx;
					}
				}
			if(base==BASE_16)
				number = stringbase2number(s, idx - i_digits, nidx, base );
			else
				number = string2number(s, idx-i_digits, nidx);
			}
		nidx = s.find_first_of("eE", nidx);
		}

	if (nidx != std::string::npos && (s[nidx] == 'e' || s[nidx] == 'E'))
		{// Parse the exponent . which is max sizeof of int. so use regular int operations
		idx = nidx + 1;
		nidx = s.length();
		sign_expo = CHAR_SIGN('+');;
		if (idx < nidx && (s[idx] == '+' || s[idx] == '-'))
			{
			sign_expo = CHAR_SIGN(s[idx]);
			idx++;
			if (idx == nidx)
				{
				throw float_precision::bad_float_syntax();
				}	// Sign but no number
			}
		else
			if (idx >= nidx)
				{
				throw float_precision::bad_float_syntax();
				}  // E but no number

		if (idx < nidx)
			epart = true;
		if (check_float_digits(s, idx, nidx, base))
			{
			// Collect exponent using base 
			//std::cout << "Exponent length=" << (nidx - idx) << std::endl;  // DEBUG HVE
			for (i = idx; i < nidx; i++)
				{
				expo_e10 *= base;
				expo_e10 += s[i]>'9' ? s[i]-'a' : s[i] - '0';
				}
			}
		if (sign_expo < 0)
			expo_e10 = -expo_e10;

		// for illegal number
		if (ipart == false && fpart == false )
			{
			throw float_precision::bad_float_syntax();
			}  // no number before the E or no number at all
		}

	// Build the float_precision number
	int_precision ip(number);
	float_precision fp(ip, p+ (size_t)log10(p), m), fppow(0,p+(size_t)log10(p),m);
	expo_e10 -= (eptype)f_digits; 
	//extern std::string formatInt(uintmax_t, int=10);
	//float_precision fp2(10); fp2 = fp; fp2.exponent(0); double fu = fp;   // DEBUG
	//std::cout << "atofp: control entry=" << log10(fu) + fp.exponent()*log10(2) << std::endl;   // DEBUG
	
	// Build the correct power adjustment if needed
 	if (expo_e10 != 0)
		{	
		// Correct for expo
		if (expo_e10 > 0)
			{
			ip = ipow(int_precision(base), int_precision(abs(expo_e10)));
			fppow = float_precision(ip, fppow.precision(), m);
			fp *= fppow;
			}
		else
			if (expo_e10 < 0)
				{
				expo_e10 = abs(expo_e10);
				ip = ipow(int_precision(base), int_precision(expo_e10));
				fppow = float_precision(ip, fppow.precision(), m);
				fp *= _float_precision_inverse(fppow);
				}
		} 
	fp.sign(sign);
	fp.precision(p);
	return fp;
	}


	//	@author Henrik Vestermark (hve@hvks.com)
	//	@date		7/Sep/2021
	//	@brief 		Convert a string decimal number into a float_precision number
	//	@return 	std::string - The decimal floating point string	
	//	@param		"str"		-	ascii string of floating point number to convert
	// @param		"p"			- The precision of the number
	// @param		"m"			- The round mode of the number
	//
	//	@todo 	
	//
	// Description:
	//   Convert ascii string into a float_precision numbers 
	//    The ascii float format is based on standard C notation
	//	
	//
	float_precision _float_precision_atofp2(const char *str, size_t p, enum round_mode m)
	{
		int sign, sign_expo;
		int expo_e10 = 0, f_digit = 0;
		std::string::size_type i, nidx, idx;
		std::string s(str);
		std::string::iterator pos;
		std::vector<iptype> number(1, 0), sum(1, 0);
		bool ipart = false, fpart = false, epart = false;

		number.reserve(s.size() + 16);
		s = remove_separators(s);
		idx = 0;
		sign = +1;
		// Parse leading sign if any
		pos = s.begin();
		if (*pos == '+' || *pos == '-')  // 
		{
			sign = CHAR_SIGN(*pos);
			pos++;
			idx = 1;
			if (pos == s.end())
				throw float_precision::bad_int_syntax();
		}

		// Determine any significant, fraction sign or exponent sign
		nidx = s.find_first_of(".eE", idx);
		if (nidx == std::string::npos) // Only digits (INTEGER) if any
		{
			int_precision ip(str);			// Construct integer to int_precision
			float_precision fp(ip, p, m);	// Construct float_precision from int_precision 
			return fp;
		} // End of Integer parsing

		  // Floating point number starts here
		  // Pick up significant beteen idx and nidx 
		if (nidx > idx) // Number of digits before the . sign or exponent Ee
		{
			ipart = true;
			// Strip leading zeros
			for (i = idx; i != nidx; i++) if (s[i] != '0') break;
			if (check_float_digits(s, i, nidx))
			{
				//decimal2number(number, s, i, nidx);
				number = string2number(s, i, nidx); // Faster than decimal2number for large integer portion of the number
			}
		}

		// Floating point representation
		if (s[nidx] == '.') // Any fraction ?
			{
			idx = nidx + 1;                      // Find start of fraction
			nidx = s.find_first_of("eE", idx); // Find end of fraction
			if (nidx == std::string::npos)
				nidx = s.length();

			if (idx < nidx)
				fpart = true;
			// Remove trailing zero digits
			for (i = nidx - 1; i >= idx; i--, nidx--) if (s[i] != '0') break;
			if (check_float_digits(s, idx, nidx))
				{
				f_digit += (int)(nidx - idx);

				//for (; idx + MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS <= nidx; idx += MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS)
				//	kilo2number(number, s, idx);  // NOT Working

				for (; idx + MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS <= nidx; idx += MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS)
					trunk2number(number, s, idx);

				if (idx != nidx)
					decimal2number(number, s, idx, nidx);
				}
			nidx = s.find_first_of("eE", idx);
			}

		if (nidx != std::string::npos && (s[nidx] == 'e' || s[nidx] == 'E'))
			{// Parse the exponent . which is max sizeof of int. so use regular int operations
			idx = nidx + 1;
			nidx = s.length();
			sign_expo = CHAR_SIGN('+');;
			if (idx < nidx && (s[idx] == '+' || s[idx] == '-'))
			{
				sign_expo = CHAR_SIGN(s[idx]);
				idx++;
				if (idx == nidx)
				{
					throw float_precision::bad_float_syntax();
				}	// Sign but no number
			}
			else
				if (idx >= nidx)
				{
					throw float_precision::bad_float_syntax();
				}  // E but no number

			if (idx < nidx)
				epart = true;
			if (check_float_digits(s, idx, nidx))
			{
				// Collect exponent using base 10
				for (i = idx; i < nidx; i++)
				{
					expo_e10 *= BASE_10;
					expo_e10 += s[i] - '0';
				}
			}
			if (sign_expo < 0)
				expo_e10 = -expo_e10;

			// for illegal number
			if (ipart == false && fpart == false)
			{
				throw float_precision::bad_float_syntax();
			}  // no number before the E or no number at all
		}

		// Build the float_precision number
		int_precision ip(number);
		float_precision fp(ip, p, m), fppow(0, p + 2, m);
		expo_e10 -= f_digit;

		// Build the correct power adjustment if needed
		if (expo_e10 != 0)
			{
			ip = ipow(int_precision(10), int_precision(abs(expo_e10)));
			fppow = float_precision(ip, p, m);
			// Correct for expo
			if (expo_e10 > 0)
				fp *= fppow;
			else
				if (expo_e10 < 0)
					fp *= _float_precision_inverse(fppow);
			}
		fp.sign(sign);
		return fp;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		10/Sep/2021
//	@brief 		Convert float_precision numbers into int_precision
//	@return		int_precision - The int_precision of the floating point
//	@param		"a"			- float_precision number to convert
//
//	@todo 	
//		Can posible be optimzed by 
//		by vector<iptype> x.insert(x.begin(), fp.x.begin(), fp.x.begin()+chunk+1)
//			ip >>= 64-within; or similar
//
// Description:
//   Convert float_precision numbers into int_precision
//		
//
int_precision _float_precision_fptoip(const float_precision *fp)
	{
	size_t i;
	int_precision ip(1);
	eptype expo = fp->exponent();
	if (expo < 0 || fp->iszero() == true )
		return int_precision(0);
	if (expo == 0) return int_precision( 1 * fp->sign() );
	// fecth expo bits from the fp number after the '.'
	size_t chunk = expo / Bitsfptype;
	int within = expo % Bitsfptype;
	for (i = 0; i < chunk; ++i)
		{
		ip <<= Bitsfptype;
		if( i+1 < fp->size() )
			ip += fp->index( i+1 );
		}
	if (within != 0)
		{
		ip <<= within;
		if (i + 1 < fp->size())
			ip += (fp->index(i + 1) >> (Bitsfptype - within));
		}

	ip.sign(fp->sign());  // Set sign
	return ip;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		10/Sep/2021
//	@brief 		Convert float_precision numbers into string (integer representation)
//	@return		std::string - The decimal floating point string	
//	@param		"a"			- float_precision number to convert
//
//	@todo 	
//
// Description:
//   Convert float_precision numbers into integer string (integer representation)
//
std::string _float_precision_fptoainteger(const float_precision *a)
	{
	int_precision ip = *a;
	return ip.toString();
	}

///////////////////////////////////////
//
// END CONVERT FLOAT PRECISION to and from ascii representation
//
///////////////////////////////////////

///////////////////////////////////////
//
// FLOATING POINT CORE FUNCTIONS. STRING
//
//   _float_precision_strip_trailing_zeros
//   _float_precision_rounding
//   _float_precision_uadd_short
//
//   Works Directly on the string class of the float number.
//		Left of the conversion to Binary. assest if this can be done differently
//
///////////////////////////////////////

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		Remove trailingnosignificant zeros from the string
//	@return 	nothing	
//	@param		"s"	-	digital string
//
//	@todo  
//
// Description:
//   Remove trailing nosignificant zeros
//
void _float_precision_strip_trailing_zeros( std::string *s )
	{
	std::string::reverse_iterator pos;
	int count;

	// Strip trailing zeros
	for( count = 0, pos = s->rbegin(); pos != s->rend() && FDIGIT( *pos ) == 0; pos++ )
         count++;
      
	s->erase( s->length() - count, count );
	if( s->length() == 0 )
		*s = FCHARACTER(0);

	return;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		Round the mantisaa to significant digits and rounding control
//	@return 	int - Return the exponent adjustment (0 or 1) 
//	@param		"m"	-	digital string
// @param		"sign"   - The sign of the number
// @param		"precision" - The digital precision
// @param		 "mode"   - Rounding mode 
//
//	@todo  
//
// Description:
//   Rounding control
//   Round the fraction to the number of precision based on the round mode 
//   Note that the mantissa number has ALWAYS been normalize prior to rounding
//   The mantissa NEVER contain a leading sign
//   Rounding Mode Positive numnber   Result    
//   Rounding to nearest              +�   
//   Rounding toward zero (Truncate)  Maximum, positive finite value   
//   Rounding up (toward +�)          +�   
//   Rounding down) (toward -�)       Maximum, positive finite value   
//
//   Rounding Mode Negative number    Result    
//   Rounding to nearest              -�   
//   Rounding toward zero (Truncate)  Maximum, negative finite value   
//   Rounding up (toward +�)          Maximum, negative finite value   
//   Rounding down) (toward -�)       -�   
//
int _float_precision_rounding( std::string *m, int sign, size_t precision, enum round_mode mode )
   {
   enum round_mode rm = mode;

   if( m->length() > precision )  // More digits than we need 
      {
      if( rm == ROUND_NEAR )
         {
         if( 2 * FDIGIT( (*m)[ precision ] ) >= BASE_10 )
            rm = ROUND_UP; //sign < 0 ? ROUND_DOWN : ROUND_UP;
         else
            rm = ROUND_DOWN; // sign < 0 ? ROUND_UP : ROUND_DOWN;
         }
      else
         if( rm == ROUND_UP && sign < 0 )
            rm = ROUND_DOWN;
         else
            if( rm == ROUND_DOWN && sign < 0 )
               rm = ROUND_UP;

      // Chuck excessive digits
      m->erase( (std::string::size_type)precision, m->length() - precision );

      if( rm == ROUND_UP ) 
         {
         size_t before;

         before = m->length();
         *m = _float_precision_uadd_short( m, 1 );
         if( m->length() > before )
            {
            if( m->length() > precision )
               m->erase( (std::string::size_type)precision, m->length() - precision );

            _float_precision_strip_trailing_zeros( m );            
            return 1;
            }
         }
      }

   _float_precision_strip_trailing_zeros( m );            
   return 0;
   }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		add a short integer to a floating point string (mantissa)
//	@return 	std::string - Return the added string
//	@param		"src1"	-	The source string
// @param		"d"  - The number to add
//
//	@todo  
//
// Description:
//   Short float Add: The digit d [0..F_RADIX] is added to the unsigned fraction string
//   Optimized 0 add or early out add is implemented
//
std::string _float_precision_uadd_short( std::string *src1, unsigned int d )
   {
   unsigned short ireg;
   std::string::reverse_iterator r1_pos, rd_pos;
   std::string des1;

   des1.reserve( src1->capacity() );
   des1 = *src1;
   if( d > BASE_10 )
      {
      throw float_precision::out_of_range();
      }

   if( d == 0 )   // Zero add
      return des1;

   ireg = (unsigned short)( BASE_10 * d );
   rd_pos = des1.rbegin();
   r1_pos = src1->rbegin();
   
   for(; r1_pos != src1->rend(); r1_pos++, rd_pos++ )
      {
      ireg = (unsigned short)( FDIGIT( *r1_pos ) + FCARRY( ireg ) ); 
      *rd_pos = FCHARACTER( (unsigned char)FSINGLE( ireg ) );
      if( FCARRY( ireg ) == 0 ) // Early out add
         break;
      }

   if( FCARRY( ireg ) != 0 )  // Insert the carry in the front of the number
      des1.insert( (std::string::size_type)0, 1, FCHARACTER( (unsigned char)FCARRY( ireg ) ) );

   return des1;
   }

//////////////////////////////////////
//
//	END of CORE Functions. STRING
//
//////////////////////////////////////


///////////////////////////////////////
//
// FLOATING POINT CORE BINARY
//
//   _float_precision_clz						-- Count trailing zeros starting in fptype
//   _float_precision_clz						-- Count trailing zeros starting in vector<fptype>
//	 _float_precision_ctz						-- Count trailing zeros starting in fptype
//	 _float_precision_ctz						-- Count trailing zeros starting in vector<fptype>
//	 _float_precision_strip_leading_zeros		-- Strip leading significant zeros
//   _float_precision_strip_trailing_zeros		-- Strip trailing zeros
//   _float_precision_normalize					-- Normalize the float number 
//   _float_precision_rounding					-- Round the number to given precision
//   _float_precision_right_shift				-- >> shift a vector<fptype> number
//   _float_precision_left_shift				-- <<  shift a vecotr<fptype> number
//   _float_precision_compare					-- Compare to vector<fptype> numbers and determine <,==,>
//   _float_precision_uadd_short				-- Add a vector<fptype> with a ftype number
//   _float_precision_uadd						-- Add two vector<fptype> numbers
//   _float_precision_usub_short				-- Subtract a ftype number from a vector<fptype>
//   _float_precision_usub						-- Subtract two vector<fptype> numbers
//   _float_precision_umul_short				-- Multiply a vector<fptype> with an ftype number
//   _float_precision_umul						-- Multiply two vector<fptype> numbers
//   _float_precision_umul_school				-- Multiply two vector<fptype> numbers using schoolbook algorithm
//   _float_precision_umul_linear				-- Multiply two vector<fptype> numbers using Schonhagen Strassen with linear convolution
//   _float_precision_umul_fourier				-- Multiply two vector<fptype> numbers using FFT
//   _float_precision_udiv_short				-- Divide vector<fptype> with a fptype number
//   _float_precision_udiv						-- Divide two vector<fptype> numbers
//   _float_precision_urem_short				-- Rem vector<fptype> with fptype number
//	 _float_precision_urem						-- Rem two vector<fptype> numbers 
//
//   Works Directly on the string class of the float number
//
///////////////////////////////////////

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		5/Sep/2021
//	@brief 		_float_precision_clz
//	@return		unsigned int	-	the count of leading zero bits in "a"
//	@param		"a"	-	fptype operand
//
//	@todo
//
// Description:
//   Count leading nosignificant zeros of the binary fptype 
//
size_t _float_precision_clz(const fptype a)
	{
	fptype x = a;
	size_t offset = 0;
	static const unsigned char lookup[16] = { 4,3,2,2,1,1,1,1,0,0,0,0,0,0 };

	if (sizeof(fptype) > 4 && x & 0xffffffff00000000u)
		x >>= 32; else offset += 32;

	if (sizeof(fptype) > 2 && x & 0xffff0000u)
		x >>= 16; else offset += 16;

	if (sizeof(fptype) > 1 && x & 0xff00u)
		x >>= 8; else offset += 8;

	if (x & 0xf0u)
		x >>= 4; else offset += 4;
	offset += lookup[(unsigned char)x];
	return offset;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Sep/2021
//	@brief 		_float_precision_clz
//	@return		unsigned int	-	the count of leading zero bits in "a"
//	@param		"mb"	-	vector<fptype> operand
//
//	@todo
//
// Description:
//   Count leading nosignificant zeros of the binary fptype 
//
size_t _float_precision_clz(const std::vector<fptype> &mb, size_t start)
	{
	size_t tot_cnt = 0, cnt;
	for (size_t i = start; i < mb.size(); ++i)
		{
		cnt = _float_precision_clz(mb[i]);
		tot_cnt += cnt;
		if (cnt != Bitsfptype ) 
			break;
		}
	return tot_cnt;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		5/Sep/2021
//	@brief 		_float_precision_ctz
//	@return		unsigned int	-	the count of trailing zero bits in "a"
//	@param		"a"	-	fptype operand
//
//	@todo
//
// Description:
//   Count trailing nosignificant zeros of the binary fptype 
//	  fptype bit word input to count zero bits on right
//   cnt will be the number of zero bits on the right,
//   so if a is 1101000 (base 2), then c will be 3
// NOTE: if 0 == a, then c = 64.
//
size_t _float_precision_ctz(const fptype a)
	{
	fptype x = a;  // sizeof(fptype) can be 8, 4, 2, or 1 
	size_t cnt;
	if (x & 0x1)
		cnt = 0;
	else
		{
		cnt = 1;
		if (sizeof(fptype) >= 8 && (x & 0xffffffffu) == 0)
			{
			x >>= 32; cnt += 32;
			}		
		if (sizeof(fptype) >= 4 && (x & 0xffffu) == 0)
			{
			x >>= 16; cnt += 16;
			}
		if (sizeof(fptype) >= 2 && (x & 0xff) == 0)
			{
			x >>= 8; cnt += 8;
			}
		if ((x & 0xf) == 0)
			{
			x >>= 4; cnt += 4;
			}
		if ((x & 0x3) == 0)
			{
			x >>= 2; cnt += 2;
			}
		cnt -= x & 0x1;
		}

	return cnt;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Sep/2021
//	@brief 		_flot_precision_ctz
//	@return		unsigned int	-	the count of trailing zero bits in "a"
//	@param		"mb"	-	vector<fptype> operand
//
//	@todo
//
// Description:
//   Count trailing nosignificant zeros of the binary fptype 
//	  fptype bit word input to count zero bits on right
//   cnt will be the number of zero bits on the right,
//   so if a is 1101000 (base 2), then c will be 3
// NOTE: if 0 == a, then c = 64.
//
size_t _float_precision_ctz(const std::vector<fptype> &mb)
	{
	size_t tot_cnt = 0, cnt;
	for (size_t i = 0; i < mb.size(); ++i)
		{
		cnt = _float_precision_ctz(mb[i]);
		tot_cnt += cnt;
		if (cnt != Bitsfptype ) break;
		}
	return tot_cnt;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		5/Sep/2021
//	@brief 		_float_precision_strip_leading_zeros
//	@return		void	-	
//	@param		"s"	-	pointer to source operand
//
//	@todo
//
// Description:
//   Remove leading nosignificant zeros of the binary number
//	This is from the start of the vector<fptype>
//
void _float_precision_strip_leading_zeros(std::vector<fptype> *s)
	{
	std::vector<fptype>::iterator pos;

	// Strip leading zeros
	for (pos = s->begin(); pos != s->end() && *pos == (fptype)0; ++pos);	// Find first not zero digit

	if (s->begin() != pos)
		s->erase(s->begin(), pos);
	if (s->empty())
		s->assign(1, (fptype)0);

	return;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/Sep/2021
//	@brief 		_float_precision_strip_trailing_zeros
//	@return		void	-	
//	@param		"s"	-	pointer to source operand
//
//	@todo
//
// Description:
//   Remove trailing nosignificant zeros of the binary number
//		this is from the top of the vector<fptype> 
//
void _float_precision_strip_trailing_zeros(std::vector<fptype> *s)
	{
	size_t i;
	std::vector<fptype>::reverse_iterator pos;

	// Strip leading zeros
	for (i = s->size() - 1, pos = s->rbegin(); i > 0 && *pos == (fptype)0; ++pos, --i);

	s->resize(i + 1);  // Keep at least one digit by default

	return;
	}
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		9/Sep/2021
//	@brief 		Right shift a string number 
//	@return 	the result of the shift	
//	@param		"src"	-	digital string
// @param		"shift" - Number of digital shifts
//
//	@todo  
//
// Description:
//   Right shift number x byniary digits by inserting 0 bits in front of the number.
//	  Doing it from index 0..src->size()
//
std::vector<fptype> _float_precision_right_shift(const std::vector<fptype> *src, const size_t shift)
	{
	size_t shiftwidth = Bitsfptype, adding, within;
	std::vector<fptype>::const_iterator pos, end;
	std::vector<fptype> des, c0(1, 0);
	fptype carry, mask;

	// Determine how many full digit shift and the last shift (last shift = shift count % Bitsfptype
	if (src->size() == 1 && src[0] == 0)  // Short cut: a zero number zero shifting right is still zero.
		return *src;
	if (shift == 0)  // Short cut: shift zero right does not change the number.
		return *src;

	adding = shift / Bitsfptype;
	des.insert(des.begin(), adding, 0);

	within = shift % Bitsfptype;
	shiftwidth -= within;
	mask = (~((fptype)0)) >> shiftwidth;  // Potential issue if shiftwidth is 64. needs to be tested
	for (carry = 0, pos = src->begin(), end=src->end(); pos != end; ++pos)
		{
		fptype n, nextcarry;
		n = *pos;
		if (within != 0)
			{
			nextcarry = n & mask;
			n >>= within;
			carry <<= shiftwidth;
			n |= carry;
			carry = nextcarry;
			}
		des.push_back(n);
		}
	if (carry != 0)
		des.push_back(carry << shiftwidth);
	if (des.size() == 0)
		des.push_back(0);
	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		9/Sep/2021
//	@brief 		Left shift a string number 
//	@return 	The result of the shift	
//	@param		"src"	-	digital string
// @param		"shift" - Number of digital binary shifts
//
//	@todo  
//
// Description:
//   Left shift number x bits. However bis that are shifted out of the first index ([0]) is discarded
//
std::vector<fptype> _float_precision_left_shift(const std::vector<fptype> *src, const size_t shift)
	{
	const unsigned int shiftwidth = Bitsfptype;
	std::vector<fptype>::const_reverse_iterator pos, end;
	std::vector<fptype> des;// , c0(1, 0);
	size_t discard, within;
	fptype carry, mask;

	// Determine how many full digit shift and the last shift (last shift = shift count % Bitsfptype
	if (src->size() == 1 && src[0] == 0)  // Short cut: a zero number zero shifting left is still zero.
		return *src;
	if (shift == 0)  // Short cut: shift zero left does not change the number.
		return *src;

	within = shift % Bitsfptype;
	mask = (~((fptype)0)) >> within;  mask = ~mask;
	discard	= shift / Bitsfptype;
	if (discard < src->size()) // Discard less  than size of src?
		{
		for (carry = 0, pos = src->rbegin(), end = src->rend(); pos != end - discard; ++pos)
			{
			fptype n, nextcarry;
			n = *pos;
			nextcarry = n & mask;
			n <<= within;
			n |= carry >> (shiftwidth - within);  // check is shiftwidth is 64 and witiin is zero and the result is still ok
			carry = nextcarry;
			des.push_back(n);
			}
		if (carry != 0&& (des.size()<src->size()-discard))
			des.push_back(carry >> (shiftwidth - within)); // check is shiftwidth is 64 and witiin is zero and the result is still ok
		reverse(des.begin(), des.end());
		}
	else
		des.push_back(0);  // Discarded more than size of src. Return 0
	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Sep/2021
//	@brief 		Normalize a floating point mantissa
//	@return 	int - Return the exponent adjustment factor due to normalization
//	@param		"m"	-	digital string
//
//	@todo  
//
// Description:
//   Normalize the mantissa
//   1) If a number does not have a leading digit != 0 then left shift until 
//   it has and adjust the exponent accordingly.
//   or 2) if a number does have more than one leading digits then right shift until it has ony one
//	  and adjust the exponent accordingly
//
eptype _float_precision_normalize(std::vector<fptype> *m)
	{
	eptype expo = 0;
	size_t shift;
	const size_t offset = Bitsfptype - 1;

	if (m->size() == 1 && *m == 0)  // m=0 is a special case
		return 0;

	shift = _float_precision_clz(*m);
	if (m->size()*Bitsfptype == shift) // All zeros (also a special case)
		{
		m->erase(m->begin()+1, m->end());  
		return 0;
		}
	if (shift != offset)
		{
		if (shift < offset)
			{
			*m = _float_precision_right_shift(m, offset - shift); 
			expo += (eptype)(offset - shift);
			}
		else
			{
			*m = _float_precision_left_shift(m, shift - offset);
			expo -= (eptype)( shift - offset );
			}
		}
	return expo;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		9/Sep/2021
//	@brief 		Round the mantisaa to significant digits and rounding control
//	@return 	int - Return the exponent adjustment (0 or 1) 
//	@param		"m"	-	digital string
// @param		"sign"   - The sign of the number
// @param		"precision" - The digital precision
// @param		"mode"   - Rounding mode 
//
//	@todo  
//
// Description:
//   Rounding control
//   Round the fraction to the number of precision based on the round mode 
//   Note that the fptype number has ALWAYS been normalize prior to rounding
//   Rounding Mode Positive numnber   Result    
//   Rounding to nearest              +�   
//   Rounding toward zero (Truncate)  Maximum, positive finite value   
//   Rounding up (toward +�)          +�   
//   Rounding down) (toward -�)       Maximum, positive finite value   
//
//   Rounding Mode Negative number    Result    
//   Rounding to nearest              -�   
//   Rounding toward zero (Truncate)  Maximum, negative finite value   
//   Rounding up (toward +�)          Maximum, negative finite value   
//   Rounding down) (toward -�) 
//		1) first check if we need to do any rounding at all
//		2) If mode == ROUND_NEAR determine if we are doing ROUND_DOWn or ROUND_UP
//		3) Discard excesive fptype digits that is beyond the precision
//		4) Discard excessive bits that is beyond precision
//		5) if Round up add one to the last bit of the precision. 
//		6) if 5) carry a digit to the most significant bit then >> shift 1 bit and return 1 for expoenent adjustment 
//		7) otherwise retun 0 for expoenent adjustment
//
int _float_precision_rounding(std::vector<fptype> *m, const int sign, const size_t precision, const enum round_mode mode)
	{
	enum round_mode rm = mode;
	const size_t pbits = (size_t)(ceil(precision*log2(BASE_10))); // Number of precision bits needed
	const size_t n =  pbits / Bitsfptype + 1;
	const size_t bn = pbits % Bitsfptype;
	size_t extra;

	if ((m->size() - 1)*Bitsfptype + 1 > pbits)  // More digits than we need 
		{
		//fptype check0; 
		size_t n1, bn1;
		switch (rm)
			{	
			case ROUND_NEAR:
				n1 = (pbits+1) / Bitsfptype + 1;
				bn1 = (pbits+1) % Bitsfptype;
				//check0 = (fptype)(1) << (Bitsfptype - bn1);  // can be removed after debug
				if (n1 < m->size() && (*m)[n1] & ((fptype)(1) << (Bitsfptype - bn1)))  // is last bit set -=> do ROUND_UP
					rm = ROUND_UP;
				else
					rm = ROUND_DOWN;
				break;
			case ROUND_UP:
				if (sign < 0)
					rm = ROUND_DOWN;
				break;
			case ROUND_DOWN:
				if (sign < 0)
					rm = ROUND_UP;
				break;
			case ROUND_ZERO:
				if (sign > 0)
					rm = ROUND_DOWN;
				else
					rm = ROUND_UP;
				break;
			}

		// Now rm is either ROUND_DOWN or ROUND_UP

		// Discard excessive fptype digits
		extra = bn == 0 ? 0 : 1;
		if(m->size()>n+extra)
			m->erase(m->begin()+n+extra,m->end());

		// Discard excessive bits within the last fptype.
		if (n < m->size())
			{
			//fptype check1 = (~(fptype)0) << (Bitsfptype - bn); // Debug can be removed after debug
			(*m)[n] &= ((~(fptype)0) << (Bitsfptype - bn));
			}

		if (rm == ROUND_UP)
			{
			const fptype bitvalue = (fptype)(1) << (Bitsfptype - bn );
			if (bitvalue == 0 )
				m->erase(m->begin() + n, m->end() );   // if we need to add to bit 0 in the previous fptype index. Erase the unneeded fptype digit
			*m = _float_precision_uadd_short( m, bitvalue );  // Add the carry
			if ((*m)[0] > (fptype)1 )
				{ // a carry was added to the most significant bit (now 2 instead of 1)
				// dont do comparison like: != 1 since it can also be zero when doing adding
				// Shift everything one to the right 
				*m = _float_precision_right_shift( m, 1 );
 				if (m->size() > n+1 )
					m->erase(m->begin()+n+1, m->end() );
				_float_precision_strip_trailing_zeros(m);
				return 1;
				}
			}
		}

	_float_precision_strip_trailing_zeros(m);
	return 0;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		9/Sep/2021, 1/Sep/2022
//	@brief 		Compare to floating point string (mantissa only)
//	@return 	int - Return the compared result. 0==same, 1==s1>s2 or -1==s1<s2
//	@param		"s1"	- First digital binary number
// @param		"s2"	- Second digital binary number
//
//	@todo  
//
// Description:
//   Compare two unsigned binary vector<fptype> numbers
//   and return 0 is equal, 1 if s1 > s2 otherwise -1
//	modified from pointers to reference
//
int _float_precision_compare(const std::vector<fptype>& s1, const std::vector<fptype>& s2)
	{
	std::vector<fptype>::const_iterator p1, p1_end, p2, p2_end;

	for (p1 = s1.begin(), p2 = s2.begin(), p1_end=s1.end(), p2_end=s2.end(); p1 != p1_end && p2 != p2_end; ++p1, ++p2)
		{
		if (*p1 > *p2 ) 
			return 1;	// s1 > s
		if (*p1 < *p2 ) 
			return -1;	// s1 < s2
		}
	// The are still the same and one or both is exhausted
	if (p1 == p1_end && p2 == p2_end)
		return 0; // is the same

	for (; p1 != p1_end; ++p1 )
		{
		if (*p1 > 0) return 1;	// s1 > s2
		}
	for (; p2 != p2_end; ++p2 )
		{
		if (*p2 > 0) return -1;  // s1 < s2
		}

	return 0;  // Same 
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		9/Sep/2021
//	@brief 		add a short integer to a floating point vector<fptype>
//	@return 	std::vector<fptype>- Return the added number
//	@param		"src"	-	The source string
// @param		"d"  - The number to add
//
//	@todo  
//
// Description:
//   Short float Add: The digit d [0..2^64] is added to the unsigned fraction string
//   Optimized 0 add or early out add is implemented
//
std::vector<fptype> _float_precision_uadd_short(const std::vector<fptype> *src, const fptype d)
	{
	fptype carry;
	std::vector<fptype>::const_reverse_iterator s_pos;
	std::vector<fptype>::reverse_iterator d_pos, d_end;
	std::vector<fptype> des;

	if (d == 0)   // Zero add
		return *src;

	carry = d;
	des = *src;		// Copy source to des1
	d_pos = des.rbegin();
	d_end = des.rend();
	s_pos = src->rbegin();

	for (; carry != 0 && d_pos != d_end; ++s_pos, ++d_pos)
		{
		*d_pos += carry;
		if (*d_pos < *s_pos)
			carry = 1;  // Set Carry
		else carry = 0;
		}

	// Exhaust the smalles of the number, so only the carry can changes the uppper radix digits
	for (; carry != 0 && d_pos != d_end; )
		{
		fptype tmp = *d_pos;
		*d_pos = tmp + carry;
		if (*d_pos < tmp)
			carry = 1;  // Set Carry
		else carry = 0;
		++d_pos;
		}

	// No more carry or end of upper radix number. 
	if (carry != 0) // If carry add the carry as a extra digit to the front of the number
		des.insert(des.begin(),1,carry);

	_int_precision_strip_trailing_zeros(&des);

	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Sep/2021
//	@brief 		std::vector<fptype> _int_precision_uadd
//	@return		std::vector<fptype>	-	the result of adding src1 and src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   Add two unsigned vector<fptype> numbers.
//   Notice src1.size() doesnt necessary needs to be of the same size as src2.size()
//	  src1[0] is the most significant elements (only containing 1bit of information), src1[src1.size()-1] the least significant
//	  same for src2.
//   Optimized: Used early out add
//
std::vector<fptype> _float_precision_uadd(const std::vector<fptype> *src1, const std::vector<fptype> *src2)
	{
	fptype carry = 0;
	std::vector<fptype> des;
	std::vector<fptype>::const_iterator pos;
	fptype tmp;
	size_t i;

	if (src1->size() >= src2->size())
		{
		des = *const_cast<std::vector<fptype> *>(src1);
		pos = src2->begin();
		i = src2->size();
		}
	else
		{
		des = *const_cast<std::vector<fptype> *>(src2);
		pos = src1->begin();
		i = src1->size();
		}

	for (; i > 0; --i )
		{ // Adding element by element for the two numbers starting with least significant vector<fptype> of the smallest number
		tmp = pos[i - 1] + carry;
		carry = tmp < pos[i - 1] ? 1 : 0;  // Next carry
		des[i-1] += tmp;
 		carry = des[i-1] < tmp ? 1 : carry;
		}

	//  The last (most significant elements) cant overflow since there is only 1 bit of significant for a normalized number
//	if (carry != 0)
//		carry = carry;  // Error

	_float_precision_strip_trailing_zeros(&des);
	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Sep/2021
//	@brief 		subtract a short integer from a floating point vector<fptype>
//	@return 	std::vector<fptype> - Return the subtracted vector<fptype> number
// @param		"result" -  If the number wraps around (d > src ) then result=1 otherwise 0
//	@param		"src"	-	The source string
// @param		"d"  - The number to subtract
//
//	@todo  
//
// Description:
//   Short Subtract: The fptype digit d [0..2^64] is subtracted from a unsigned vector<fptype> number
//   if src1 < src2 return -1 (wrap around) otherwise return 0 (no wrap around)
//
std::vector<fptype>_float_precision_usub_short(int *result, const std::vector<fptype> *src1, const fptype d)
	{
	fptype r, borrow=0;
	std::vector<fptype> ::const_reverse_iterator pos, end;
	std::vector<fptype> des;

	if (d == 0) // Nothing to subtract
		{
		*result = 0;
		return *src1;
		}

	pos = src1->rbegin();
	end = src1->rend();
	r = *pos - (d + borrow);
	borrow = *pos < (d + borrow) ? 1 : 0;
	des.push_back(r);
	for (++pos; borrow>0 && pos != end; ++pos)
		{
		r = *pos - borrow;
		borrow = *pos < borrow ? 1 : 0;
		des.push_back(r);
		}
	_float_precision_strip_trailing_zeros(&des);
	reverse(des.begin(), des.end());
	*result = borrow > 0 ? -1 : 0;
	return des;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		12/Sep/2021
//	@brief 		subtract two floating point string
//	@return 	std::vector<fptype> - Return the subtracted vector<fptype> number
// @param		"result" -  If the number wraps around (d > src1 ) then result=1 otherwise 0
//	@param		"src1"	-	The first source string
// @param		"src2"  - The second source string
//
//	@todo  
//
// Description:
//   Subtract two unsigned vector<fptype> numbers
//   src1 or src2 need not be of the same size
//   if src1 < src2 return -1 (wrap around) otherwise return 0 (no wrap around)
//
std::vector<fptype> _float_precision_usub(int *result, const std::vector<fptype> *src1, const std::vector<fptype> *src2)
	{
	fptype r, borrow = 0;
	std::vector<fptype>::const_reverse_iterator pos1, pos2;
	std::vector<fptype> des;
	size_t s1len=src1->size(), s2len=src2->size(), icur=std::max(s1len, s2len);

	if (s1len > s2len)
		des.reserve(src1->capacity());  // Reserver space to avoid time consuming reallocation
	else
		des.reserve(src2->capacity());  // Reserver space to avoid time consuming reallocation
	pos1 = src1->rbegin();
	pos2 = src2->rbegin();
	
	for ( ; icur > 0; --icur )
		{
		if (icur<=s1len && icur <= s2len )
			{
			r = *pos1 - (*pos2 + borrow);
			borrow = *pos1 < (*pos2 + borrow) ? 1 : 
					 *pos1 == 0 ? borrow : 0;      // if borrow was not paid then propagate it to next fptype subtraction
			++pos1; ++pos2;
			}
		else
			if (icur <= s1len )
				{
				r = *pos1 - (borrow);
				borrow = *pos1 < borrow ? 1 : 0;
				++pos1;
				}
			else
				{// icur <=s2len
				r = 0 - (*pos2 + borrow);
				borrow = 1;
				++pos2;
				}
		des.push_back(r);
		}
	reverse(des.begin(), des.end());
	_float_precision_strip_trailing_zeros(&des);
	*result = borrow>0 ? -1 : 0;
	return des;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		13/Sep/2021
//	@brief 		std::vector<fptype> _float_precision_umul_short
//	@return 	std::vector<fptype> - 	the result of the short multiplication
//	@param      "src1"	-	Source string to multiply short number
//	@param      "d"	   -	Number to multiply   
//
//	@todo
//
// Description:
//   Short Mul: The digit d [0..2^64] is multiplied to the unsigned vector<fptype> number
//   Optimized Multiply with zero yields zero, Multiply with one return the original.
//
std::vector<fptype> _float_precision_umul_short(const std::vector<fptype> *src1, const fptype d)
	{
	fptype carry = 0;
	std::vector<fptype>::const_reverse_iterator pos, end;
	std::vector<fptype> des;
	std::vector<fptype> tmp(2);

	if (d == 0)  // Multiply by zero is zero.
		{
		des.push_back(0);
		return des;
		}

	if (d == 1)  // Multiply by one dont change the src1.
		{
		des = *src1;
		_float_precision_strip_trailing_zeros(&des);
		return des;
		}

	des.reserve(src1->capacity());  // Reserver space to avoid time consuming reallocation   
	pos = src1->rbegin();
	end = src1->rend();

	for (; pos != end; ++pos)
		{
		tmp = _precision_umul64(*pos, d);
		tmp[0] += carry;
		carry = (tmp[0] < carry) ? tmp[1] + 1 : tmp[1];
		des.push_back(tmp[0]);
		}

	if (carry != 0)
		des.push_back(carry);
	reverse(des.begin(), des.end()); 
	_float_precision_strip_trailing_zeros(&des);

	return des;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		13/Sep/2021
//	@brief 		std::vector<fptype>  _float_precision_umul
//	@return		std::vector<fptype> -	the result of multiplying src1 and src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
// Multiply two unsigned vector<fptype>. Brach out to the different multiplication algorithm based on size of the operands for optimized performance
//
std::vector<fptype> _float_precision_umul(const std::vector<fptype> *src1, const std::vector<fptype> *src2)
	{
	const size_t s1len = src1->size(), s2len = src2->size();
	std::vector<fptype> des;

	if (*src1->begin() == 0 || *src2->begin() == 0)  // zero can only arise if the number is zero and therefore no need to check the size
		{
		des.assign(1, 0);
		}
	else
		{  // Both s1[0] and s2[0] == 1 check size to dermine if it is 1 or any power of 2
		if (s1len == 1)
			{
			des = *src2;  // Mutiply with 1 or any true power of 2
			}
		else
			if (s2len == 1)
				{
				des = *src1; // Mutiply with 1 or any true power of 2
				}
			else
				{// src1 and src2 size > 1. 
				if (s2len == 2 || s1len + s2len < 20)	// Measured for best performance .
					{
					des = _float_precision_umul_school(src1, src2);
					}
				else
					if (s1len+s2len < 6000)
						{
						des = _float_precision_umul_linear(src1, src2);
						}
					else
						{
						des = _float_precision_umul_fourier(src1, src2);
						}
				}
		}
	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		13/Sep/2021
//	@brief 		std::vector<fptype>  _float_precision_umul_school
//	@return		std::vector<fptype> -	the result of multiplying src1 and src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
// Multiply two unsigned vector<fptype>.
//	Not used anymore since the complexity is o(n^2)
//
std::vector<fptype> _float_precision_umul_school(const std::vector<fptype> *src1, const std::vector<fptype> *src2)
	{
	fptype carry;
	std::vector<fptype> des, tmp;
	std::vector<fptype>::const_reverse_iterator pos1, pos2;

	des.insert(des.begin(), src1->size() + src2->size(), 0);
	size_t i = src2->size(), j = src1->size(), k;
	tmp.assign(2, 0);
	for (pos2 = src2->rbegin(); i > 0; ++pos2, --i)
		{
		carry = 0; j = src1->size();
		for (pos1 = src1->rbegin(), k = i + j; j>0 && k >= i; ++pos1, --k, --j)
			{
			if (*pos2 == 1)
				{
				tmp[0] = *pos1; tmp[1] = 0;
				}
			else
				if (*pos1 == 1)
					{
					tmp[0] = *pos2; tmp[1] = 0;
					}
				else
					tmp = _precision_umul64(*pos1, *pos2);
			tmp[0] += carry;
			tmp[1] += tmp[0] < carry ? 1 : 0;
			des[k - 1] += tmp[0];
			carry = tmp[1];
			carry += des[k - 1] < tmp[0] ? 1 : 0;  //can't overflow by just adding 1 since tmp[1] will be sufficient less than maximum number
			}
		if (carry != 0)  // Note: no overflow for last carry since number can't be bigger than the sum of the two digits
			des[k - 1] = carry;
		}
	_float_precision_strip_leading_zeros(&des);
	_float_precision_strip_trailing_zeros(&des);
	return des;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		13/Sep/2021
//	@brief 		std::vector<fptype> _float_precision_umul_fourier
//	@return 	std::vector<fptyp> -	the result of multplying src1 and src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//		Multiply two unsigned binary numbers
//		Optimized: Used FFT algorithm to performed the multiplication
//		Plus added multi threading speeding up the calculation
//		Since we convert the binary numbers into float we have to ensure proper accuracy in the calculation.
//		In numerical recipies in C (2nd edition, chaper 20.6) they state that using double the equations that need to be fulfilled for accuracy
//		1byte binary:
//			log2(256^2)+log2(N)+safeguard*log2(log2(N))< digits in double which is 53. using safeguard as 1 we get for N=10^9 digits
//			16+30+4.9=50.9  which should be just Ok for 1 byte binary digits.
//		2byte binary:
//			log2(256^2^2)+log2(N)+safeguard*log2(log2(N))< digits in double which is 53. using safeguard as 1 we get for N=10^5 digits
//			32+16.6+4.05=52.65 Only 10^5 decimal digits is usualy not enough for arbitrary precsion so we are using a max of 1 byte.	
//
std::vector<fptype> _float_precision_umul_fourier(const std::vector<fptype> *src1, const std::vector<fptype> *src2, int nbits)
	{
	std::vector<fptype> des;
	size_t n, l, l1, l2, j;
	double cy;
	int bits = nbits==0?8:nbits;
	int radix = bits==8 ? 256: 16;
	size_t sz = sizeof(fptype);
	std::vector<double> va, vb;

	l1 = src1->size();
	l2 = src2->size();
	des.reserve(l1 + l2 + 16);  // Ensure enough space to hold the Multiplication result to avoid reallocation of des
	l = l1 < l2 ? l2 : l1;
	// Since we split the 64bit numbers into chunk of 8bit to ensure we have enough accuray when using double 
	l *= sizeof(fptype);  // Convert to byte size 
	if (l > 6'000'000*sizeof(fptype) || bits==4)  
		{
		bits = 4; l <<= 1; radix = 16; sz *= 2; // use 2^4 instead of 2^8
		//std::cout << "float_umul_fourier" << " do 4bit" << std::endl;  // DEBUG
		}
	for (n = 1; n < l; n <<= 1)  ;
	n <<= 1;

#ifdef HVE_THREAD
	// Using parallel sections below speeds up the performance of the two calls to _int_real_Fourier() with a factor of 1.8 
	if (nbits == 0|| l1+l2>10'000)
		{// Starting thread using lambda expressions
		// L1, l2, va, vb by reference since it is used after the thread has terminated
		std::thread first( [&, n, bits]() 
			{std::vector<fptype>::const_iterator pos, end;
			size_t i;
			va.resize(n);
			for (i = 0, pos = src1->begin(), end = src1->end(); pos != end; ++pos)
				i += convertbinary2double(&va[i], *pos, i == 0, bits);
			l1 = i; // L1 now Number of bytes or nibbles
			_vector_real_fourier(va, n, 1); // FFT va
			} );
		
		std::thread second([&, n, bits]() 
			{std::vector<fptype>::const_iterator pos, end;
			size_t i;
			vb.resize(n);
			for (i= 0, pos = src2->begin(), end = src2->end(); pos != end; ++pos)
				i += convertbinary2double(&vb[i], *pos, i == 0, bits);
			l2 = i; // L2 now Number of bytes or nibbles
			_vector_real_fourier(vb, n, 1); // FFT vb
			});

		first.join();
		second.join();
		}
	else
#endif
		{
		std::vector<fptype>::const_iterator pos, end;
		va.resize(n);
		vb.resize(n);
		// Start with most significant fptype e.g. src1[0]
		for (l1 = 0, pos = src1->begin(), end = src1->end(); pos != end; ++pos)
			l1 += convertbinary2double(&va[l1], *pos, l1 == 0, bits);
		// L1 now Number of bytes or nibbles
		// Start with most significant fptype e.g. src2[0]
		for (l2 = 0, pos = src2->begin(), end = src2->end(); pos != end; ++pos)
			l2 += convertbinary2double(&vb[l2], *pos, l2 == 0, bits);
		// L2 now number of bytes or nibbles
 		_vector_real_fourier(va, n, 1); // FFT va
		_vector_real_fourier(vb, n, 1); // FFT vb
		}

	// Do the multiplication in the frequence domain
	vb[0] *= va[0];
	vb[1] *= va[1];
	for (j = 2; j < n; j += 2)
		{
		double t;
		vb[j] = (t = vb[j])*va[j] - vb[j + 1] * va[j + 1];
		vb[j + 1] = t*va[j + 1] + vb[j + 1] * va[j];
		}

	_vector_real_fourier(vb, n, -1); // Reverse FFT vb
	for (cy = 0, j = 0; j <= n - 1; ++j)
		{
		double t;
		t = vb[n - 1 - j] / (n >> 1) + cy + 0.5;
		cy = (unsigned long)(t / radix);  // Byte Radix 2^8 or 2^4
		vb[n - 1 - j] = t - cy * radix;
		}

	// Now collect then back into a vector<fptype> format
	l1 += l2 - 1;				// max number of bytes or nibbles plus the carry
	l2 = l1 / sz;				// Number of full 64bit digits
	for (l = l2; l > 0; --l)	// do the full 64bit integers first starting backwards from b
		{
		fptype num;
		size_t inx = l1 - sz *(l2 - l + 1);
		num = convertdouble2binary(&vb[inx], sz, 0, bits);
		des.push_back(num);
		}
	l2 = l1 % sz;			// Number of remaing 8bits or 4bits digits
	if (l2>0 || cy != 0)	// do the the last 64bit integers from b
		{
		fptype num;
		num = convertdouble2binary(&vb[0], l2, cy,bits);
		des.push_back(num);
		}
	reverse(des.begin(), des.end());
	_float_precision_strip_trailing_zeros(&des);
	
	return des;
	}

	//	@author Henrik Vestermark (hve@hvks.com)
	//	@date		25/Dec/2021
	//	@brief 		std::vector<fptype> _float_precision_schonhage_strassen_linear_umul
	//	@return		std::vector<fptype>	-	the result of multiplying src1 and src2
	//	@param		"lhs"	-	First unsigned source argument
	//	@param		"rhs"	-	Second unsigned source argument
	//
	//	@todo
	//
	// Description:
	//   Multiply two unsigned decimal strings, using the Schonhage-Strassen (linear convolution) method
	//   Using the full advantages of Half bit size of fptype in the linear convolution
	//
	std::vector<fptype> _float_precision_umul_linear(const std::vector<fptype> *lhs, const std::vector<fptype> *rhs)
		{
		size_t i, j;
		size_t l_length = lhs->size(), r_length = rhs->size();
		const unsigned int HalfBitsfptype = Bitsfptype >> 1;  // same as / 2
		const fptype mask = (~(fptype)0) >> HalfBitsfptype;
		const uintmax_t radix = (uintmax_t)1 << HalfBitsfptype;
		std::vector<fptype> des;
		std::vector<fptype>::const_iterator pos, end;
		std::vector<uintmax_t> linearconvolution(2 * (l_length + r_length), 0);  // initialize it with zero
		std::vector<fptype> ua(l_length * 2), ub(r_length * 2);
		uintmax_t nextCarry = 0; 

		// Convert to half fptype from vector<fptype> and notice we dont stored in reverse order as we did for integers
		// by first converting lhs onto ua and then rhs into ub
		// e.g. lhs=a0+a1*R+a2*R^2,...an-1*R^n-1, a0 can be subdivied into half fptype  from fptype by mapping each mNumber number into 2 half fptype numbers 
		// the function convertbinary2Halfiptype() does this job per mNumber fptype number
		for (i = 0, pos = lhs->begin(), end = lhs->end(); pos != end; ++pos)
			i += convertbinary2Halfiptype(&ua[i], *pos, i == 0);
		l_length = i;  // l_length now in half iptype's  instead of iptype's
		for (j = 0, pos = rhs->begin(), end = rhs->end(); pos != end; ++pos)
			j += convertbinary2Halfiptype(&ub[j], *pos, j == 0);
		r_length = j;  // l_length now in half iptype's instead of iptype's

		// do the linear convolution
		for (i = 0; i < r_length; ++i)
			for (j = 0; j < l_length; ++j)
				{
				uintmax_t m = (uintmax_t)ub[r_length - 1 - i] * (uintmax_t)ua[l_length - 1 - j];
				linearconvolution[i + j] += m;
				if (linearconvolution[i + j] < m) // carry
					{
					// Propagate carry
#ifdef HVE_DEBUG
					//std::cout << "Propagating overflow at indx=" << i + j << std::endl;
#endif
					for (size_t k = i + j + 1; ; ++k)
						{
						linearconvolution[k] += radix;	// Add carry

						if (linearconvolution[k] >= radix) // Continue adding carry ?
							break;
#ifdef HVE_DEBUG
					//	else
					//		std::cout << "Propagating overflow at indx=" << k << " Radix=" << radix << " Val=" << linearconvolution[k] << std::endl;
#endif
						}
					}
				}

 		des.reserve(r_length + l_length + 2);
		for (i = 0; i < l_length + r_length - 1; ++i)
			{	
			linearconvolution[i] += nextCarry;
#ifdef HVE_DEBUG
	//		if (linearconvolution[i] < nextCarry)
	//			std::cout << "Final overflow at index=" << i  << " NextCarry=" << nextCarry << " Index value="<< linearconvolution[i] << std::endl;
#endif
			
			nextCarry = linearconvolution[i] >> HalfBitsfptype;  //  same as / radix;
			linearconvolution[i] &= mask;  // same as %= radix;
			}
		if (nextCarry != 0)
			linearconvolution[i++] = nextCarry & mask; // same as % radix;

		//linearconvolution now holds the result with [0] as the most significant byte number and [i-1] as the least sinificant byte as Halfbitsiptype numbers
		// Now convert then back into a vector<iptype> format. i is the number of HalfBitsiptype's in the result
		// do the full 64bit integers first starting from least significant HalfBitsiptype
		for (j = 0; j < i; j += 2)
			{
			fptype num;
			num = convertHalfiptype2binary(&linearconvolution[j], 2);
			des.push_back(num);
			}
	reverse(des.begin(), des.end());
	_float_precision_strip_trailing_zeros(&des);
	return des;
	}

///////////////////////////////////////////////////////////////////////////////////////
//
//
//	Specialized function for squaring using fourier and linear
//
//
///////////////////////////////////////////////////////////////////////////////////////

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		23/Apr/2022
//	@brief 		std::vector<fptype> _float_precision_umul2_fourier
//	@return 	std::vector<fptyp> -	the result of multplying src1 and src2
//	@param		"src1"	-	unsigned source argument
//
//	@todo
//
// Description:
//		squaring the unsigned binary number
//		Optimized: Used FFT algorithm to performed the multiplication
//		Plus added multi threading speeding up the calculation
//		Since we convert the binary numbers into float we have to ensure proper accuracy in the calculation.
//		In numerical recipies in C (2nd edition, chaper 20.6) they state that using double the equations that need to be fulfilled for accuracy
//		1byte binary:
//			log2(256^2)+log2(N)+safeguard*log2(log2(N))< digits in double which is 53. using safeguard as 1 we get for N=10^9 digits
//			16+30+4.9=50.9  which should be just Ok for 1 byte binary digits.
//		2byte binary:
//			log2(256^2^2)+log2(N)+safeguard*log2(log2(N))< digits in double which is 53. using safeguard as 1 we get for N=10^5 digits
//			32+16.6+4.05=52.65 Only 10^5 decimal digits is usualy not enough for arbitrary precsion so we are using a max of 1 byte.	
//
std::vector<fptype> _float_precision_umulsq_fourier(const std::vector<fptype> *src, int nbits)
	{
	std::vector<fptype> des;
	size_t n, l, l1, l2, j;
	double cy;
	int bits = nbits == 0 ? 8 : nbits;
	int radix = bits == 8 ? 256 : 16;
	size_t sz = sizeof(fptype);
	std::vector<double> va;

	l1 = src->size();
	des.reserve(l1 + l1 + 16);  // Ensure enough space to hold the Multiplication result to avoid reallocation of des
	l = l1;
	// Since we split the 64bit numbers into chunk of 8bit to ensure we have enough accuray when using double 
	l *= sizeof(fptype);  // Convert to byte size 
	if (l > 6'000'000 * sizeof(fptype) || bits == 4)
		{
		bits = 4; l <<= 1; radix = 16; sz *= 2; // use 2^4 instead of 2^8
		}
	for (n = 1; n < l; n <<= 1);
	n <<= 1;

	std::vector<fptype>::const_iterator pos, end;
	va.resize(n);
	// Start with most significant fptype e.g. src1[0]
	for (l1 = 0, pos = src->begin(), end = src->end(); pos != end; ++pos)
		l1 += convertbinary2double(&va[l1], *pos, l1 == 0, bits);
	// L1 now Number of bytes or nibbles
	_vector_real_fourier(va, n, 1); // FFT va

	// Do the multiplication in the frequence domain
	va[0] *= va[0];
	va[1] *= va[1];
	for (j = 2; j < n; j += 2)
		{
		double t = va[j], u = va[j + 1];
		va[j] = t*t - u*u;
		va[j + 1] = 2*t*u;
		}

	_vector_real_fourier(va, n, -1); // Reverse FFT vb
	for (cy = 0, j = 0; j <= n - 1; ++j)
		{
		double t;
		t = va[n - 1 - j] / (n >> 1) + cy + 0.5;
		cy = (unsigned long)(t / radix);  // Byte Radix 2^8 or 2^4
		va[n - 1 - j] = t - cy * radix;
		}

	// Now collect then back into a vector<fptype> format
	l1 += l1 - 1;				// max number of bytes or nibbles plus the carry
	l2 = l1 / sz;				// Number of full 64bit digits
	for (l = l2; l > 0; --l)	// do the full 64bit integers first starting backwards from b
		{
		fptype num;
		size_t inx = l1 - sz *(l2 - l + 1);
		num = convertdouble2binary(&va[inx], sz, 0, bits);
		des.push_back(num);
		}
	l2 = l1 % sz;			// Number of remaing 8bits or 4bits digits
	if (l2>0 || cy != 0)	// do the the last 64bit integers from b
		{
		fptype num;
		num = convertdouble2binary(&va[0], l2, cy, bits);
		des.push_back(num);
		}
	reverse(des.begin(), des.end());
	_float_precision_strip_trailing_zeros(&des);

	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		24/Jan/2022
//	@brief 		std::vector<fptype> _float_precision_schonhage_strassen_linear_umul
//	@return		std::vector<fptype>	-	the result of multiplying src1 and src2
//	@param		"lhs"	-	unsigned source argument
//
//	@todo
//
// Description:
//   square the the unsigned decimal string, using the Schonhage-Strassen (linear convolution) method
//   Using the full advantages of Half bit size of fptype in the linear convolution
//
std::vector<fptype> _float_precision_umulsq_linear(const std::vector<fptype> *src )
	{
	size_t i, j;
	size_t l_length = src->size();
	const unsigned int HalfBitsfptype = Bitsfptype >> 1;  // same as / 2
	const fptype mask = (~(fptype)0) >> HalfBitsfptype;
	const uintmax_t radix = (uintmax_t)1 << HalfBitsfptype;
	std::vector<fptype> des;
	std::vector<fptype>::const_iterator pos, end;
	std::vector<uintmax_t> linearconvolution(2 * (l_length + l_length), 0);  // initialize it with zero
	std::vector<fptype> ua(l_length * 2);
	uintmax_t nextCarry = 0;

	// Convert to half fptype from vector<fptype> and notice we dont stored in reverse order as we did for integers
	// by first converting lhs onto ua and then rhs into ub
	// e.g. lhs=a0+a1*R+a2*R^2,...an-1*R^n-1, a0 can be subdivied into half fptype  from fptype by mapping each mNumber number into 2 half fptype numbers 
	// the function convertbinary2Halfiptype() does this job per mNumber fptype number
	for (i = 0, pos = src->begin(), end = src->end(); pos != end; ++pos)
		i += convertbinary2Halfiptype(&ua[i], *pos, i == 0);
	l_length = i;  // l_length now in half iptype's  instead of iptype's

	// do the linear convolution
	for (i = 0; i < l_length; ++i)
		for (j = 0; j < l_length; ++j)
			{
			uintmax_t m = (uintmax_t)ua[l_length - 1 - i] * (uintmax_t)ua[l_length - 1 - j];
			linearconvolution[i + j] += m;
			if (linearconvolution[i + j] < m) // carry
				{
				// Propagate carry
				for (size_t k = i + j + 1; ; ++k)
					{
					linearconvolution[k] += radix;	// Add carry
					if (linearconvolution[k] >= radix) // Continue adding carry ?
						break;
					}
				}
			}

	des.reserve(l_length + l_length + 2);
	for (i = 0; i < l_length + l_length - 1; ++i)
		{
		linearconvolution[i] += nextCarry;
		nextCarry = linearconvolution[i] >> HalfBitsfptype;  //  same as / radix;
		//	if (nextCarry != 0)
		//		nextCarry = nextCarry;
		linearconvolution[i] &= mask;  // same as %= radix;
		}
	if (nextCarry != 0)
		linearconvolution[i++] = nextCarry & mask; // same as % radix;

	//linearconvolution now holds the result with [0] as the most significant byte number and [i-1] as the least sinificant byte as Halfbitsiptype numbers
	// Now convert then back into a vector<iptype> format. i is the number of HalfBitsiptype's in the result
	// do the full 64bit integers first starting from least significant HalfBitsiptype
	for (j = 0; j < i; j += 2)
		{
		fptype num;
		num = convertHalfiptype2binary(&linearconvolution[j], 2);
		des.push_back(num);
		}
	reverse(des.begin(), des.end());
	_float_precision_strip_trailing_zeros(&des);
	return des;
	}

///////////////////////////////////////////////////////////////////////////////////////
//
//
//	END squaring functions
//
//
///////////////////////////////////////////////////////////////////////////////////////

	 
// Short Division: The fptype digit d  is divide up into the unsigned fptype vector
//
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		13/Sep/2021
//	@brief 		std::vector<fptype>	_float_precision_udiv_short
//	@return 	std::vector<fptype>	- The result of the short division
//	@param      "src1"				- Source string to divide with the short number
//	@param      "d"					- Number to divide
// @param		"remaind"			- The remaind of the short division
//
//	@todo
//
// Description:
//   Short divide: The fptype digit d [0..2^32] is divided up in the unsigned vector<fptype> 
//	  Notice only up to int 32bit can be handle as short div.
//   Divide with zero throw an exception
//
std::vector<fptype> _float_precision_udiv_short(fptype *remaind, const std::vector<fptype> *src1, const fptype d)
	{
	const unsigned int shifts = 4 * sizeof(fptype);
	const fptype mask = (~((fptype)0)) >> shifts;
	fptype ir;
	std::vector<fptype>::const_iterator s_pos, s_end;
	std::vector<fptype> des;

	if (d == 0)
		{
		throw float_precision::divide_by_zero();
		}

	if (d == 1)  // Divide by one dont change the src1.
		{
		des = *const_cast<std::vector<fptype> *> (src1);
		_float_precision_strip_trailing_zeros(&des);
		*remaind = 0;
		return des;
		}

	des.reserve(src1->capacity());  // Reserver space to avoid time consuming reallocation
	s_pos = src1->begin();
	s_end = src1->end();
	for (ir = 0; s_pos != s_end; ++s_pos)
		{
		fptype n, qh, ql;
		/*if (ir == 0)
		{// Shortcut when carry is zero
		des.push_back(*s_pos / d );
		ir = *s_pos % d;
		}
		else*/
		{
			n = *s_pos >> shifts;
			n |= ir << shifts;
			qh = n / d;	ir = n % d;
			n = *s_pos & mask;
			n |= ir << shifts;
			ql = n / d;	ir = n % d;
			n = (qh << shifts) | ql;
			des.push_back(n);
		}
		}

	//reverse(des.begin(), des.end());
	_float_precision_strip_trailing_zeros(&des);
	*remaind = ir;
	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		13/Sep/2021
//	@brief 		std::vector<fptype> _float_precision_udiv
//	@return		std::vector<fptype>-	the result of disivison
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   Divide two unsigned binary numbers
//   Optimized: Used early out add and multiplication w. zero
//
std::vector<fptype> _float_precision_udiv(const std::vector<fptype> *src1, const std::vector<fptype> *src2)
	{
	int plusdigit, plusbits, wrap, i;
	fptype underflow;
	std::vector<fptype> des, tmp, quotient, divisor;

	des.push_back(0);
	divisor = *const_cast<std::vector<fptype> *> (src1);
	if (src2->size() == 1 && (src2->front() >> 32) == 0) // Make short div if denominator <= 32 bit integer.
		return _float_precision_udiv_short(&underflow, &divisor, src2->front());

	plusdigit = (int)divisor.size() - (int)src2->size();
	if (plusdigit < 0)  //src1 / src2 == 0
		return des;

	plusbits = (int)_float_precision_clz(src2->back()) - (int)_float_precision_clz(divisor.back());
	plusbits = plusdigit * Bitsfptype + plusbits;
	for (i = 0; plusbits >= 1; ++i)
		{
		tmp = _int_precision_ushiftleft(src2, plusbits);
		if (_int_precision_compare(&divisor, &tmp) < 0)
			{ // Too much reduce with one power of radix
			--plusbits; continue;
			}
		divisor = _int_precision_usub(&wrap, &divisor, &tmp);
		quotient.clear();
		quotient.insert(quotient.begin(), (plusbits / (Bitsfptype)) + 1, 0);
		quotient[quotient.size() - 1] = (fptype)(1) << (plusbits % (Bitsfptype));
		des = _int_precision_uadd(&des, &quotient);
		}

	for (wrap = 0; wrap == 0; )
		{
		divisor = _int_precision_usub(&wrap, &divisor, src2);
		if (wrap == 0) // src1 was indeed > src2
			des = _int_precision_uadd_short(&des, 1);
		}

	_int_precision_strip_trailing_zeros(&des);
	return des;
	}

// Short Remainder: The fptype digit d [1..2^64] is divide up into the unsigned vector<fptype> and the remaing is returned
//
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		13/Sep/2021
//	@brief 		std::vector<fptype> _float_precision_urem_short
//	@return 	std::vector<fptype> - 	the result of the short remainder
//	@param      "src1"	-	Source string to divide with the short number
//	@param      "d"	   -	Number to divide
//
//	@todo
//
// Description:
//   Short remainder: The fptype digit d [0..2^64] is divided up in the unsigned vector<fptype> and the remaing is retuened
//   Divide with zero throw an exception
//   if d==1 then result == 0, for d==2,4,5,8,10 we only test the last few digits to get the result. This speed up rem for large integers with small rem value
//   since we dont have to run through every digits in src1
//
std::vector<fptype> _float_precision_urem_short(const std::vector<fptype> *src1, const fptype d)
	{
	const unsigned int shifts = 4 * sizeof(fptype);
	const fptype mask = (~((fptype)0)) >> shifts;
	fptype ir;
	std::vector<fptype>::const_reverse_iterator s_pos, s_end;
	std::vector<fptype> des;

	if (d == 0)
		{
		throw float_precision::divide_by_zero();
		}

	if (d == 1)  // Remainer is always 0 for d==1
		{
		des.push_back(0);
		return des;
		}

	// Short cut
	ir = *src1->begin();
	switch (d)
		{
		case 2: des.push_back(ir % 2); return des; break;
		case 4: des.push_back(ir % 4); return des; break;
		case 5: des.push_back(ir % 5); return des; break;
		case 8: des.push_back(ir % 8); return des; break;
		case 10: des.push_back(ir % 10); return des; break;
		default:;   // No Short cut
		}

	s_pos = src1->rbegin();
	s_end = src1->rend();
	for (ir = 0; s_pos != s_end; ++s_pos)
		{
		fptype n;
		if (ir == 0)
			{// Shortcut when carry is zero
			ir = *s_pos % d;
			}
		else
			{
			n = *s_pos >> shifts;
			n += ir << shifts;
			ir = n % d;
			n = *s_pos & mask;
			n += ir << shifts;
			ir = n % d;
			}
		}

	des.push_back(ir);
	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		13/Sep/2021
//	@brief 		std::vector<fptype> _float_precision_urem
//	@return		std::vector<fptype>	-	the remaing result of divide src1 with src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   Find the remainder when divide two unsigned vector<fptype> numbers
//   Optimized: Used early out add and multiplication w. zero
//
std::vector<fptype> _float_precision_urem(const std::vector<fptype> *src1, const std::vector<fptype> *src2)
	{
	int wrap;
	std::vector<fptype> des, tmp;

	des.push_back(0);
	if (src2->size() == 1 && (src2->front() >> 32) == 0) // Make short rem 
		{
		fptype rem;
		_float_precision_udiv_short(&rem, src1, src2->front());
		des[0] = rem;
		return des;
		}

	tmp = _float_precision_udiv(src1, src2);
	tmp = _float_precision_umul(&tmp, src2);
	des = _float_precision_usub(&wrap, src1, &tmp);

	_float_precision_strip_trailing_zeros(&des);

	return des;
	}

//////////////////////////////////////
//
//	END of CORE Functions. BINARY
//
//////////////////////////////////////




///////////////////////////////////////
//
// END FLOATING POINT CORE FUNCTIONS
//
//
///////////////////////////////////////


///////////////////////////////////////
//
// FLOAT PRECISION FUNCTIONS
//
///////////////////////////////////////


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		4/Oct/2021
//	@brief 		Calculate the inverse of a 
//	@return 	float_precision -	Return 1/a
//	@param      "a"	-	The float_precision number to inverse
//
//	@todo  
//
// Description:
//   Inverse of V
//  
//	Optimized for inverting a true power of 2.
//	Using Brent suggestion with small performance enchancement 
//
//	Hybrid 3rd order method for calculation of the inverse of a float_precision number
//	with improve stopping criteria to avoid unnecessary calculation
//	It determine if the 3rd order or 2nd order method is the most efficient way to calculate the inverse
//
float_precision _float_precision_inverse(const float_precision& a)
	{
	const size_t extra = 5;
	const size_t precision = a.precision();
	const eptype expo = a.exponent();
	const intmax_t limit = -(intmax_t)((precision + 1)*log2(10)) - 1;
	const float_precision c1(1), c2(2);
	size_t digits, loopcnt = 1;
	double fx = precision + extra < 16 ? 1.0 : log(precision + extra) - log(16);
	bool halley = ceil(fx / log(2)) / ceil(fx / log(3)) > 1.58 ? true : false;
	const int step = halley == true ? 3 : 2;
	float_precision r, s, x, y(a);

	if (a.iszero() == true)
		{
		throw float_precision::divide_by_zero();
		}
	// if a is a true power of 2 then we dont need to iterate but just reverse the exponenent and return
	if (a.size() == 1)
		{
		y.exponent(-y.exponent());
		return y;
		}

	// find the inverse of y without exponent and adjust for exponent later
	y.exponent(0);	// y is in the interval [1..2[
					// Do iteration using extra digits higher precision
	x.precision(precision + extra);
	// Get a initial guess using ordinary floating point
	fx = 1 / (double)y;
	x = float_precision(fx);
	digits = halley == true ? 48 : 32;
	// cout << " Iteration type: " << (halley == true ? "Halley" : "Newton") << endl;  // DEBUG
	//double tloop=clock();										// DEBUG
	// Now iterate using 3rd order x=x+x(1-xy)+x(1-yx)^2
	for (digits = std::min(digits, precision); ; digits = std::min(precision + extra, digits * step), ++loopcnt)
		{
		//tloop = clock();  // DEBUG
		// Increase precision by a factor of two for the working variable s & x.
		s.precision(digits);
		x.precision(digits);
		if (halley == true)
			{// Use Halley 3rd order iteration
			// Only half the precision for r as suggested by Richard P. Brent
			r.precision(digits / 2 + 1);
			s = y;		// Build r
			s *= x;
			r = c1;
			r -= s;		// r=(1-yx)
			s = x;		// Build s
			s *= r;		// s=x(1-yx)
			x += s;		// x=x+x(1-yx) 
			s *= r;		// s=x(1-yx)^2
			x += s;		// x=x+x(1-yx)+x(1-yx)^2

			if (digits == precision + extra)
				{// Reach maximum precision
				fx = (double)limit / r.exponent();
				//cout << "LoopH=" << loopcnt << " Digits=" << digits << ((digits == precision + extra) ? 'T' : 'F') << " Limit=" << formatInt(limit) << " r expo=" << formatInt(r.exponent()) << " ratio=" << fx << " s expo=" << formatInt(s.exponent()) << endl;	  // DEBUG
				if (r.iszero() || fx < 2.5)
					{
					//	cout << " Fast break=" << formatInt(2*r.exponent()) << endl;
					break;
					}
				if (fx <= 4.0)
					halley = false;	// Close enough to use 2nd order in the last iteration
				}
			}
		else
			{// Use Newton 2nd order iteration
		 // Do a regular Newton step requires less work than a full 3rd order step
			s = y;					// y
			s *= x;					// xy
			s = c2 - s;             // 2-xy
			x *= s;					// x=x(2-xy)
			if (digits == precision + extra)
				{// Reach Maximum precision
				s.precision(precision + 1);	// round to final precision
				std::vector<fptype> *mp = s.pointer();
				size_t offset = _float_precision_clz(*mp, 1) + 1 + s.exponent();
				fx = -(double)limit / offset;
				//cout << "LoopN=" << loopcnt << " Digits=" << digits << ((digits == precision + extra) ? 'T' : 'F') << " Limit=" << formatInt(limit) << " r expo=-" << formatInt(offset) << " ratio=" << fx << endl;  // DEBUG
				if (s == c1 || fx < 1.9)	// break if no improvement
					break;
				s.precision(precision + extra);
				}
			}
		//	tloop = clock() - tloop;	// DEBUG
		//	cout << "\tLoop=" << loopcnt << " Digits="<< formatInt(digits) << " time=" << tloop << " msec" <<  " r expo="<<formatInt(r.exponent()) << endl;	// DEBUG
		}
	//	tloop = clock() - tloop;	// DEBUG
	//	cout << "\tFinal Loop=" << loopcnt <<" Digits=" << formatInt(digits) << " time=" << tloop << " msec" << endl;	// DEBUG

	// Reapply exponent, mode, and precision
	x.adjustExponent(-expo);
	x.mode(a.mode());
	x.precision(precision + 1);
	return x;
	}

/*  OLD Inverse method
float_precision _float_precision_inverse(const float_precision& a)
	{
	const size_t extra = 5;
	const size_t precision = a.precision();
	const eptype expo = a.exponent();
	const intmax_t limit = -(intmax_t)((precision + 1)*log2(10)) - 1;
	const float_precision c1(1);
	size_t digits, loopcnt = 1;
	double fx;
	float_precision r, s, x, y(a);

	if (a.iszero() == true)
		{
		throw float_precision::divide_by_zero();
		}
	// if a is a true power of 2 then we dont need to iterate but just reverse the exponenent and return
	if (a.size() == 1)
		{
		y.exponent(-y.exponent());
		return y;
		}

	// find the inverse of y without exponent and adjust for exponent later
	y.exponent(0);	// y is in the interval [1..2[
	// Do iteration using extra digits higher precision
	x.precision(precision + extra);
	// Get a initial guess using ordinary floating point
	fx = 1 / (double)y;
	x = float_precision(fx);

	// Now iterate using 3rd order x=x+x(1-xy)+x(1-yx)^2
	for (digits = std::min((size_t)48, precision); ; digits = std::min(precision + extra, digits * 3), ++loopcnt)
	{
		// Increase precision by a factor of two for the working variable s & x.
		s.precision(digits);
		x.precision(digits);
		// Only half the precision for r as suggested by Richard P. Brent
		r.precision(digits / 2 + 1);

		s = y;	s *= x;	r = c1;	r -= s;		//r = c1 - y * x;		// (1-yx)
		s = x; 	s *= r;		//s = x * r;			// x(1-yx)
		x += s;				// x = x + x(1 - yx) 
		if (2 * r.exponent() > limit)
			{ s *= r;	x += s; }	// x=x+x(1-yx)+x(1-yx)^2
		if (digits == precision + extra && (r.iszero() || 2 * r.exponent() < limit))
			break;
		}

	// Reapply exponent, mode, and precision
	x.adjustExponent(-expo);
	x.mode(a.mode());
	x.precision(precision + 1);
	return x;
	}
	*/

// Float Precision support functions

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/Oct/2021
//	@brief 		Calculate sqrt(x)
//	@return 	float_precision -	Return sqrt(a)
//	@param      "x"	-	The sqrt argument
//
//	@todo  
//
// Description:
//   sqrt(V)
//   Equivalent with the same standard C function call
//   Separate exponent. e.g. sqrt(V*2^x)=2^x/2*sqrt(V)
// his is a hybrid version that use Halley (3rd order) or Newton (2nd order method)
// It determine which of the two methods is the most usefull to use
//	The function has been improved using iterative deepening creating
//	a speed up with a factor of 3 over the classic method.
//  New hybrid sqrt mthod using Both Halleys 3rd order and Newton nd order convergence rate with dynamic precision
//
float_precision sqrt(const float_precision& a)
	{
	const unsigned int extra = 2;
	const size_t precision = a.precision();
	const eptype expo = a.exponent();
	const intmax_t limit = -(intmax_t)((precision + 1)*log2(10)) - 1;
	const float_precision c1(1), c2(2), c15(15), c3(3), c10(10);
	eptype expo_sq;
	size_t digits, loopcnt;
	double fx = precision + extra < 16 ? 1.0 : log(precision + extra) - log(16);
	bool halley = ceil(fx / log(2)) / ceil(fx / log(3)) > 1.58 ? true : false;
	const int step = halley == true ? 3 : 2;
	float_precision r, x, y(a), z;

	if (a.iszero() || a == c1)  // Simple square root
		return a;
	if (a.sign() < 0)
		{
		throw float_precision::domain_error();
		}

	if (y.size() == 1 && (expo & 0x1) == 0)
		{// True power of 2 and the expoenent is even
		y.exponent(y.exponent() >> 1); // Half the exponent and return the result.
		return y;
		}
	expo_sq = expo / 2;
	y.exponent(expo - 2 * expo_sq);
	// Do iteration using 2 digits higher precision
	y.precision(precision + extra);
	r.precision(precision + extra);
	x.precision(precision + extra);
	// Get a initial guess using ordinary floating point
	fx = (double)y;				// Convert to double	
	// set the initial guess with at approx 16 correct digits
	fx = 1 / sqrt(fx);
	x = float_precision(fx);
	digits = halley ? 48 : 32;
	//cout << "Iterations type=" << (halley == true ? "Halley" : "Newton") << endl;	// DEBUG

	// Now iterate using Halley 3rd order or Newton 2nd order iteration
	// Notice y is the original number to squareroot which has the full precision
	for (loopcnt = 1, digits = std::min(digits, precision); ; digits = std::min(precision + extra, digits * step), ++loopcnt)
		{
		r.precision(digits);
		x.precision(digits);
		if (halley == true)
			{// Use Halley 3rd order iteration
			z.precision(digits);
			r = y;						// y 
			r *= x.square();			// yx^2
			z = c3; z *= r; z = c10 - z;
			r *= z;						//r *= c10 - c3*r or r=yx^2(10-3yx^2)
			r = c15 - r;				// 15-yx^2*(10-3*yx^2)
			r.adjustExponent(-3);		// r=r/8   	
			x *= r;						// x=x/8(15-yx^2*(10-3*yx^2)
			if (digits == precision + extra) // Reach final iteration step in regards to precision
				{// Reach maximum precision
				r.precision(precision + 1);	// round to final precision
				std::vector<fptype> *mp = r.pointer();
				size_t offset = _float_precision_clz(*mp, 1) + 1 + r.exponent();
				fx = -(double)limit / offset;
				if (offset == 1)// r==c1
					fx = 1.0;
				//cout << "LoopH=" << loopcnt << " Digits=" << formatInt(digits) << ((digits == precision + extra) ? 'T' : 'F') << " Limit=" << formatInt(limit) << " r expo=-" << formatInt(offset) << " ratio=" << fx << endl;  // DEBUG
				if (r == c1 || fx < 2.5)	// break if no improvement
					break;
				if (fx <= 4.0)
					halley = false;		// Close enough to only use 2nd order in the last iteration
				}
			}
		else
			{// Use Newton iterations
			// so we start by assigning it to r, rounding it to the precision of r
			r = y;						// y
			r *= x.square();			// yx^2
			r = c3 - r;					// 3-yx^2
			r.adjustExponent(-1);		// (3-yx^2)/2  	
			x *= r;						// x=x(3-yx^2)/2
			if (digits == precision + extra) // Reach final iteration step in regards to precision
				{// Reach maximum precision
				r.precision(precision + 1);	// round to final precision
				std::vector<fptype> *mp = r.pointer();
				size_t offset = _float_precision_clz(*mp, 1) + 1 + r.exponent();
				fx = -(double)limit / offset;
				if (offset == 1)// r==c1
					fx = 1.0;
				//cout << "LoopN=" << loopcnt << " Digits=" << formatInt(digits) << ((digits == precision + extra) ? 'T' : 'F') << " Limit=" << formatInt(limit) << " r expo=-" << formatInt(offset) << " ratio=" << fx << endl;  // DEBUG
				if (r == c1 || fx < 1.9)	// break if no improvement
					break;
				r.precision(precision + extra);
				}
			}
		}

	x *= y;
	x.adjustExponent(expo_sq);
	// Round to same precision as argument and mrounding mode
	x.mode(a.mode());
	x.precision(precision);
	return x;
	}
/*
float_precision sqrt(const float_precision& a)
	{
	const unsigned int extra = 2;
	const size_t precision = a.precision();
	const eptype expo = a.exponent();
	eptype expo_sq;
	size_t digits, loopcnt;
	double fv;
	float_precision r, x, y(a);
	const float_precision c1(1), c15(15), c3(3), c10(10);

	if (a.iszero() || a == c1)  // Simple square root
		return a;
	if (a.sign() < 0)
		{
		throw float_precision::domain_error();
		}

	if (y.size() == 1 && (expo & 0x1) == 0)
		{// True power of 2 and the expoenent is even
		y.exponent(y.exponent() >> 1); // Half the exponent and return the result.
		return y;
		}
	expo_sq = expo / 2;
	y.exponent(expo - 2 * expo_sq);
	// Do iteration using 2 digits higher precision
	y.precision(precision + extra);
	r.precision(precision + extra);
	x.precision(precision + extra);
	// Get a initial guess using ordinary floating point
	fv = (double)y;				// Convert to double	
								// set the initial guess with at approx 16 correct digits
	fv = 1 / sqrt(fv);
	x = float_precision(fv);

	// Now iterate using Halley
	// Notice y is the original number to squareroot which has the full precision
	for (loopcnt = 1, digits = std::min((size_t)48, precision); ; digits = std::min(precision + extra, digits * 3), ++loopcnt)
		{
		r.precision(digits);
		x.precision(digits);
		r = y * x * x;				// yx^2
		r = (c15 - r*(c10 - c3*r)); // 15-yx^2*(10-3*yx^2)
		r.adjustExponent(-3);		// r=r/8   	
		x *= r;						// x=x/8(15-yx^2*(10-3*yx^2)
		if (digits == precision + extra) // Reach final iteration step in regards to precision
			{
			r.precision(precision + 1);	// round to final precision
			if (r == c1)	// break if no improvement
				break;
			}
		}

	x *= y;
	x.adjustExponent(expo_sq);
	// Round to same precision as argument and mrounding mode
	x.mode(a.mode());
	x.precision(precision);
	return x;
	}
	*/


////////////////////////////////////////////////////////////////////////////////////////
//
// FLOAT PRECISION
//    Universal Constants LN2, LN10, e, _SQRT2, _SQRT3, _INVSQRT2, _INVSQRT3 and PI
//
///////////////////////////////////////////////////////////////////////////////////////


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		17/Jul/2022
//	@brief 		Calculate Binary Splitting of e
//	@return 	Nothing  -	Return p and q as reference arguments
//	@param      "a"	-	   Lower range argument
//	@param      "b"	-	   Upper range argument
//	@param      "p"	-	   Returned argument p
//	@param      "q"	-	   Returned argument q
//
//	@todo  
//
// Description:
//   Use Binary splitting method for e
//   e^1=p(k)/q(k)), where k is the number of Taylor terms needed	
//	Speed up by calculating p & q for range [a..b] <= 4
//
static void binarysplittingE(const uintmax_t a, const uintmax_t b, int_precision& p, int_precision& q)
	{
	int_precision pp, qq;
	uintmax_t mid;

	if (b - a == 1)
		{// No overflow using 64bit arithmetic
		p = int_precision(1);
		q = int_precision(b);
		return;
		}
	if (b - a == 2 /* && b <= 4'294'967'296ull*/)
		{// No overflow using 64bit arithmetic if b<=4'294'967'296
		p = int_precision(b + 1);
		q = int_precision(b*(b - 1));
		return;
		}
	if (b - a == 3 && b <= 2'642'245ull)
		{// No overflow using 64bit arithmetic if b<=2'642'245
		uintmax_t b1 = b*(b - 1);
		p = int_precision(b1 + b + 1);
		q = int_precision(b1*(b - 2));
		return;
		}
	if (b - a == 4 && b <= 65'536ull)
		{// No overflow using 64bit arithmetic if b<=65'536
		uintmax_t b1 = b*(b - 1), b2 = b1*(b - 2);
		p = int_precision(b2 + b1 + b + 1);
		q = int_precision(b2*(b - 3));
		return;
		}
	mid = (a + b) / 2;
	binarysplittingE(a, mid, p, q);  // interval [a..mid]
	binarysplittingE(mid, b, pp, qq);// interval [mid..b]
	// Reconstruct interval [a..b] and return p & q
	p = p*qq + pp;
	q *= qq;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		17/Jul/2022
//	@brief 		Calculate  the number of Taylor terms needed using Stirling approximation 
//	@return 	The number of needed Taylor terms for e
//	@param      "digits"	- The precision in decimal digits
//
//	@todo  
//
// Description:
//   Use Binary splitting method for e
//   Firt calculate the needed number of Taylor terms and then feed that into
//	the binary splitting method for e
// 
static uintmax_t stirling_approx(uintmax_t digits)
	{
	double xnew, xold;
	const double test = (digits + 1) * log(10);
	// Stirling approximation of k!~Sqrt(2*pi*k)(k/e)^k.
	// Taken ln on both side you get: k*(log((k)-1)+0.5*log(2*pi*m);
	// Use Newton method to find in less that 4-5 iteration
	for (xold = 5, xnew = 0; ; xold = xnew)
		{
		double  f = xold*(log(xold) - 1) + 0.5*log(2 * 3.141592653589793 * xold);
		double f1 = 0.5 / xold + log(xold);
		xnew = xold - (f - test) / f1;
		if ((uintmax_t)ceil(xnew) == (uintmax_t)ceil(xold))
			break;
		}
	return (uintmax_t)ceil(xnew);
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		17/Jul/2022
//	@brief 		Calculate  e using the binary splitting method
//	@return 	e			- e for the requested precision
//	@param      "digits"	- The precision in decimal digits
//
//	@todo  
//
// Description:
//   Use Binary splitting method for e
//   First calculate the needed number of Taylor terms and then feed that into
//	the binary splitting method for e
// 
static float_precision computeEdigits(const uintmax_t digits)
	{
	uintmax_t k;
	int_precision p, q;
	float_precision fp, fq;

	fp.precision(digits + 1);
	fq.precision(digits + 1);
	k = stirling_approx(digits);
	if (k < 2)
		k = 2;  // Minimum 2 terms otherwise it cant split
#ifdef HVE_THREAD
	if (digits < 100'000)
		{// Do no threading for digits < 100'000
		binarysplittingE(0, k, p, q);
		}
	else
		{// digits >= 100'000
		// Four threads
		int_precision pp, qq, p3, q3, p4, q4;
		std::thread first([=, &p, &q]()
		{binarysplittingE(0, k/4, p, q);});		// interval [a..k/4]

		std::thread second([=, &pp, &qq]()
		{binarysplittingE(k/4, k/2, pp, qq);});		// interval [k/4..k/2]

		std::thread third([=, &p3, &q3]()
		{binarysplittingE(k/2, 3*k/4, p3, q3);});	// interval [k/2..3k/4]
		
		std::thread fourth([=, &p4, &q4]()
		{binarysplittingE(3*k/4, k, p4, q4);});		// interval [3k/4..k]
		
		first.join();
		second.join();
		// Reconstruct [0..k/2]
		p = p*qq + pp;
		q *= qq;
		third.join();
		fourth.join();
		// Reconstruct [k/2..k]
		p3 = p3*q4 + p4;
		q3 *= q4;
		// Reconstruct [0..k]
		p = p*q3 + p3;
		q *= q3;
		}
#else		
	binarysplittingE(0, k, p, q);
#endif
	
	fp = float_precision(p + q, digits + 1); fq = float_precision(q, digits + 1);
	fp /= fq;
	fp.precision(digits);
	return fp;
	}



// Spigot LN(X/Y) where x and y are integers and x>0 && x>y as conditions

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		28/Jan/2017
//	@brief 		Calculate transcendetal constant of ln(x/y)
//	@return		std::string -	return the constant as a standard std::string
//	@param		"x"	-	The nominator of the number x
// @param		"y"	-	The Denominator of the number x
//	@param		"digits"	-	Number of digits
// @param		"no_dig"	-	Number of digits calculated per loop
//
//	@todo
//
// Description:
//
// 64 bit version of spigot algorithm for LN(x/y) fraction 
// It has automatic 64bit integer overflow detection in which case the result start with the string "Overflow...."
// A Column: x-1,x-1,x-1,...,x-1
// B Column: x,x,x,x,x,...,x
// Initialization values: (x-1)/(x(n+1))...
// The function is declare static since it only serve as a sub function for the function _float_table()
//
static std::string spigot_lnxy_64(const unsigned int x, const unsigned int y, const size_t digits, int no_dig = 1)
	{
	static unsigned long f_table[] = { 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000 };
	bool first_time = true;				// First iteration of the algorithm
	bool overflow_flag = false;			// 64bit integer overflow flag
	//char buffer[32];
	std::string ss;						// The std::string that holds the ln(x)
	size_t dig;
	unsigned int car;
	size_t no_terms;				// No of terms to complete as a function of digits
	unsigned long f;					// New base 1 decimal digits at a time
	unsigned long dig_n;				// dig_n holds the next no_dig digit to add
	unsigned long long carry;
	unsigned long long tmp_n, tmp_dn;
	ss.reserve(digits + 16);
	int factor;

	if (x < y || x < 1)
		{ throw float_precision::domain_error(); }

	if (no_dig > 8) no_dig = 8;			// Ensure no_dig<=8
	if (no_dig < 1) no_dig = 1;			// Ensure no_dig>0
	// Since we do it in trunks of no_dig digits at a time we need to ensure digits is divisble with no_dig.
	dig = (digits / no_dig + (digits%no_dig>0 ? 1 : 0)) * no_dig;
	dig += no_dig;						// Extra guard digits
										// Calculate the number of terms needed
	factor = (int)ceil(10 * log(0.5) / log((double)(x - y) / (double)x));
	no_terms = (unsigned int)(factor * dig / 3 + 3);
	// Allocate the needed accumulators
	unsigned long long *acc_n = new unsigned long long[no_terms + 1];
	unsigned long long *acc_dn = new unsigned long long[no_terms + 1];
	f = f_table[no_dig];				// Load the initial f
	carry = 0;							// Set carry to 0
	//Loop for each no_dig
	for (size_t i = 0; i <= dig && overflow_flag == false; i += first_time == true ? 1 : no_dig, first_time = false)
		{
		// Calculate new number of terms needed
		no_terms = (unsigned int)(factor * (dig-i) / 3 + 3);
		// Loop for each no_terms
		for (size_t j = no_terms; j>0 && overflow_flag == false; --j)
			{
			if (first_time == true)
				{// Calculate the initialize value
				tmp_dn = (j + 1) * x;
				tmp_n = (x - y);
				}
			else
				{
				tmp_n = acc_n[j];
				tmp_dn = acc_dn[j];
				}
			if (tmp_n > (ULLONG_MAX) / f)
				overflow_flag = true;
			tmp_n *= f;		// Scale it
			// Check for 64bit overflow. Not very likely 
			if (carry > 0 && tmp_dn > (ULLONG_MAX - tmp_n) / carry)
				overflow_flag = true;
			tmp_n += carry * tmp_dn;
			carry = (tmp_n / (x * tmp_dn));
			carry *= (x - y);
			acc_n[j] = tmp_n % (tmp_dn * x);
			acc_dn[j] = tmp_dn;
			}

		if (first_time == true)
			{
			tmp_n = (x - y) * f;
			if (carry > 0 && tmp_n > (ULLONG_MAX - carry * x))
				overflow_flag = true;

			acc_n[0] = (tmp_n + carry*x);
			acc_dn[0] = x;
			dig_n = (unsigned)(acc_n[0] / (f*acc_dn[0]));
			}
		else
			{
			if (acc_n[0] > (ULLONG_MAX - carry * acc_dn[0]) / f)
				overflow_flag = true;
			dig_n = (unsigned)((acc_n[0] * f + carry * acc_dn[0]) / (f*acc_dn[0]));
			}
		car = (unsigned)(dig_n / f);
		dig_n %= f;
		// Add the carry to the existing number for digits calculate so far.
		if (car > 0)
			{
			for (size_t j = ss.length(); car > 0 && j > 0; --j)
				{
				unsigned int dd;
				dd = (ss[j - 1] - '0') + car;
				car = dd / 10;
				ss[j - 1] = dd % 10 + '0';
				}
			}
		ss += uitostring10((uintmax_t)dig_n, first_time==true? 1:no_dig );  //need to use uitostring10 instead of snprintf for performance
		//(void)snprintf(buffer, sizeof(buffer), "%0*lu", first_time == true ? 1 : no_dig, dig_n);
		//ss += std::string(buffer);
		if (first_time == true)
			acc_n[0] %= f*acc_dn[0];
		else
			{
			acc_n[0] = acc_n[0] * f + carry *acc_dn[0];
			acc_n[0] %= f  * acc_dn[0];
			}
		carry = 0;
		}

	ss.insert(1, ".");// add a . come after the first digit to create 2.30...
	if (overflow_flag == false)
		ss.erase(digits + 1); // Remove the extra digits that we didnt requested.
	else
		ss = std::string("Overflow:") + ss;

	delete[] acc_n;
	delete[] acc_dn;
	return ss;
	}

// End Spigot LN(X/Y)

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		19/Aug/2022
//	@brief 		Calculate he Chudnovsky PI using the Binary Splittting method
//	@return		float_precision -	return the constant pi
//	@param		"a"	-	The lower index a
//	@param		"b"	-	The upper index b
//	@param		"p"	-	The P(k) parameter in the binary splitting method
//	@param		"q"	-	The Q(k) parameter in the binary splitting method
//	@param		"r" -	The R(k) parameter in the binary splitting method
//
//	@todo
//
// Description:
//
// Implement the binary splitting method for Chudnovsky PI method
// It call the function recursively until the index (a+1==b) 
// The function is declare static since it only serve as a sub function for the function in _float_table()
//
static void binarysplittingChudnovskiPI(const uintmax_t a, const uintmax_t b, int_precision& p, int_precision& q, int_precision& r)
	{
	int_precision pp, qq, rr;
	uintmax_t mid;

	if (b - a == 1)
		{// No overflow using 64bit arithmetic
		if(b<635'130)
			r = int_precision( (2 * b - 1)*(6 * b - 5)*(6 * b - 1));
		else
			{// No overflow if b<=715'827'883 ~ 10B digits
			r = int_precision((6 * b - 5)*(6 * b - 1)); r *= int_precision(2*b - 1);
			}
		p = int_precision(13'591'409ull + 545'140'134ull * b); p *= r;
		if (b & 0x1)
			p.change_sign();
		if (b <= 2'642'245)
			q = int_precision(b*b*b);  // No overflow in 64bit environment
		else
			{
			q = int_precision(b*b); q *= int_precision(b);
			}
		q *= int_precision(10'939'058'860'032'000ull);
		return;
		}
	mid = (a + b) / 2;
	binarysplittingChudnovskiPI(a, mid, p, q, r);	// interval [a..mid]
	binarysplittingChudnovskiPI(mid, b, pp, qq, rr);// interval [mid..b]
	// Reconstruct interval [a..b] and return p, q & r
	p = p*qq + pp*r;
	q *= qq;
	r *= rr;
	}

//End Chudnovsky Binary Splitting 

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/26/2021
//	@brief 		Lookup or generate "fixed" constant ln2, PI log10 etc
//	@return 	float_precision	-	return the new table lookup value
//	@param		"tt"	-	Which table type lookup is needed
//	@param		"precision"	-	Number of significant digits
//
//	@todo
//
// Description:
//   Dynamic tables for "fixed" constant like ln(2), ln(10), e, PI and 1/Sqrt(2), sqrt(2)
//   If a higher precision is requested we create it and return otherwise 
//   we just the "constant" at a higher precision which eventually will be
//   rounded to the destination variables precision 
//	Added 1/sqrt(2) and sqrt(2) as constants
//
float_precision _float_table( enum table_type tt, size_t precision )
   {
   static float_precision ln2( 0, 0, ROUND_NEAR );
   static float_precision ln10( 0, 0, ROUND_NEAR );
   static float_precision pi( 0, 0, ROUND_NEAR );
   static float_precision e( 0, 0, ROUND_NEAR);
   static float_precision invsqrt2(0, 0, ROUND_NEAR);
   static float_precision invsqrt3(0, 0, ROUND_NEAR);
   float_precision res(0, precision, ROUND_NEAR);

	switch( tt )
		{
		case _EXP1:
			if (e.precision() >= precision)
				res = e;
			else
				{// Using Binary Splitting Algorithm for exp(1) Calculation with threading
				std::string ss;
				size_t prec = std::max((size_t)20U, precision + 2);
				e.precision(prec);
				e = computeEdigits( prec );		// Calculate exp(1)
				res = e;
				}
			break;
		case _LN2:
			if (ln2.precision() >= precision)
				res = ln2;
			else
				{ // Using Spigot algorithm for LN2 Calculation
				std::string ss;
				size_t prec = std::max((size_t)20U, precision+2);
				ln2.precision( prec );
				ss = spigot_lnxy_64(2, 1, prec, 4);	// The result as a string in BASE_10
				ln2 = float_precision(ss, prec);		// Convert to float_precision
				ln2.precision(std::max((size_t)20U, precision));
				res = ln2;								// Save the result
				}
			break;
		case _LN10:
			if( ln10.precision() >= precision )
				res = ln10;
			else
				{ // Using Spigot Algorithm for LN10. LN(10)=3*ln(2)+ln(10/8)
				std::string ss;
				size_t prec=std::max((size_t)20U, precision + 2);
				ln10.precision(prec);
				ss = spigot_lnxy_64(2, 1, prec, 4);		// The result as a string in BASE_10
				ln10 = float_precision(ss, prec);		// Convert to float_precision
				ln10 *= float_precision(3);
				ss = spigot_lnxy_64(10, 8, prec, 4);	// The result as a string in BASE_10
				ln10 += float_precision(ss, prec);		// Convert and add to float_precision ln10
				ln10.precision(std::max((size_t)20U, precision));
				res = ln10;								// Save the result
				}
			break;
		case _PI: 
			if( pi.precision() > precision )
			res = pi;
			else
				{// Binary splitting with recursion
				size_t prec = std::max((size_t)20U, precision + 2);
				const uintmax_t k = (uintmax_t)ceil(prec*log(10) / log(151931373056000ull)); //(uintmax_t)ceil(precision / 14.18); // 
				int_precision p, q, r;
#ifdef HVE_THREAD
				if (k <= 700)
					{// No threads when k<=700
					binarysplittingChudnovskiPI(0, k, p, q, r);
					}
				else
					{
					int_precision pp, qq, rr;
					int_precision p3, q3, r3;
					int_precision p4, q4, r4;

					std::thread first([=, &p, &q, &r]()
					{binarysplittingChudnovskiPI(0, k / 4, p, q, r); });
					std::thread second([=, &pp, &qq, &rr]()
					{binarysplittingChudnovskiPI(k / 4, k / 2, pp, qq, rr); });
					std::thread third([=, &p3, &q3, &r3]()
					{binarysplittingChudnovskiPI(k / 2, 3 * k / 4, p3, q3, r3); });
					std::thread fourth([=, &p4, &q4, &r4]()
					{binarysplittingChudnovskiPI(3 * k / 4, k, p4, q4, r4); });
					// Wait for the thread to finish
					first.join();
					second.join();
					// Reconstruct interval [a..(a+b)/2] and update p, q & r
					p = p*qq + pp*r;
					q *= qq;
					r *= rr;

					// Wait for the thread to finish
					third.join();
					fourth.join();
					// Reconstruct interval [(a+b)/2..b] and update p3, q3 & r3
					p3 = p3*q4 + p4*r3;
					q3 *= q4;
					r3 *= r4;
					// Reconstruct interval [a..b] and update p, q & r
					p = p*q3 + p3*r;
					q *= q3;
					r *= r3;
					}

#else
				pi.precision(precision+1);
				binarysplittingChudnovskiPI(0, k, p, q, r);
#endif
				pi.precision(prec);
				pi = (float_precision(4270934400ull) * float_precision(q, prec)) / (float_precision(p, prec) + float_precision(13591409ull)*float_precision(q, prec));
				pi *= 1 / sqrt(float_precision(10005, prec));
				pi.precision(precision);
				res = pi;
				}
			break;
		case _INVSQRT2:
			if (invsqrt2.precision() > precision)
			res = invsqrt2;
			else
				{
				unsigned int loopcnt = 0;
				const unsigned int extra = 2+5;
				const float_precision c1(1), c3(3);
				const intmax_t limit = -(intmax_t)((precision + 1)*log2(10)) - 1;
				size_t digits;
				float_precision r;
	 
				digits = invsqrt2.precision();  // Get current precision
				if (invsqrt2.iszero() == true)  // First time, do initialization 
					{
					digits = 16;
					invsqrt2.precision(std::max(precision, digits));  // Ensure minimum as 16 decimal digits
					// New. Get a initial guess using ordinary floating point
					// set the initial guess with at approx 16 correct digits
					invsqrt2 = float_precision(1.0 / sqrt(2.0),digits); // Ensure same precision as standard IEEE754
					}
				precision = std::max(precision, digits); // Keep maxium precision already calculated
				invsqrt2.precision(precision);	
#ifdef HVE_DEBUG
				std::cout << "INVSQRT2 Max precision=" << precision+extra << " start Precision=" << digits*2 << std::endl;  // Debug
				int tadd = clock();
#endif
				// Now iterate using Netwon x=0.5x(3-2x^2), where x=invsqrt2
				for (digits *= 2; ; digits = std::min(precision + extra, digits * 2), ++loopcnt)
					{
					// Increase precision by a factor of two for the working variable s r & u. 
					r.precision(digits);
					invsqrt2.precision(digits);
					// Notice 2 is the original number to squareroot which has the full precision 
					// so we start by assigning it to r, rounding it to the precision of r
					r = invsqrt2.square();		// x^2
					r.adjustExponent(+1);		// 2x^2
					r = c3 - r;					// 3-2x^2
 					r.adjustExponent(-1);		// r *= c05; 
					invsqrt2 *= r;				// (3-2x^2)/2
#ifdef HVE_DEBUG
					tadd = clock() - tadd;
					std::cout << "\tINVSQRT2 Iteration=" << loopcnt << " Precision="<< digits << " Error exponent=" << ((r - c1).exponent()) / log2(10) << " Time=" << tadd / CLOCKS_PER_SEC << std::endl;  // Debug
					tadd = clock();
#endif
					if (digits == precision + extra) // Reach final iteration step in regards to precision
						{/*
  						r.precision(precision + 1);	// round to final precision
						if (r == c1)	// break if no improvement
							break;
						r.precision(precision + extra);
						*/
						r.precision(precision + 1);	// round to final precision
						std::vector<fptype> *mp = r.pointer();
						size_t offset = _float_precision_clz(*mp, 1) + 1 + r.exponent();
						double ratio = -(double)limit / offset;

						if (r == c1 || ratio < 1.9)			// break if no improvement
							break;
						r.precision(precision + extra);
						}
					}

				// Round to same precision as argument
				invsqrt2.precision(precision);
				res = invsqrt2;
#ifdef HVE_DEBUG
				std::cout << "INVSQRT2 finish with precision=" << precision << std::endl;  // Debug
#endif
				}
			break;
		case _SQRT2:  // use 2/sqrt(2)=sqrt(2)
			res = _float_table(_INVSQRT2, precision); //float_precision(2) * _float_table(_INVSQRT2, precision);
			res.adjustExponent( +1);	// Faster way to multiply by 2
			break;
		case _INVSQRT3:
			if (invsqrt3.precision() > precision)
				res = invsqrt3;
			else
				{
				unsigned int loopcnt = 0;
				const unsigned int extra = 2 + 5;
				const intmax_t limit = -(intmax_t)((precision + 1)*log2(10)) - 1;
				const float_precision c1(1), c3(3);
				size_t digits;
				float_precision r;

				digits = invsqrt3.precision();  // Get current precision
				if (invsqrt3.iszero() == true)  // First time, do initialization 
					{
					digits = 16;
					invsqrt3.precision(std::max(precision, digits));  // Ensure minimum as 16 decimal digits
					// New. Get a initial guess using ordinary floating point
					// set the initial guess with at approx 16 correct digits
					invsqrt3 = float_precision(1.0 / sqrt(3.0), digits); // Ensure same precision as standard IEEE754
					}
				precision = std::max(precision, digits); // Keep maxium precision already calculated
				invsqrt3.precision(precision);
				// Now iterate using Netwon Un=0.5U(3-2U^2), where U=invsqrt2
				for (digits *= 2; ; digits = std::min(precision + extra, digits * 2), ++loopcnt)
					{
					// Increase precision by a factor of two for the working variable s r & u. 
					r.precision(digits);
					invsqrt3.precision(digits);
					// Notice V is the original number to squareroot which has the full precision 
					// so we start by assigning it to r, rounding it to the precision of r
					r = c3;						// 3
					r *= invsqrt3.square();		// 3x^2
					r = c3 - r;					// 3-3x^2
 					r.adjustExponent(-1);		// r *= c05;	// (3-3x^2)/2
					invsqrt3 *= r;				// x=x(3-3x^2)/2
					if (digits == precision + extra) // Reach final iteration step in regards to precision
						{
						/*
						r.precision(precision + 1);	// round to final precision
						if (r == c1)	// break if no improvement
							break;
						r.precision(precision + extra);
						*/
						r.precision(precision + 1);	// round to final precision
						std::vector<fptype> *mp = r.pointer();
						size_t offset = _float_precision_clz(*mp, 1) + 1 + r.exponent();
						double ratio = -(double)limit / offset;

						if (r == c1 || ratio < 1.9)			// break if no improvement
							break;
						r.precision(precision + extra);
						}
					}

				// Round to same precision as argument
				invsqrt3.precision(precision);
				res = invsqrt3;
				}
			break;
		case _SQRT3:  // use 3/sqrt(3)=sqrt(3)
			res = float_precision(3) * _float_table(_INVSQRT3, precision);
			break;
		case _ONETENTH:
			{// Generate the constant 0.1 with the requested percision
			size_t n = (size_t)ceil(precision * log2(10) / Bitsfptype);
			res.precision(precision + 1);
			std::vector<fptype> v(0);
			for (v.push_back(1); n > 0; --n)
				v.push_back(0x9999999999999999ull);
			res.number(v);
			res.precision(precision); // Force normalize and rounding
			res.exponent(-4);
			break;
			}
		default:
			ln2=pi= e = float_precision(0);
			e.precision(0);
			pi.toPrecision(0);
			ln2.precision(0);
			break;
		}

   return res;
   }

///////////////////////////////////////
//
// FLOAT PRECISION FUNCTIONS
//    Exp(), Log(), Log10(), Log2()
//
///////////////////////////////////////


//	Experimental new exp() using sinh() and sqrt()
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/Aug/2013, 10/May/2022
//	@brief 		Calculate exp(x)
//	@return 	float_precision -	Return exp(x)
//	@param      "x"	-	   The argument
//
//	@todo  
//
// Description:
//   Use a the identity that exp(x)=sinh(x)+sqrt(1+sinh(x)^2)
//	  This has proven to be faster than the standard taylor series for exp()
//   exp(x) == 1 + x + x^2/2!+x^3/3!+....
//	  A test is made for x is an integer in which case we do pow(e,x) which is more than 400 times faster
//    sine e is calculated using the spigot algorithm using pure 64bit integer arithmetic
//
float_precision exp(const float_precision& x)
	{
	size_t precision = x.precision() + 2 + (size_t)ceil(log10(x.precision()));
	float_precision v(x);
	const float_precision c1(1);

	v.precision(precision);
	if (v.sign() < 0)
		v.change_sign();

	if (floor(v) == v)  // v is an Integer
		{// use the 100 times faster shortcut exp(v)=exp(1)^v
		v = _float_table(_EXP1, precision);
		v = pow(v, abs(x));
		}
	else
		{
		v = sinh(v);
		//v.precision(precision+10);  // Double the precision to avoid loss of significant when performaing 1+v*v
		v += sqrt(c1 + v.square());
		v.precision(precision);
		}

	if (x.sign() < 0)
		v = _float_precision_inverse(v);
	// Round to same precision as argument and rounding mode
	v.mode(x.mode()); v.mode(x.mode());
	v.precision(x.precision());
	return v;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		20/May/2022
//	@brief 		Calculate log(x) using the AGM method
//	@return 	float_precision -	Return log(x)
//	@param      "x"	-	   The argument
//
//	@todo  
//
// Description:
//	Use the AGM methof to calculate log(x) for digits > 4,0000 this function exceed the performance
//	of the standard Taylor serie version even when argument reduction and coefficient scalling is added to later method.
//
static float_precision logAGM(const float_precision& x)
	{
	const size_t guard = 5;
	const size_t precision = x.precision() + (size_t)ceil(log10(x.precision())) + guard;
	const uintmax_t s = (uintmax_t)ceil(precision*log(10) / (2 * log(2)) + 1 - log((double)x) / log(2));
	const uintmax_t slost = (uintmax_t)ceil(log((double)s) / log(10.0));  // Loss of precision
	const float_precision c1(1);
	float_precision logx, z(x), agm, ln2;

	if (x <= float_precision(0))
		{
		throw float_precision::domain_error();
		}

	// Adjust to working precision
	agm.precision(precision);
	logx.precision(precision);
	z.precision(precision);
	ln2.precision(precision);

	if (z < c1)
		z = 1 / z; // Now z >= 1
	z.adjustExponent(s - 2);		// z/=4

#ifdef HVE_THREAD
	// Threaded version
	// First st thread calcuate PI
	std::thread first([=, &logx]()
	{ logx = _float_table(_PI, precision); });
	// Second thread calculate ln(2)
	std::thread second([=, &ln2]()
	{ ln2 = _float_table(_LN2, precision + slost); });
	// Third thread calculate AGM(z,1)
	std::thread third([=, &agm, &z]()
	{ agm = AGM(z, 1);
	agm.adjustExponent(+1);  //2*AGM
	});
	// Wait for th thread 1 & 3 finish
	first.join();
	third.join();

	logx *= z;		//PI*x^(s-2)
	logx /= agm;  // PI*x^(s-2)/AGM

	// Wait for ln(2) thread to finish
	second.join();
	// Increase precision to avoid loos of significants when subtrating two large numbers
	logx.precision(precision + slost);
	logx -= float_precision(s, precision + slost) * ln2;
#else
	// Non threaded version
	logx = _float_table(_PI, precision);
	logx *= z;			//PI*x^(s-2)
	agm = AGM(z, 1);	// PI*x^(s-2)/AGM
	agm.adjustExponent(+1);  //2*AGM
	logx /= agm;
	// Increase precision to avoid loos of significants when subtrating two large numbers
	logx.precision(precision + slost);
	logx -= float_precision(s, precision + slost) *_float_table(_LN2, precision + slost);
#endif

	// Round to same precision as argument and rounding mode
	if (x < c1)
		logx.change_sign();
	logx.mode(x.mode());
	logx.precision(x.precision());
	return logx;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/Oct/2021, 20/May/2022
//	@brief 		Calculate log(x) using Taylor series
//	@return 	float_precision -	Return log(x)
//	@param      "x"	-	   The argument
//
//	@todo  
//
// Description:
//	 Simplified and improved with the new binary storage
//   Use a taylor series until their is no more change in the result
//   Equivalent with the same standard C function call
//   ln(x) == 2( z + z^3/3 + z^5/5 ...
//   z = (x-1)/(x+1)
//	Argument reduction is done to 1.001 level
//	Coefficients rescaling is done for group=5
//
static float_precision logTaylor(const float_precision& x)
	{
	const int group = 5;
	size_t precision = x.precision() + 2 + (size_t)ceil(log10(x.precision()));
	eptype expo;
	int i, k, no_reduction;
	size_t loopcnt = 1;
	double zd;
	float_precision logx, z(x), zsq, terms;
	const float_precision c1(1);

	if (x <= float_precision(0))
		{
		throw float_precision::domain_error();
		}

	expo = z.exponent();	// Get original exponent
	z.exponent(0);			// Set exponent to zero getting a z number beween [1..2)

	// Check for augument reduction and increase precision if necessary
	zd = (double)z;
	no_reduction = (int)ceil(log(log(zd) / log(1.001)) / log(2));
	no_reduction = std::max(no_reduction, 0);
	precision += no_reduction;
	precision = std::max((size_t)20, precision);

	// adjust precision to allow correct rounding of result
	z.precision(precision);
	zsq.precision(precision);
	terms.precision(precision);
	logx.precision(precision);

	// In order to get a fast Taylor series result we need to get the fraction closer to one
	// The fraction part is [1...1.1) (base 10) at this point
	// Repeat a series of no_reduction square root 
#define HVE_USE_NROOT
#ifdef HVE_USE_NROOT
	if (no_reduction > 0) { z = nroot(z, 2 << (no_reduction - 1)); }
	k = no_reduction;
#else
	for (k = 0; k < no_reduction; ++k)
		z = sqrt(z);
#endif

	// number now at [1...1.1). Setup the iteration
	z = (z - c1) / (z + c1);
	zsq = z.square();
	logx = z;
	if (group == 1)
		{
		// Iterate using taylor series ln(x) == 2(z + z^3/3 + z^5/5 ... )
		for (i = 3;; i += 2, ++loopcnt)
			{
			z *= zsq;
			terms = z / float_precision(i);
			if (logx + terms == logx)
				break;
			logx += terms;
			}
		}
	else
		{
		std::vector<float_precision> vn(group);
		std::vector<float_precision> cn(group);
		int j, l;

		for (i = 0; i < group; ++i)
			{
			vn[i].precision(precision);
			cn[i].precision(precision);
			if (i == 0) vn[i] = zsq;
			if (i > 0) vn[i] = vn[i - 1] * zsq;
			}
		// Now iterate 
		for (i = 3; ; )
			{
			for (j = 0; j < group; ++j)
				{
				cn[j] = c1;
				for (l = 0; l < group; ++l)
					if (j != l)
						cn[j] *= i + 2 * l;
				}
			for (j = 0, terms = 0; j < group; ++j)
				terms += cn[j] * vn[j];
			terms *= z / (cn[0] * float_precision(i, precision));
			i += 2 * group;				// Update term count
			loopcnt += group;
			if (logx + terms == logx)	// Reach precision
				break;					// yes terminate loop
			logx += terms;				// Add taylor terms to result
			if (group > 1)
				z *= vn[group - 1];		// ajust z to last Taylor term
			}
		}


	// Adjust result from the reduction by multiply it with 2^(k+1)
	logx *= float_precision(pow(2.0, (double)(k + 1)));
	if (expo != 0)  // Adjust for original exponent
		{// Ln(x^y) = Ln(x) + Ln(2^y) = Ln(x) + y * ln(2) 
		logx += float_precision(expo) * _float_table(_LN2, precision + 1);
		}

	// Round to same precision as argument and rounding mode
	logx.mode(x.mode());
	logx.precision(x.precision());
	return logx;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		20/May/2022
//	@brief 		Calculate log(x)
//	@return 	float_precision -	Return log(x)
//	@param      "x"	-	   The argument
//
//	@todo  
//
// Description:
//	 Call either LogTaylor() or logAGM() depending on the precision of x
//
float_precision log(const float_precision& x)
	{
	if (x.precision() < 4'000)
		return logTaylor(x);
	return logAGM(x);
	}



//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		Calculate log10(x)
//	@return 	float_precision -	Return log10(x)
//	@param      "x"	-	   The argument
//
//	@todo  
//
// Description:
//   Log10. Use the equation log10(x)=log(x)/log(10)
//   Equivalent with the same standard C function call
//
float_precision log10( const float_precision& x )
	{
	size_t precision = x.precision();  
	float_precision res( 0, precision + 1 );

	if( x <= float_precision(0) ) 
		{ throw float_precision::domain_error(); }

	res = x;
	res = log( res ) / _float_table( _LN10, precision + 1 );
   
	// Round to same precision as argument and rounding mode
	res.mode( x.mode() );
	res.precision( x.precision() );  
	return res;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		20-May-2022
//	@brief 		Calculate log2(x)
//	@return 	float_precision -	Return log2(x)
//	@param      "x"	-	   The argument
//
//	@todo  
//
// Description:
//   Log10. Use the equation log2(x)=log(x)/log(2)
//   Equivalent with the same standard C function call
//
float_precision log2(const float_precision& x)
	{
	size_t precision = x.precision();
	float_precision res(0, precision + 1);

	if (x <= float_precision(0))
		{
		throw float_precision::domain_error();
		}

	res = x;
	res = log(res) / _float_table(_LN2, precision + 1);

	// Round to same precision as argument and rounding mode
	res.mode(x.mode());
	res.precision(x.precision());
	return res;
	}


///////////////////////////////////////
//
// FLOAT PRECISION FUNCTIONS
//    Special functions: 
//		pow()
//		fmod()
//		floor()
//		ceil()
//		modf()
//		fabs()
//		ldexp()
//		frexp()
//
///////////////////////////////////////


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		4/Oct/2021
//	@brief 		Calculate pow(x,y)
//	@return 	float_precision -	Return pow(x,y)
//	@param      "x"	- The argument
// @param      "y" - The power argument
//
//	@todo  
//
// Description:
//   x^y == exp( y * ln( x ) ) ); in general, however if y is an integer then we use the ipow() algorithm instead.
//   Update to use the new method .toInteger()
// 
float_precision pow( const float_precision& x, const float_precision& y )
   {
   float_precision i, res(1);
   eptype expo;
   bool yinteger=false;

   // add two extra guard digits to avoid loss of precision when performing  exp( y * ln(x) ) )
   res.precision( x.precision()+2 );  

   expo = y.exponent();
   if( expo >= 0 )
      {
      i.precision( y.precision() );
      i = y;
	  i.toInteger(); // now i is the integer part of y. 
	  // Check that y is a true integer, with a max range of a 32 bit integer
	  if( y == i && abs(i) <= float_precision( LLONG_MAX ) )
		  yinteger = true;
      }
   
   if( yinteger == false ) // y is not an integer so do x^y= exp^(y*log(x)) the regular way
      {
	  res = x;
      res = log( res ) * y;
      res= exp( res );
      }
   else
		{ // raise to the power of y when y is an integer. Use optimized method.
		int sign = i.sign();
		if( sign < 0 )
			i.change_sign();
		float_precision p(x);

		for( uintmax_t n = (uintmax_t)i; n > 0; n >>= 1 ) 
			{
			if( ( n & 0x1 ) != 0 ) 
				res *= p;		// Odd
			if(n>1)	p *= p;		// Square it						 
			}
		if( sign < 0 )
			res = _float_precision_inverse( res );
		 }
   
   res.precision(x.precision());
   return res;
   }
 
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		27/Sep/2021
//	@brief 		Calculate fmod(x,y)
//	@return 	float_precision -	Return fmod(x,y)
//	@param      "x"	- The argument
// @param      "y"   - The argument
//
//	@todo  
//
// Description:
//   float precision. fmod remainder of x/y
//   Equivalent with the standard C function fmod
//   x = i * y + f or f = x - i * y; and i = integer(x/y)
//	  Revised to use the method.toInteger() for faster calculation
//
float_precision fmod( const float_precision& x, const float_precision& y )
   {
   float_precision i, f;
   
   f.precision( x.precision() );
   i.precision( x.precision() );
   i = x / y;
   if( i.exponent() < 0 )
		f = x;
   else
		{
		i.toInteger();
		f = x - i * y;
		}	

   return f;
   }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		5/Oct/2021
//	@brief 		Calculate floor(x)
//	@return 	float_precision -	Return floor(x)
//	@param      "x"	- The argument
//
//	@todo  
//
// Description:
//   Float Precision floor
//   Equivalent with the same standard C function floor()
//	  Rounds x downward, returning the largest integral value that is not greater than x.
//
float_precision floor( const float_precision& x )
	{
	float_precision f(0, x.precision() );
	const float_precision c1(1);

	if( x.exponent() < 0 ) // is number less than |1|
		{
		if( x.sign() < 0 )
			f = -c1;
		}
	else
		{
		f = x;
		f.toInteger();
		if (f.sign() < 0 && (x - f).iszero() == false )
			f -= c1;
		}
 
	return f;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/Oct/2021
//	@brief 		Calculate ceil(x)
//	@return 	float_precision -	Return ceil(x)
//	@param      "x"	- The argument
//
//	@todo  
//
// Description:
//   Float Precision ceil
//   Equivalent with the same standard C function ceil()
//   Rounds x upward, returning the smallest integral value that is not less than x.
//
float_precision ceil( const float_precision& x )
	{
	float_precision f(0, x.precision() );
	const float_precision c1(1);

	if( x.exponent() < 0 ) // is number less than |1|
		{
		if( x.sign() > 0 )
			f = c1;
		}
	else
		{
		f = x;
		f.toInteger();
		if (f.sign() > 0 && (x - f).iszero() == false)
			f += c1;
		}

	return f;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/Oct/2021
//	@brief 		Calculate trunc(x)
//	@return 	float_precision -	Return trunc(x)
//	@param      "x"	- The argument
//
//	@todo  
//
// Description:
//   Float Precision trunc
//   Equivalent with the same standard C function trunc()
//   Rounds x towards zero.
//
float_precision trunc(const float_precision& x)
	{
	float_precision f(0, x.precision());

	if (x.exponent() < 0) // is number less than |1|
		{
		if (x.sign() > 0)
			f = float_precision(0);
		}
	else
		{
		f = x;
		f.toInteger();
		}

	return f;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/Oct/2021
//	@brief 		Calculate round(x)
//	@return 	float_precision -	Return round(x)
//	@param      "x"	- The argument
//
//	@todo  
//
// Description:
//   Float Precision round
//   Equivalent with the same standard C function round()
//
float_precision round(const float_precision& x)
	{
	float_precision f(x);
	const float_precision c05(0.5);

	f += (f.sign() < 0 ? -c05 : c05 );
	f = trunc(f);

	return f;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		30/Sep/2021
//	@brief 		Split a float number into integer part and fraction
//	@return 	float_precision -	Fraction part of x
//	@param      "x"	- The argument
// @param      "intptr" - Float_precision pointer to integer part of x
//
//	@todo  
//
// Description:
//   Float Precision fmod
//   Split a Floating point number into integer part and fraction
//   Equivalent with the same standard C function call
//	  Use a modified version of the function fmod(x,c1)
//	  Converted to binary version
//   Notice that if x < 0 then BOTH the return value and intptr becomes negative
//	  even when inptr==0
//
float_precision modf( const float_precision& x, float_precision *intptr )
	{
	float_precision i, f;
	eptype expo;

	f.precision(x.precision());
	i.precision(x.precision());
	intptr->precision(x.precision());
	f = x;
	expo = f.exponent();
	if (expo < 0)
		{
		i = float_precision(0);
		i.sign(x.sign());
		}
	else
		{
		i=f.toFraction();
		//i.toInteger();
		//f = x - i;
		}
	*intptr = i;
	return f;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/25/2012
//	@brief 		Calculate abs(x)
//	@return 	float_precision -	Return absolute value of x
//	@param      "x"	- The argument
//
//	@todo  
//
// Description:
//   Float Precision abs()
//   Equivalent with the same standard C function call fabs()
//
float_precision abs( const float_precision& x )
   {
   float_precision f(0, x.precision() );

   f = x;
   if( f.sign() < 0 )
      f.change_sign();

   return f;
   }



//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		Calculate ldexp((x)
//	@return 	float_precision -	Return ldexp(x)
//	@param      "x"	- The argument
// @param      "exp" - exponent argument
//
//	@todo  
//
// Description:
//   The ldexp function returns the value of x * 2^exp
//
float_precision ldexp( const float_precision& x, eptype exp )
   {
   if( exp == 0 )
      return x;
   if( exp > 0 && exp <= 63 )
      return x * float_precision( 1U << exp );

   return x * pow( float_precision( 2 ), float_precision( exp ) );
   }



//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		Calculate frexp(x,expptr)
//	@return 	float_precision -	Return mantissa part of number x
//	@param      "x"	- The argument
// @param      "expptr" - Pointer to the exponent part of number
//
//	@todo  
//
// Description:
//   The frexp()
//   The frexp function breaks down the floating-point value (x) into a mantissa (m) and an exponent (n), 
//   such that the absolute value of m is greater than or equal to 1/RADIX and less than RADIX, and x = m*Radix^n. 
//   The integer exponent n is stored at the location pointed to by expptr. 
//
float_precision frexp( float_precision& x, eptype *expptr )
   {
   if( x.iszero() )
      *expptr = 0;

   *expptr = x.exponent();
   x.exponent( 0 );

   return x;
   }



///////////////////////////////////////
//
// TRIGONOMETRIC FUNCTIONS
//
//   asin()
//	 acos()
//	 atan()
//   atan2()
//   sin()
//   cos()
//   tan()
//
///////////////////////////////////////

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/26/2013, 24/Jun/2022
//	@brief 		Calculate asin(x)
//	@return 	float_precision -	Return asin(x)
//	@param      "x"	-	   The argument
//
//	@todo  
//
// Description:
//   Use a Taylor series until their is no more change in the result
//   asin(x) == x + x^3/(2*3)+(1*3)x^5/(2*4*5)+(1*3*5)x^7/(2*4*6*7)....
//   Use argument reduction via the identity arcsin(x)=2arcsin(x)/(sqrt(2)+sqrt(1-x*x))
//	and use Taylor terms grouping of 5 Taylor terms at a time.
//	asin(0)==0
//
float_precision asin(const float_precision& x)
	{
	const int group = 5;
	size_t precision = x.precision() + 2 + (size_t)ceil(log10(x.precision()));
	intmax_t k, i;
	size_t loopcnt = 1;
	int sign;
	float_precision r, asinx, v(x), vsq, lc, uc, terms;
	const float_precision c1(1), c2(2);

	if (x > c1 || x < -c1)
		{
		throw float_precision::domain_error();
		}

	if (v.iszero())
		return v;
	sign = v.sign();
	if (sign < 0)
		v.change_sign();

	// Automatically calculate optimal reduction factor as a power of two
	k = 2 * (intmax_t)ceil(log(2)*log(precision));
	k = std::min((intmax_t)30, k);  // Top of a maximum 30 reduction. only relevant if the multipier for k>2
									// Adjust k for final value of v when v is small (less than 1). we know it is in the interval between [0..1]
									// This indicate that the exponent is in the range [-inf..0] 
	k += v.exponent();				// Avoid uncessary argument reduction if v is small
	k = std::max((intmax_t)0, k);
	
	// Adjust the precision
	precision += k / 4;
	r.precision(precision);
	asinx.precision(precision);
	v.precision(precision);
	vsq.precision(precision);
	lc.precision(precision);
	uc.precision(precision);
	terms.precision(precision);

	// Now use the identity arcsin(x)=2arcsin(x/(sqrt(2)*sqrt(1+sqrt(1-x*x)))
	// k number of times
	r = _float_table(_SQRT2, precision);
	for (i = 0; i < k; ++i)
		v /= r * sqrt(c1 + sqrt(c1 - v.square()));

	vsq = v.square();
	r = v;
	asinx = v;
	if (group == 1)
		{
		// Now iterate using taylor expansion
		for (i = 3;; i += 2, ++loopcnt)
			{
			if (i < (uintmax_t)1 << 32)
				{// Multiplication fit into 64bit
				uc = float_precision((i - 2) * (i - 2));
				lc = float_precision(i * i - i);
				}
			else
				{
				uc = float_precision(i - 2); uc = uc.square();
				lc = float_precision(i - 1) * float_precision(i);
				}
			r *= uc * vsq / lc;
			if (asinx + r == asinx)
				break;
			asinx += r;
			}
		}
	else
		{
		std::vector<float_precision> vn(group);  // vn[0] is not used
		std::vector<float_precision> un(group);  //
		std::vector<float_precision> ln(group);

		for (i = 0; i < group; ++i)
			{
			un[i].precision(precision); ln[i].precision(precision);  vn[i].precision(precision);
			if (i == 1) vn[1] = vsq;
			if (i > 1) vn[i] = vn[i - 1] * vsq;
			}
		// Now iterate 
		for (i = 3;; )
			{
			// Recalulate the coefficients
			intmax_t j; uintmax_t tmp;
			for (j = group - 1; j >= 0; --j)
				{
				if (j == group - 1)
					{
					tmp = (i - 2); tmp *= tmp;
					uc = float_precision(tmp); un[j] = uc;
					tmp = (i - 1 + j * 2); tmp *= tmp + 1;
					lc = float_precision(tmp); ln[j] = lc;
					}
				else
					{
					tmp = i - 4 + (group - j) * 2; tmp *= tmp;
					uc = float_precision(tmp);
					un[j] = uc;
					tmp = i - 1 + j * 2; tmp *= tmp + 1;
					lc = float_precision(tmp);
					ln[j] = lc;
					un[j] *= un[j + 1];	ln[j] *= ln[j + 1];
					}
				}

			ln[0] = ln[0].inverse();
			// Adding from smalles to largest number
			uc = terms = vn[group - 1] * un[0];
			for (j = group - 1; j >= 2; --j)
				terms += (un[group - j] * ln[j]) * vn[j - 1];
			terms += un[group - 1] * ln[1];
			i += 2 * group;
			loopcnt += group;
			r *= vsq * ln[0];
			terms *= r;
			if (asinx + terms == asinx)
				break;
			asinx += terms;
			if (group > 1)
				r *= uc;		// ajust r to last Taylor term
			}
		}

	// Revere argument reduction
	if (k > 0)
		asinx.adjustExponent(+k); // asinx*=2^k

	// Round to same precision as argument and rounding mode
	asinx.mode(x.mode());
	asinx.precision(x.precision());

	if (sign < 0)
		asinx.change_sign();
	return asinx;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		acos
//	@return 	float_precision	-	return acos(x)
//	@param      "x"	-	float_precision argument
//
//	@todo 
//
// Description:
//    Use Arccos(x)=PI/2 - Arcsin(x) or ArcCos(x)=0.5*PI-ArcSin(x)
//
float_precision acos(const float_precision& x)
	{
	size_t precision;
	float_precision y;
	const float_precision c1(1);

	if (x > c1 || x < -c1)
		{
		throw float_precision::domain_error();
		}

	precision = x.precision();
	y = _float_table(_PI, precision);
	y.adjustExponent(-1);
	y -= asin(x);

	// Round to same precision as argument and rounding mode
	y.mode(x.mode());
	y.precision(precision);
	return y;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005, 24/Jun/2022
//	@brief 		atan
//	@return		float_precision	-	return atan(x)
//	@param      "x"	-	float_precision argument
//
//	@todo 
//
// Description:
//   Use the taylot series. ArcTan(x) = x - x^3/3 + x^5/5 ...
//   With moderate use of argument reduction  using the identity. ArcTan(x)=2*ArcTan(x/(1+sqrt(1+x^2))) 
//	And Taylor terms grouping of 5 Taylor terms at a time
//	atan(0)==0
//
float_precision atan(const float_precision& x)
	{
	const int group = 5;
	size_t precision = x.precision() + 2 + (size_t)ceil(log10(x.precision()));
	intmax_t k, i;
	size_t loopcnt = 1;
	float_precision r, atanx, v(x), vsq, terms;
	const float_precision c1(1);

	if (v.iszero())
		return v;
	// Automatically calculate optimal reduction factor as a power of two
	k = 2 * (intmax_t)ceil(log(2)*log(precision));
	//k = std::min((intmax_t)30, k);  // Top of a maximum 30 reduction. only relevant if the multipier for k>2
	if (v.exponent() >= 0)
		++k;	// We only need one reduction to get x below 1
	else
		k += v.exponent(); // Avoid uncessary argument reduction if v is small 
	k = std::max((intmax_t)0, k);
	
	// Adjust the precision
	if (k > 0)
		precision += k / 4; ;
	r.precision(precision);
	atanx.precision(precision);
	v.precision(precision);
	vsq.precision(precision);
	terms.precision(precision);

	// Argument reduction. Transform the solution to ArcTan(x)=2*ArcTan(x/(1+sqrt(1+x^2)))
	for (i = k; i>0; --i)
		v = v / (c1 + sqrt(c1 + v.square()));

	vsq = v.square();
	r = v;
	atanx = v;
	if (group == 1)
		{
		// Now iterate using taylor expansion
		for (i = 3;; i += 2, ++loopcnt)
			{
			v *= vsq;
			v.change_sign();
			r = v / float_precision(i);
			if (atanx + r == atanx)
				break;
			atanx += r;
			}
		}
	else
		{
		std::vector<float_precision> vn(group);  // vn[0] is not used
		std::vector<float_precision> cn(group + 1);

		for (i = 0; i < group; ++i)
			{
			cn[i].precision(precision); vn[i].precision(precision);
			if (i == 1) vn[1] = vsq;
			if (i > 1) vn[i] = vn[i - 1] * vsq;
			}
		cn[group].precision(precision);
		// Now iterate 
		for (i = 3;; )
			{
			// Recalulate the coefficients
			intmax_t j, m;
			for (j = 0, cn[group] = c1; j < group; ++j)
				{
				cn[j] = c1;
				cn[group] *= float_precision(i + 2 * j);
				for (m = 0; m < group; ++m)
					{
					if (m == j)  continue;
					cn[j] *= float_precision(i + 2 * m);
					}
				if ((i + 2 * j) / 2 & 0x1)
					cn[j].change_sign();
				}

			cn[group] = cn[group].inverse();
			// Adding from smalles to largest number
			terms = vn[group - 1] * cn[group - 1];
			for (j = group - 1; j >= 2; --j)
				terms += cn[j - 1] * vn[j - 1];
			terms += cn[0];
			i += 2 * group;
			loopcnt += group;
			r *= vsq;
			terms *= r * cn[group];
			if (atanx + terms == atanx)
				break;
			atanx += terms;
			if (group > 1)
				r *= vn[group - 1];		// ajust r to last Taylor term
			}
		}

	atanx.adjustExponent(k);  // multiply with 2^k

	// Round to same precision as argument and rounding mode
	atanx.mode(x.mode());
	atanx.precision(x.precision());
	return atanx;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		atan2
//	@return 	float_precision	-	return the angle (in radians) from the X axis to a point (y,x).
// @param		"y"   -  float_precision y-axis
//	@param      "x"	-	float_precision x-axis
//
//	@todo 
//
// Description:
//   use atan() to calculate atan2()
//
float_precision atan2( const float_precision& y, const float_precision& x )
   {
   const size_t precision=x.precision()+2;
   float_precision atan2x;
   const float_precision c0(0);

   if( x.iszero() && y.iszero() )
      return x;

   atan2x.precision( precision );
   if( x.iszero() )
      {
      atan2x = _float_table( _PI, precision );
	  atan2x.adjustExponent(-1);		// atanx *= c05;
	  if (y < c0)
		  atan2x.change_sign();
      }
   else
      if( y.iszero() )
         {
         if( x < c0 )
            atan2x = _float_table( _PI, precision );
         else
            atan2x = c0;
         }
      else
         {
         atan2x = atan( y / x );
         if( x < c0  && y < c0 )
            atan2x -= _float_table( _PI, precision );

         if( x < c0 && y >= c0 )
            atan2x += _float_table( _PI, precision );
         }

	// Round to same precision as argument and rounding mode
   atan2x.mode( x.mode() );
   atan2x.precision( x.precision() );  

   return atan2x;
   }

 
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005 & 11-Jun-2022
//	@brief 		sin
//	@return 	float_precision	-	return sin(x)
//	@param      "x"	-	float_precision argument
//
//	@todo 
//
// Description:
//   Use the taylor series. Sin(x) = x - x^3/3! + x^5/5! ...
//   1) However first reduce x to between 0..2*PI 
//   2) Then reduced further to between 0..PI using sin(x+PI)=-Sin(x) or an extra reduction of argument
//   3) Finally reduced it to a small number using 8*ceil(log(2)*log(precision)) as reduction factor,
//		using the trisection identity: sin(3x)=3*sin(x)-4*sin(x)^3
//   4) Then Do the taylor using a coefficient scaling of 5 Taylor terms at a time
//	 5) Then reverse the reduction factor
//   The argument reduction is used to reduced the number of Taylor iterations 
//   and to minimize round off erros and calculation time
//	sin(0)==0
//
float_precision sin(const float_precision& x )
	{
	const int group = 5;
	size_t precision = x.precision() + 2 + (size_t)ceil(log10(x.precision()));
	intmax_t k;
	uintmax_t i;
	int sign;
	uintmax_t loopcnt = 1;
	float_precision r, sinx, v(x), vsq, terms;
	const float_precision c3(3), c4(4);

	if (v.iszero())
		return v;
	// Check for augument reduction and increase precision if necessary
	// Automatically calculate optimal reduction factor as a power of two
	k = 8 * (intmax_t)ceil(log(2)*log(precision));
	
	// Now use the trisection identity sin(3x)=sin(x)(3+4Sin^2(x))
	// until argument has been reduced 2/3*k times. Converting power of 2 to power of 3.
	k = (intmax_t)ceil(2.0*k / 3);
	precision += k / 4;
	r.precision(precision);
	sinx.precision(precision);
	v.precision(precision);
	vsq.precision(precision);
	terms.precision(precision);

	sign = v.sign();
	if (sign < 0)
		v.change_sign();

	// Check that argument is larger than 2*PI and reduce it if needed to the range [0..2*PI]. 
	// No need for high perecision. we just need to figure out if we need to Calculate PI with a higher precision
	if (v > float_precision(2 * 3.14159265))
		{
		// Reduce argument to between 0..2PI
		sinx = _float_table(_PI, precision);
		sinx.adjustExponent(+1);  // same as sinx*= c2;
		if (abs(v) > sinx)
			{
			r = v / sinx;
			(void)modf(r, &r);
			v -= r * sinx;
			}
		if (v < float_precision(0))
			v += sinx;
		}

	// Reduced it further to between 0..PI
	// However avoid calculating PI is not needed.
	// No need for high perecision. we just need to figure out if we need to Calculate PI with a higher precision
	if (v > float_precision(3.14159265))
		{
		if (sinx.iszero())  // PI not call before. .
			{	// Then increase redcution factor with one reducing sinx from interval [3.14..6.28] to max [1.05..2.10]
			++k;
			}
		else
			{	// We dont need to worry that we called it a second time since it will be cached from the first calculation
			sinx = _float_table(_PI, precision); 
			if (v > sinx)
				{
				v -= sinx;
				sign *= -1;  // Change sign
				}
			}
		}

	// Adjust k for final value of v when v is small (less than 1). we know it is in the interval between [0..PI]
	// This indicate that the exponent is in the range [-inf..1]
	// Avoid uncessary argument reduction if v is small 
	k += v.exponent();
	k = std::max((intmax_t)0, k);

	// Now use the trisection identity sin(3x)=3*sin(x)-4*sin(x)^3
	// k times k.  Where k is the number of reduction factor based on the needed precision of the argument.
	r = c3;
	r = pow(r, float_precision(k));  // Since r and k is an integer this wil be faster
	v /= r;
	vsq = v.square();
	r = v;
	sinx = v;

	if (group == 1)
		{
		// Now iterate using taylor expansion
		for (i = 3; ; i += 2, ++loopcnt)
			{
			r *= vsq / float_precision(i*(i - 1));
			r.change_sign();
			if (sinx + r == sinx)
				break;
			sinx += r;
			}
		}
	else
		{ 
		std::vector<float_precision> vn(group);  // vn[0] is not used
		std::vector<float_precision> cn(group);  // 

		for (i = 0; i < group; ++i)
			{
			cn[i].precision(precision); vn[i].precision(precision);
			if (i == 1) vn[1] = vsq;
			if (i > 1) vn[i] = vn[i - 1] * vsq;
			}
		// Now iterate 
		for (i = 3; ; )
			{
			intmax_t j;
			for (j = group - 1; j >= 0; --j)
				{
				if (j == group - 1)
					{
					cn[j] = float_precision((i + 2 * j - 1)*(i + 2 * j), precision);
					if ((i / 2 + j - 1) & 0x1)  // Odd
						cn[j].change_sign();
					}
				else
					{
					cn[j] = -cn[j + 1] * float_precision((i + 2 * j - 1)*(i + 2 * j), precision);
					}
				}

			cn[0] = abs(cn[0]).inverse();
			// Adding from smalles to largest number
			terms = vn[group - 1];
			if ((i / 2 + group - 1) & 0x1)
				terms.change_sign();
			for (j = group - 1; j >= 2; --j)
				terms += cn[j] * vn[j - 1];
			terms += cn[1];

			r *= vsq*cn[0];
			terms *= r;
			i += 2 * group;				// Update term count
			loopcnt += group;
			if (sinx + terms == sinx)	// Reach precision
				break;					// yes terminate loop
			sinx += terms;				// Add taylor terms to result
			if (group > 1)
				r *= vn[group - 1];		// ajust r to last Taylor term
			}
		}

	// Reverse argument reductions
	for (; k > 0; --k)
		sinx *= (c3 - c4 * sinx.square());

	// Round to same precision as argument and rounding mode
	sinx.mode(x.mode());
	sinx.precision(x.precision());
	if (sign < 0)
		sinx.change_sign();
	return sinx;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		24/Jun/2022
//	@brief 		cos
//	@return 	float_precision	-	return cos(x)
//	@param      "x"	-	float_precision argument
//
//	@todo 
//
// Description:
//   Use the cos(x)=sqrt(1+sin(x)^2)
//	This is faster than a Taylor series of cos(x) since sin(x) usually is faster and more accurate thatn the corresponding Taylor series of cos(x)
//
float_precision cos(const float_precision& x)
	{
	//size_t precision = x.precision() + 2 + (size_t)ceil(log10(x.precision()));
	float_precision cosx(x), sinx(x);
	const float_precision c1(1);
	double d;

	sinx = sin(x);
	cosx = sqrt(c1 - sinx.square());
	d = x; d = cos(d);
	if( d < 0 )
		cosx = -cosx;
	cosx.mode(x.mode());
	cosx.precision(x.precision());
	return cosx;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005, 14/Jun/2022
//	@brief 		tan
//	@return 	float_precision	-	return tan(x)
//	@param      "x"	-	float_precision argument
//
//	@todo 
//
// Description:
//   Use the identity tan(x)=Sin(x)/Sqrt(1-Sin(x)^2)
//   However first reduce x to between 0..2*PI so we can checkfor domain error
//	tan(0)==0
//   
//
float_precision tan( const float_precision& x )
   {
   const size_t precision=x.precision() + 2;
   float_precision tanx, r, v(x), pi;
   const float_precision c1(1), c3(3);

   if (v.iszero())
	   return v;
	// Increase working precision
   tanx.precision( precision );
   v.precision( precision );
   pi.precision( precision );
  
   // Check that argument is larger than 2*PI and reduce it if needed. 
   pi = _float_table( _PI, precision );
   tanx = pi;
   tanx.adjustExponent(+1);		// 2*PI
   if( abs( v ) > tanx )
      {
      r = v / tanx; 
      (void)modf( r, &r ); 
      v -= r * tanx;
      }
   if( v < float_precision( 0 ) )
      v += tanx;
    
   pi.adjustExponent(-1);	// pi *= 0.5;
   if( v == pi || v ==  pi * c3 )
      { throw float_precision::domain_error(); }

   tanx = sin( v ); 
   if( v < pi || v > pi * c3 ) 
      tanx /= sqrt( c1 - tanx.square() );
   else
      tanx /= -sqrt( c1 - tanx.square() );
   
   // Round to same precision as argument and rounding mode
   tanx.mode( x.mode() );
   tanx.precision( x.precision() );  

   return tanx;
   }

///////////////////////////////////////
//
// END TRIGONOMETRIC FUNCTIONS
//
///////////////////////////////////////

///////////////////////////////////////
//
// HYPERBOLIC FUNCTIONS
// 
//	sinh()
//	cosh()
//	tanh()
//	arcsinh()
//	arccosh()
//	arctanh()
//
///////////////////////////////////////

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/21/2013, 10/May/2022, 14/Jul/2022
//	@brief 		Calculate Sinh(x)
//	@return 	float_precision -	Return Sinh(x)
//	@param      "x"	-	   The argument
//
//	@todo  
//
// Description:
//   Use a taylor series until their is no more change in the result
//   sinh(x) == x + x^3/3!+x^5/5!+....
//   Use argument reduction via sinh(3x)=sinh(x)(3+4sinh^2(x))
//	and coefficient scalling (grouping of 5 terms has been added)
// sinh(0)==0
//
float_precision sinh(const float_precision& x)
	{
	const int group = 5;
	size_t precision = x.precision() + 2 + (size_t)ceil(log10(x.precision()));
	size_t loopcnt = 2;
	intmax_t k;
	uintmax_t i;
	float_precision r, sinhx, v(x), vsq, terms;
	const float_precision c1(1), c3(3), c4(4);

	if (x.iszero())
		return v;
	if (x.sign() < 0)
		v.change_sign();	// sinh(-x)==-sinh(x)

	// Automatically calculate optimal reduction factor as a power of two
	k = 8 * (intmax_t)ceil(log(2)*log(precision));
	k += v.exponent() + 2;
	k = std::max((intmax_t)0, k);
	// Now use the trisection identity sinh(3x)=sinh(x)(3+4Sinh^2(x))
	// until argument has been reduced 2/3*k times. Converting power of 2 to power of 3.
	k = (intmax_t)ceil(2.0*k / 3);

	// Adjust the precision
	precision += k / 4;
	v.precision(precision);
	r.precision(precision);
	sinhx.precision(precision);
	vsq.precision(precision);
	terms.precision(precision);
	r = c3;
	r = pow(r, float_precision(k));  // Since r and k is an integer this wil be faster
	v /= r;
	vsq = v.square();
	r = v;
	sinhx = v;

	if (group == 1)
		{
		// Now iterate using taylor expansion
		for (i = 3;; i += 2, ++loopcnt)
			{
			r *= vsq / float_precision(i*(i - 1));
			if (sinhx + r == sinhx)
				break;
			sinhx += r;
			}
		}
	else
		{
		std::vector<float_precision> vn(group);  // vn[0] is not used
		std::vector<float_precision> cn(group);

		for (i = 0; i < group; ++i)
			{
			cn[i].precision(precision); vn[i].precision(precision);
			if (i == 1) vn[i] = vsq;
			if (i > 1) vn[i] = vn[i - 1] * vsq;
			}
		// Now iterate 
		for (i = 3; ; )
			{
			int j;
			for (j = group - 1; j >= 0; --j)
				{
				if (j == group - 1)
					cn[j] = float_precision((i + 2 * j - 1)*(i + 2 * j));
				else
					cn[j] = cn[j + 1] * float_precision((i + 2 * j - 1)*(i + 2 * j));
				}
			// Adding from smalles to largest number
			terms = vn[group - 1];
			for (j = group - 1; j >= 2; --j)
				terms += cn[j] * vn[j - 1];
			terms += cn[1];
			r *= vsq / cn[0];
			terms *= r;
			i += 2 * group;				// Update term count
			loopcnt += group;
			if (sinhx + terms == sinhx)	// Reach precision
				break;					// yes terminate loop
			sinhx += terms;				// Add taylor terms to result
			if (group > 1)
				r *= vn[group - 1];		// ajust r to last Taylor term
			}
		}

	// Reverse argument reduction
	for (; k > 0; k--)
		sinhx *= (c3 + c4*sinhx.square());

	// Round to same precision as argument and rounding mode
	sinhx.mode(x.mode());
	sinhx.precision(x.precision());
	if (x.sign() < 0)
		sinhx.change_sign();
	return sinhx;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/21/2013, 12/May/2022, 14/Jul/2022
//	@brief 		Calculate Cosh(x)
//	@return 	float_precision -	Return Cosh()
//	@param      "x"	-	   The argument
//
//	@todo  
//
// Description:
//   Use a taylor series until their is no more change in the result
//   cosh(x) == 1 + x^2/2!+x^4/4!+....
//   Use argument reduction via cosh(3x)=cosh(x)(4cosh^2(x)-3)	
//   and coefficient rescalling (grouping of 5 terms has been added)
//	 cosh(0)==1
//
float_precision cosh(const float_precision& x)
	{
	const int group = 5;
	size_t precision = x.precision() + 2 + (size_t)ceil(log10(x.precision()));
	size_t loopcnt = 1;
	intmax_t k;
	uintmax_t i;
	float_precision r, coshx, v(x), vsq, terms;
	const float_precision c1(1), c3(3), c4(4);

	if (x.iszero())
		return v=c1;		// Preserve the precision
	if (x.sign() < 0)
		v.change_sign();  // cosh(-x) = cosh(x)

	// Automatically calculate optimal reduction factor as a power of two
	k = 8 * (intmax_t)ceil(log(2)*log(precision));
	k += v.exponent() + 2;
	// Now use the trisection identity cosh(3x)=cosh(x)(4Cosh^2(x)-3)
	// until argument has been reduced 2/3*k times. Converting power of 2 to power of 3.
	k = (intmax_t)ceil(2.0*k / 3);

	// Adjust the precision
	precision += k;
	r.precision(precision);
	coshx.precision(precision);
	v.precision(precision);
	vsq.precision(precision);
	terms.precision(precision);

	r = c3;
	r = pow(r, float_precision(k));  // Since r and k is an integer this wil be faster
	v /= r;
	vsq = v.square();
	r = c1;
	coshx = r;

	if (group == 1)
		{
		// Now iterate using taylor expansion
		for (i = 2;; i += 2, ++loopcnt)
			{
			r *= vsq / float_precision(i*(i - 1));
			if (coshx + r == coshx)
				break;
			coshx += r;
			}	
		}
	else
		{
		std::vector<float_precision> vn(group);  // vn[0] is not used
		std::vector<float_precision> cn(group);

		for (i = 0; i < group; ++i)
			{
			cn[i].precision(precision); vn[i].precision(precision);
			if (i == 1) vn[i] = vsq;
			if (i > 1) vn[i] = vn[i - 1] * vsq;
			}
		// Now iterate 
		for (i = 2; ; )
			{
			intmax_t j;
			for (j = group - 1; j >= 0; --j)
				{
				if (j == group - 1)
					cn[j] = float_precision((i + 2 * j - 1)*(i + 2 * j));
				else
					cn[j] = cn[j + 1] * float_precision((i + 2 * j - 1)*(i + 2 * j));
				}

			// Adding from smalles to largest number
			terms = vn[group - 1];
			for (j = group - 1; j >= 2; --j)
				terms += cn[j] * vn[j - 1];
			terms += cn[1];
			r *= vsq / cn[0];
			terms *= r;
			i += 2 * group;				// Update term count
			loopcnt += group;
			if (coshx + terms == coshx)	// Reach precision
				break;					// yes terminate loop
			coshx += terms;				// Add taylor terms to result
			if (group > 1)
				r *= vn[group - 1];		// ajust r to last Taylor term
			}
		}

	// Reverse argument reduction
	for (; k > 0; k--)
		coshx *= (c4 * coshx.square() - c3);

	// Round to same precision as argument and rounding mode
	coshx.mode(x.mode());
	coshx.precision(x.precision());
	return coshx;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/21/2013
//	@brief 		Calculate Tanh(x)
//	@return 	float_precision -	Return Tanh()
//	@param      "x"	-	   The argument
//
//	@todo  
//
// Description:
//	tanh = ( exp(x) - exp(-x) ) / ( exp( x) + exp(-x) )=(e^(2x)-1/(e^(2x)+1)
//	tanh(0)==0
// 
//
float_precision tanh( const float_precision& x )
   {
   float_precision v(x), vsq;
   const float_precision c1(1);

   v.precision( x.precision() + 1 );
   vsq.precision( x.precision() + 1 );
   v = exp( v );
   vsq= v.square();
   v = (vsq-c1)/(vsq+c1);

   // Round to same precision as argument and rounding mode
   v.mode( x.mode() );
   v.precision( x.precision() );  
   return v;
   }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/25/2013
//	@brief 		Calculate ArcSinh(x)
//	@return 	float_precision -	Return ArcSinh()
//	@param      "x"	-	   The argument
//
//	@todo  
//
// Description:
//	ArcSinh=Ln(x+Sqrt(x^2+1))
//	asinh(0)==0
// 
float_precision asinh( const float_precision& x )
   {
   float_precision v(x);
   const float_precision c1(1);

   if (v.iszero())
	   return v;
   v.precision( x.precision() + 1 );
   v = log(v+sqrt(v.square()+c1));
   
   // Round to same precision as argument and rounding mode
   v.mode( x.mode() );
   v.precision( x.precision() );  
   return v;
   }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/25/2013
//	@brief 		Calculate ArcCosh(x)
//	@return 	float_precision -	Return ArcCosh()
//	@param      "x"	-	   The argument
//
//	@todo  
//
// Description:
//	ArcCosh=Ln(x+Sqrt(x^2-1))
// 
//
float_precision acosh( const float_precision& x )
   {
   float_precision v(x);
   const float_precision c1(1);

   if( x < c1 )
      { throw float_precision::domain_error(); }
   
   v.precision( x.precision() + 1 );
   v = log(v+sqrt(v.square()-c1));
   
   // Round to same precision as argument and rounding mode
   v.mode( x.mode() );
   v.precision( x.precision() );  
   return v;
   }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/25/2013
//	@brief 		Calculate ArcTanh(x)
//	@return 	float_precision -	Return ArcTanh()
//	@param      "x"	-	   The argument
//
//	@todo  
//
// Description:
//	ArcTanh=0.5*Ln((1+x)/(1-x))
//	Atanh(0)==0
// 
//
float_precision atanh( const float_precision& x )
   {
   float_precision v(x);
   const float_precision c1(1);

   if( x >= c1 || x <= -c1 )
      { throw float_precision::domain_error(); }
   if (x.iszero())
	   return v;
   v.precision( x.precision() + 1 );
   v = log((c1+v)/(c1-v));
   v.adjustExponent(-1);   // v *= c05;
   
   // Round to same precision as argument and rounding mode
   v.mode( x.mode() );
   v.precision( x.precision() );  
   return v;
   }

///////////////////////////////////////
//
// END HYPERBOLIC FUNCTIONS
//
///////////////////////////////////////


///////////////////////////////////////
//
// FAST integer division and remaining using Floating point arithmetic
//
///////////////////////////////////////

#if false

int_precision _int_precision_fastdiv( const int_precision &s1, const int_precision &s2 )
	{
	size_t ss;
	int_precision r2;
	float_precision f1, f2, rf;

	ss = std::max(s1.size(), s2.size());
	ss = (int)(ceil(Bitsiptype*ss / log2(BASE_10)));
	f1.precision( ss+2);
	f2.precision( ss+2);
	rf.precision( ss+2 );
	f1=float_precision(s1, ss+ 2);
	f2=float_precision( s2, ss+ 2);
	rf=f1/f2;
	r2 = (int_precision)rf; 
	return r2;
	}

int_precision _int_precision_fastrem( const int_precision &s1, const int_precision &s2 )
	{
	size_t ss;
	int_precision r2;
	float_precision f1, f2, rf;

	ss = std::max(s1.size(), s2.size());
	ss = (int)(ceil(Bitsiptype*ss / log2(BASE_10)));
	f1.precision( ss+2);
	f2.precision( ss+2);
	rf.precision( ss+2 );
	f1=float_precision(s1, ss+ 2);
	f2=float_precision( s2, ss+ 2);
	rf=f1/f2;
	r2 = (int_precision)rf; 
	r2=s1-s2*r2;
	return r2;
	}
#endif

////////////////////////////////////////
//
//
// Auxillary functions:
//		AGM, nrooth
//
////////////////////////////////////////

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		20/May/2022
//	@brief 		Calculate agm(x,y)
//	@return 	float_precision -	Return agm(x,y)
//	@param      "x"	-	   The argument
//	@param      "y"	-	   The argument
//
//	@todo  
//
// Description:
//	Calculate the Arithmetic-Geometric mean, through iteration
// 
//
static float_precision AGM(const float_precision& a, const float_precision& b)
	{
	const int guard = 0;
	size_t precision = std::max(a.precision(), b.precision()) + guard;
	size_t loopcnt;
	float_precision x(a), y(b), xnew(a), ynew(b);
	float_precision diff;

	x.precision(precision);
	y.precision(precision);
	xnew.precision(precision);
	ynew.precision(precision);
	for (loopcnt = 1;; ++loopcnt)
		{
		xnew = 0.5*(x + y);
		ynew = sqrt(x*y);
		diff = xnew - ynew;
		if (diff.iszero() || xnew == x || ynew == y)
			break;
		x = xnew;
		y = ynew;
		}

	xnew.precision(precision - guard);
	return xnew;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		12/Nov/2015, 24/July/2022
//	@brief 		Calculate n root of x  (x^(1/n)
//	@return 	float_precision -	Return nroot(a)
//	@param      "x"	-	The nroot argument
//
//	@todo  
//
// Description:
//   nroot(V)
//   The nth root of x^(1/n) No Equivalent standard C function call
//   Separate exponent. e.g. nroot(V*10^x)=10^x/2*nroot(V)
//   Un=U(1/n)(n+1-VU^n)
//   Then Un == 1/nroot(V). and nroot(V) = 1/Un
// The function has been improved using Newton with iterative deepening creating
//	a speed up with a factor of 3 over the classic Newton method.
// This is a much much faster option instead of the traditional pow() function x^y
// and that is why it has been added as a separate function
//
float_precision nroot(const float_precision& a, const uintmax_t n)
	{
	const size_t extra = 2;
	const size_t precision = a.precision();
	const eptype expo = a.exponent();
	const float_precision c1(1);
	eptype expo_sq;
	size_t digits;
	double fv, fv2;
	float_precision r, x, y(a), fn(n);

	if (a.iszero() || a == c1 || n == 1)
		return a;
	if (a.sign() < 0)
		throw float_precision::domain_error();

	y.precision(precision + extra);
	expo_sq = expo / (eptype)n;
	y.exponent(expo - (eptype)n * expo_sq);
	// Do iteration using guard digits higher precision
	x.precision(precision + extra);
	fn.precision(precision + extra);

	// Get a initial guess using ordinary floating point
	fv = y;
	/*	if (expo - 2 * expo_sq > 0)
	fv *= 2.0;
	else
	if (expo - 2 * expo_sq < 0)
	fv *= 0.5;*/
	fv = pow(fv, 1.0 / n);  // set the initial guess with at approx 16 correct digits
	fv = 1 / fv;
	x = float_precision(fv);
	fn = 1 / fn;
	// Now iterate using Netwon  x=x*(-yx^n+(n+1))/n
	for (digits = std::min((size_t)32, precision); ; digits = std::min(precision + extra, digits * 2))
		{
		// Increase precision by a factor of two
		r.precision(digits);
		x.precision(digits);
		float_precision p(x);
		float_precision res(1, digits);
		// Do x^n
		for (uintmax_t i = n; i > 0; i >>= 1)
			{
			if ((i & 0x1) != 0)
				res *= p;  // Odd
			if (i>1)
				p *= p;
			}
		// Notice y is the original number to nroot which has the full precision 
		r = float_precision(n + 1) - y*res; // (n+1)-yx^n
		r *= fn;							// (-yx^n+(n+1))/n
		x *= r;								// x=x*(-yx^n+(n+1))/n
		if (digits == precision + extra)	 // Reach final iteration step in regards to precision
			{
			r.precision(precision + 1);			// round to final precision
			if (r == c1)		// break if no improvement
				break;
			}
		}

	fv2 = x;
	x = 1 / x;			// n root of x is now 1/x;
	x.exponent(x.exponent() + expo_sq);
	// Round to same precision as argument and rounding mode
	x.mode(a.mode());
	x.precision(precision);
	return x;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		12/Nov/2015
//	@brief 		Calculate n root of x  (x^(1/n)
//	@return 	float_precision -	Return nroot(a)
//	@param      "x"	-	The nroot argument
//
//	@todo  
//
// Description:
//   nroot(V)
//   The nth root of x^(1/n) No Equivalent standard C function call
//   Separate exponent. e.g. nroot(V*10^x)=10^x/2*nroot(V)
//   Un=U(1/n)(n+1-VU^n)
//   Then Un == 1/nroot(V). and nroot(V) = 1/Un
// The function has been improved using Newton with iterative deepening creating
//	a speed up with a factor of 3 over the classic Newton method.
// This is a much much faster option instead of the traditional pow() function x^y
// and that is why it has been added as a separate function
//
/*
float_precision nroot(const float_precision& x, uintmax_t n)
	{
	const size_t extra = 2;
	size_t precision;
	size_t digits;
	eptype expo, expo_sq;
	double fv;
	float_precision r, u, v, fn(n);
	const float_precision c1(1);

	if (x.iszero() || x == c1 || n == 1)
		return x;
	if (x.sign() < 0)
		{
		throw float_precision::domain_error();
		}
	precision = x.precision();
	v.precision(precision + extra);
	v = x;
	expo = v.exponent();
	expo_sq = expo / 2;
	v.exponent(expo - 2 * expo_sq);
	r.precision(precision + extra); // Do iteration using 2 digits higher precision
	u.precision(precision + extra);
	fn.precision(precision + extra);

	// Get a initial guess using ordinary floating point
	fv = v;
	fv = pow(fv, 1.0 / n);  // set the initial guess with at approx 16 correct digits
	fv = 1 / fv;

	u = float_precision(fv);
	fn = 1 / fn;
	// Now iterate using Netwon  Un=U*(-VU^n+(n+1))/n
	for (digits = std::min((size_t)32, precision); ; digits = std::min(precision + extra, digits * 2))
		{
		// Increase precision by a factor of two
		r.precision(digits);
		u.precision(digits);
		float_precision p(u);
		float_precision res(1, digits);
		// DO U^N
		for (uintmax_t i = n; i > 0; i >>= 1)
			{
			if ((i & 0x1) != 0)
				res *= p;  // Odd
			if (i>1)
				p *= p;
			}
		// Notice V is the original number to nroot which has the full precision 
		// so we start by assigning it to r, rounding it to the precision of r
		r = v*res;						// VU^n
		r = float_precision(n + 1) - r; // (n+1)-VU^n
		r *= fn;						// (-VU^n+(n+1))/n
		u *= r;							// U=U*(-VU^n+(n+1))/n
		if (digits == precision + extra) // Reach final iteration step in regards to precision
			{
			r.precision(precision + 1);	// round to final precision
			if (r == c1)		// break if no improvement
				break;
			}
		}

	u = 1 / u;			// n root of u is now 1/u;
	u.exponent(u.exponent() + expo_sq);
	// Round to same precision as argument and rounding mode
	u.mode(x.mode());
	u.precision(precision);

	return u;
	}
*/

///////////////////////////////////////
//
// FLOATING POINT FUNCTIONS
//
///////////////////////////////////////