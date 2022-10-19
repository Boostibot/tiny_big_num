#ifndef INC_INTERVAL
#define INC_INTERVAL

/*
 *******************************************************************************
 *
 *
 *                       Copyright (c) 2002-2022
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
 * Module name     :   intervalprecision.h
 * Module ID Nbr   :
 * Description     :   Interval Precision arithmetic template class
 *                     Works with the int_precision and float_precision classes
 * --------------------------------------------------------------------------
 * Change Record   :
 *
 * Version	Author/Date		Description of changes
 * -------  -----------		----------------------
 * 01.01	HVE/020209		Initial release
 * 01.02    HVE/030421		Optimized the * operator to reduce the 8 multiplications to 2.
 * 01.03	HVE/JUN-26-2014	Added is_empty(), contains_zero() method to the class
 * 01.04	HVE/JUN-27-2014	Corrected several errors in in cin >> template function
 * 01.05	HVE/JUN-28-2014	Added is_class() method for getting the interval classification
 *							and width() method for the interval width
 * 01.06	HVE/JUN-30-2014	An error was corrected for interval subtraction of float_preicsion numbers
 *							Also added the method bool contain() for test if a float or interval is included in the interval
 * 01.07	HVE/JUL-6-2014	Corrected an error in /= for the software emulation of of float & double
 * 01.08	HVE/JUL-13-2014	Added Hardware support for interval arithmetic when applicable. Also fix several errors in the
 *							implementation of sqrt, log, log10, exp and pow functions. Also added new method is_class(), is_empty()
 * 01.09	HVE/JUL-15-2014	Added support for Sin(), Cos() and Tan() interval trigonometric functions.
 * 01.10	HVE/JUL-17-2014	Added support for atan() interval trigonometric function
 * 01.11	HVE/JUL-22-2014 Found a bug that floating point was not reset to near (default by IEEE754) after a hardware supported multiplication
 * 01.12	HVE/JUL-22-2014	Added support for asin() interval trigonometric function
 * 01.13	HVE/JUL-29-2014	Added support for interval versions of LN2, LN10 and PI
 * 01.14	HVE/AUG-10-2014	Added support for mixed mode arithmetic for interval +,- classes
 * 01.15	HVE/JUN-20-2015	Fixed and un-declare variable x when compiling with no interval hardware support
 * 01.16	HVE/Jul-07-2019	Moved Hardware support up prior to the template class definition to make the code more portable and added <iostream> header
 * 01.17    HVE/Jul-07-2019 Make the code more portable to a GCC environment
 * 01.18	HVE/24-Mar-2021 Updated license info
 * 01.19	HVE/4-Jul-2021	Added software interval runding via towsum and twoproduct functions and other functions
 * 01.20	HVE/5-Jul-2021	Replace deprecreated headers with current headers
 * 01.21	HVE/15-Jul-2021 Decpreated hardware support for interval arithmetic since this was not portable and didnt takes advantages of the latest 
 *							Intel instructions set. instead if used only software emaulation of intervals.
 * 01.22	HVE/29-Jul-2021 Corrected bugs in all the trigonometir functions and added interval version of hyperbolic functions
 *							sinh(), cosh(), tanh(), asinh(), acosh(), atanh().
 * 01.23	HVE/30-Jul-2021 Added intervalsection(), unionsection(), boolean precedes(), interior()
 *
 * End of Change Record
 * --------------------------------------------------------------------------
*/


/* define version string */
static char _VinterP_[] = "@(#)intervalprecision.h 01.23 -- Copyright (C) Henrik Vestermark";

#include <cfloat>
#include <cmath>
#include <algorithm>
#include <iostream>


/// The four different interval classification
/// # ZERO			a=0 && b=0
/// # POSITIVE		a>=0 && b>0
/// # NEGATIVE		a<0 && b<=0
/// # MIXED			a<0 && b>0
enum int_class { NO_CLASS, ZERO, POSITIVE, NEGATIVE, MIXED };

//
// Interval class
// Realisticly the class Type can be float, double, int_precision, & float_precision
// Since float and double are done unsing the Intel cpu (H/W) and using "specilization"
// the int_precision and float_precision is done using the arbitrary precision packages
// Since their is no way to specific portable ROUND_DOWN and ROUND_UP mode the template class
// heavily use specilization. For any other type besides float, double, int_precision and float_precision
// the operations is not defined
//
template<class _IT> class interval {
   _IT low, high;
   public:
      typedef _IT value_type;

      // constructor. zero, one or two arguments for type _IT
      interval()							{ low = _IT(0); high = _IT(0); }
	  interval( const _IT& d )				{ low = _IT(d); high = _IT(d); }
	  interval( const _IT& l, const _IT& h) { if( l < h ) { low =l; high = h; } else { low = h; high = l; } }
	  // Constrcutor for mixed type _IT != _X (base types). Allows auto construction of e.g. interval<float_precision> x(float)
	 template <class _X> interval(const _X& x) { low = _IT(x); high = _IT(x); }

      // constructor for any other type to _IT. Both up and down conversion possible
      template<class X> interval( const interval<X>& a ) : low((_IT)a.lower()), high((_IT)a.upper()) 
	  /*: low(_IT(a.lower())), high( _IT(a.upper())) */
	  {
		/*  int is = sizeof(_IT), js = sizeof(X);
		  cout << "_IT=" << typeid(_IT).name() << " X=" << typeid(X).name() << endl;
		  fpdown(); low = _IT(a.lower());
		  fpup(); high = _IT(a.upper());
		  fpnear();
		  if (typeid(_IT) == typeid(double) && typeid(X)==typeid(float))  // upscaling float to double
			  {
			  double u = high;
			  u += u * 0.5f * FLT_EPSILON;interval<float_precision>(ip1.lower(),ip1.upper())
			 // high = _IT(u);
			  }
		
		  if (a.lower() < a.upper()) 
			{ low = _IT(a.lower());
			high = _IT(a.upper()); }
		  else { low = _IT(a.upper()); high = _IT(a.lower()); }
		  */
	  }

      // Coordinate functions
      _IT upper() const					{ return high; }
      _IT lower() const					{ return low; }
      _IT upper( const _IT& u )			{ return ( high = u ); }
      _IT lower( const _IT& l )			{ return ( low = l ); }

      _IT center() const				{ return ( high + low ) / _IT(2); }
      _IT radius() const				{ _IT r; r =( high - low ) / _IT(2); if( r < _IT(0) ) r = -r; return r; }
	  _IT width() const					{ _IT r; r = high - low; if (r < _IT(0)) r = -r; return r; }

	  bool contain_zero() const			{ return low <= _IT(0) && _IT(0) <= high; }  // Obsolete. use contains() instead.
	  bool contain( const _IT& f=_IT(0)){ return low <= f && f <= high;  }
	  bool contain(const interval<_IT>& i) { return low <= i.lower() && i.upper() <= high; }
	  bool is_empty() const				{ return high < low; }

	  enum int_class is_class() const	{
										if (low == _IT(0) && high == _IT(0)) return ZERO;
										if (low >= _IT(0) && high > _IT(0)) return POSITIVE;
										if (low < _IT(0) && high <= _IT(0)) return NEGATIVE;
										if (low < _IT(0) && high > _IT(0)) return MIXED;
										return NO_CLASS;
										}
	  
	  // Conversion methods. Safer and less ambiguios than overloading implicit/explivit conversion operators
//	  std::string toString() const		{ return _float_precision_ftoa(this); }

	  // Operators
	  operator short() const			{ return (short)((high + low) / _IT(2)); }				// Conversion to short
	  operator int() const				{ return (int)( ( high + low ) / _IT(2) ); }			// Conversion to int
	  operator long() const				{ return (long)((high + low) / _IT(2)); }				// Conversion to long
	  operator unsigned short() const	{ return (unsigned short)((high + low) / _IT(2)); }		// Conversion to unsigned short
	  operator unsigned int() const		{ return (unsigned int)((high + low) / _IT(2)); }		// Conversion to unsigned int
	  operator unsigned long() const	{ return (unsigned long)((high + low) / _IT(2)); }		// Conversion to unsigned long
	  operator double() const			{ return (double)( ( high + low ) / _IT(2) ); }			// Conversion to double
	  operator float() const			{ return high == low? (float)low : (float)((high + low) / _IT(2)); }				// Conversion to float
	  operator int_precision() const	{ return int_precision(((high + low) / _IT(2))); }		// Conversion to int_precision
	  operator float_precision() const	{ return (float_precision)((high + low) / _IT(2)); }	// Conversion to float_precision

      _IT *ref_lower()					{ return &low; }
      _IT *ref_upper()					{ return &high; }

      // Essential operators
      interval<_IT>& operator= ( const interval<_IT>& );
      interval<_IT>& operator+=( const interval<_IT>& );
      interval<_IT>& operator-=( const interval<_IT>& );
      interval<_IT>& operator*=( const interval<_IT>& );
      interval<_IT>& operator/=( const interval<_IT>& );
	  interval<_IT>& operator&=( const interval<_IT>& );
	  interval<_IT>& operator|=( const interval<_IT>& );
	  interval<_IT>& operator^=( const interval<_IT>& );

	  // Exception class
	  class bad_int_syntax {};
	  class bad_float_syntax {};
	  class out_of_range   {};
	  class divide_by_zero {};
	  class domain_error   {};
	  class base_error		{};
   };


// Unary and Binary arithmetic
// Arithmetic + Binary and Unary
template <class _IT, class _X> inline interval<_IT> operator+( const interval<_IT>&, const _X&);
template <class _IT, class _X> inline interval<_IT> operator+( const _X&, const interval<_IT>&);
inline interval<float_precision> operator+(const interval<float_precision>&, const float_precision&);	// Specialization for interval<float_precision> and float_precision
inline interval<float_precision> operator+(const float_precision&, const interval<float_precision>&);	// Specialization for ifloat_precision and interval<float_precision>

template<class _IT> interval<_IT> operator+( const interval<_IT>&, const interval<_IT>& );
template<class _IT> interval<_IT> operator+( const interval<_IT>& );									// Unary

// Arithmetic - Binary and Unary
template <class _IT, class _X> inline interval<_IT> operator-(const interval<_IT>&, const _X&);
template <class _IT, class _X> inline interval<_IT> operator-(const _X&, const interval<_IT>&);
inline interval<float_precision> operator-(const interval<float_precision>&, const float_precision&);	// Specialization for interval<float_precision> and float_precision
inline interval<float_precision> operator-(const float_precision&, const interval<float_precision>&);	// Specialization for ifloat_precision and interval<float_precision>

template<class _IT> interval<_IT> operator-( const interval<_IT>&, const interval<_IT>& );
template<class _IT> interval<_IT> operator-( const interval<_IT>& );									// Unary

// Arithmetic * Binary
template <class _IT, class _X> inline interval<_IT> operator*(const interval<_IT>&, const _X&);
template <class _IT, class _X> inline interval<_IT> operator*(const _X&, const interval<_IT>&);
inline interval<float_precision> operator*(const interval<float_precision>&, const float_precision&);	// Specialization for interval<float_precision> and float_precision
inline interval<float_precision> operator*(const float_precision&, const interval<float_precision>&);	// Specialization for ifloat_precision and interval<float_precision>

template<class _IT> interval<_IT> operator*( const interval<_IT>&, const interval<_IT>& );

// Arithmetic / Binary
template <class _IT, class _X> inline interval<_IT> operator/(const interval<_IT>&, const _X&);
template <class _IT, class _X> inline interval<_IT> operator/(const _X&, const interval<_IT>&);
inline interval<float_precision> operator/(const interval<float_precision>&, const float_precision&);	// Specialization for interval<float_precision> and float_precision
inline interval<float_precision> operator/(const float_precision&, const interval<float_precision>&);	// Specialization for ifloat_precision and interval<float_precision>

template<class _IT> interval<_IT> operator/( const interval<_IT>&, const interval<_IT>& );

template<class _IT> interval<_IT> operator&( const interval<_IT>&, const interval<_IT>& );
template<class _IT> interval<_IT> operator|( const interval<_IT>&, const interval<_IT>& );
template<class _IT> interval<_IT> operator^( const interval<_IT>&, const interval<_IT>& );

// Boolean Comparison Operators
template<class _IT> bool operator==(const interval<_IT>&, const interval<_IT>&);
template<class _IT> bool operator!=(const interval<_IT>&, const interval<_IT>&);

// Other functions
template<class _IT> interval<_IT> abs(const interval<_IT>&);
template<class _IT> interval<_IT> intersection(const interval<_IT>&, const interval<_IT>&);
template<class _IT> interval<_IT> unionsection(const interval<_IT>&, const interval<_IT>&);
template<class _IT> bool interior(const interval<_IT>&, const interval<_IT>&);
template<class _IT> bool precedes(const interval<_IT>&, const interval<_IT>&);

// Manifest Constants like PI, LN2 and LN10
inline interval<float> int_pifloat();
inline interval<float_precision> int_pi(const unsigned int);
inline interval<float> int_ln2float();
inline interval<float_precision> int_ln2(const unsigned int);
inline interval<float> int_ln10float();
inline interval<float_precision> int_ln10(const unsigned int);

// Elementary functions
inline interval<float> sqrt(const interval<float>&);
inline interval<double> sqrt(const interval<double>&);
inline interval<float_precision> sqrt(const interval<float_precision>&);

inline interval<float> log( const interval<float>& );
inline interval<double> log(const interval<double>&);
inline interval<float_precision> log(const interval<float_precision>&);

inline interval<float> log10(const interval<float>&);
inline interval<double> log10(const interval<double>&);
inline interval<float_precision> log10(const interval<float_precision>&);

inline interval<float> exp(const interval<float>&, const float);
inline interval<double> exp(const interval<double>&, const double);
inline interval<float_precision> exp(const interval<float_precision>&);

inline interval<float> pow(const interval<float>&, const float );
inline interval<double> pow( const interval<double>&, const double );
inline interval<float_precision> pow( const interval<float_precision>& );

// Trigonometric functions
inline interval<float> sin(const interval<float>&);
inline interval<double> sin(const interval<double>&);
inline interval<float_precision> sin(const interval<float_precision>&);

inline interval<float> cos(const interval<float>&);
inline interval<double> cos(const interval<double>&);
inline interval<float_precision> cos(const interval<float_precision>&);

inline interval<float> tan(const interval<float>&);
inline interval<double> tan(const interval<double>&);
inline interval<float_precision> tan(const interval<float_precision>&);

inline interval<float> asin(const interval<float>&);
inline interval<double> asin(const interval<double>&);
inline interval<float_precision> asin(const interval<float_precision>&);

inline interval<float> acos(const interval<float>&);
inline interval<double> acos(const interval<double>&);
inline interval<float_precision> acos(const interval<float_precision>&);

inline interval<float> atan(const interval<float>&);
inline interval<double> atan(const interval<double>&);
inline interval<float_precision> atan(const interval<float_precision>&);

// Hyperbolic functions
inline interval<float> sinh(const interval<float>&);
inline interval<double> sinh(const interval<double>&);
inline interval<float_precision> sinh(const interval<float_precision>&);

inline interval<float> cosh(const interval<float>&);
inline interval<double> cosh(const interval<double>&);
inline interval<float_precision> cosh(const interval<float_precision>&);

inline interval<float> tanh(const interval<float>&);
inline interval<double> tanh(const interval<double>&);
inline interval<float_precision> tanh(const interval<float_precision>&);

inline interval<float> asinh(const interval<float>&);
inline interval<double> asinh(const interval<double>&);
inline interval<float_precision> asinh(const interval<float_precision>&);

inline interval<float> acosh(const interval<float>&);
inline interval<double> acosh(const interval<double>&);
inline interval<float_precision> acosh(const interval<float_precision>&);

inline interval<float> atanh(const interval<float>&);
inline interval<double> atanh(const interval<double>&);
inline interval<float_precision> atanh(const interval<float_precision>&);

// Low level support functions
inline float tofloat(const double &, const enum round_mode );
inline double add_down(const double &, const double &);
inline double add_up(const double &, const double &);
inline double sub_down(const double &, const double &);
inline double sub_up(const double &, const double &);
inline double mul_down(const double &, const double &);
inline double mul_up(const double &, const double &);
inline double div_down(const double &, const double &);
inline double div_up(const double &, const double &);
inline double sqrt_down(const double &);
inline double sqrt_up(const double &);

//////////////////////////////////////////////////////////////////////////////////////
///
/// Interval constants like PI, LN2 and LN10, E
///
//////////////////////////////////////////////////////////////////////////////////////

// Load manifest constant PI, LN, LN10, E, SQRT2
//
static const interval<double> PI(3.1415926535897931, 3.1415926535897936);
static const interval<double> LN2(0.69314718055994529, 0.69314718055994540);
static const interval<double> LN10(2.3025850929940455, 2.3025850929940459);
static const interval<double> E(2.7182818284590451, 2.7182818284590455);
static const interval<double> SQRT2(1.4142135623730947, 1.4142135623730951);

//
// Output Operator <<
//
template<class _Ty> inline std::ostream& operator<<( std::ostream& strm, interval<_Ty>& a ) { return strm << "[" << a.lower() << "," << a.upper() << "]"; }

// Input operator >>
//
template<class _Ty> inline std::istream& operator>>( std::istream& strm, interval<_Ty>& c )
   {
   _Ty l, u; char ch;
   if( strm >> ch && ch != '[')
      strm.putback(ch), strm >> l, u = l;
	else
      if( strm >> l >> ch && ch != ',')
	      if( ch == ']')
	         u = l;
	      else
            strm.putback( ch ); // strm.setstate(std::ios::failbit);
	   else
         if( strm >> u >> ch && ch != ']')
	         strm.putback( ch ); //, strm.setstate(ios_base::failbit);

   if(!strm.fail())
	   c = interval<_Ty>( l, u );

   return strm;
   }


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//    Essential Operators =,+=,-=,*=,/=
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////


// Assignment operator. Works for all class types
//
template<class _IT> inline interval<_IT>& interval<_IT>::operator=( const interval<_IT>& a )
   {
   low = a.lower();
   high = a.upper();
   return *this;
   }

// += operator. Works all other classes.
// Please note that this is for all integer classes. interval<int>, interval<long>, interval<int_precision>
// were there os no loss of precision
//
template<class _IT> inline interval<_IT>& interval<_IT>::operator+=( const interval<_IT>& a )
   {
   low += a.lower();
   high += a.upper();
   return *this;
   }

// Specilization for float_precision and +=
//
 template<> inline interval<float_precision>& interval<float_precision>::operator+=( const interval<float_precision>& a )
   {
   low.mode( ROUND_DOWN );
   low += a.lower();
   high.mode( ROUND_UP );
   high += a.upper();
   return *this;
   }


// Specialization for float and +=
// That can work with both managed an unmanaged application.
// HARDWARE_SUPPORT indicate if we use the build in ROUND_DOWN or ROUND_UP Mode when performing the operation
// For now Hardware support we use software emulation to do the interval arithmetic in software
//
template<> inline interval<float>& interval<float>::operator+=( const interval<float>& a )
	{
	low = tofloat(add_down(low, a.lower()), ROUND_DOWN );
	high = tofloat(add_up(high, a.upper()), ROUND_UP );
	return *this;
	}

// Specialization for double and +=
// That can work with both managed an unmanaged application.
// HARDWARE_SUPPORT indicate if we use the build in ROUND_DOWN or ROUND_UP Mode when performing the operation
// For now Hardware support we covert to float_precision to do the interval arithmetic in software
//
template<> inline interval<double>& interval<double>::operator+=( const interval<double>& a )
	{
	low = add_down(low, a.lower());
	high = add_up(high, a.upper());
	return *this;
	}

// -= operator. Works all other classes.
// Please note that this is for all integer classes. interval<int>, interval<long>, interval<int_precision>
// were there is no loss of precision
//
template<class _IT> inline interval<_IT>& interval<_IT>::operator-=( const interval<_IT>& a )
   {
   low -= a.high;
   high -= a.low;
   return *this;
   }

// Specilization for float_precision and -=
//
template<> inline interval<float_precision>& interval<float_precision>::operator-=( const interval<float_precision>& a )
   {
   low.mode( ROUND_DOWN );
   low -= a.upper();
   high.mode( ROUND_UP );
   high -= a.lower();
   return *this;
   }

// Specialization for float and -=
// That can work with both managed an unmanaged application.
// HARDWARE_SUPPORT indicate if we use the build in ROUND_DOWN or ROUND_UP Mode when performing the operation
// For now Hardware support we covert to float_precision to do the interval arithmetic in software
//
template<> inline interval<float>& interval<float>::operator-=( const interval<float>& a )
	{
	interval<float> b(-a.upper(), -a.lower());
	low = tofloat(add_down(low, b.lower()), ROUND_DOWN);
	high = tofloat(add_up(high, b.upper()), ROUND_UP);
	return *this;
	}

// Specialization for double and -=
// That can work with both managed an unmanaged application.
// HARDWARE_SUPPORT indicate if we use the build in ROUND_DOWN or ROUND_UP Mode when performing the operation
// For now Hardware support we do software emaultion to do the interval arithmetic in software
//
template<> inline interval<double>& interval<double>::operator-=( const interval<double>& a )
	{
	interval<double> b(-a.upper(), -a.lower());
	low = add_down(low, b.lower());
	high = add_up(high, b.upper());
	return *this;
	}

// Works all other classes.
// Please note that this is for all interger classes. interval<int>, interval<long>, interval<int_precision>
// were there is no loss of precision
// Instead of doing the mindless low = MIN(low*a.high, low*a.low,high*a.low,high*a.high) and
// high = MAX(low*a.high, low*a.low,high*a.low,high*a.high) requiring a total of 8 multiplication
//
//   low, high, a.low, a.high    result
//    +     +     +     +        +  +  [ low*a.low, high*a.high ]
//    +     +     -     +        -  +  [ high*a.low, high*a.high ]
//    +     +     -     -        -  -  [ high*a.low, low*a.high ]
//    -     +     +     +        -  +  [ low*a.high, high*a.high ]
//    -     +     -     +        -  +  [ MIN(low*a.high,high*a.low), MAX(low*a.low,high*a.high) ]
//    -     +     -     -        -  -  [ high*a.low, low*a.low ]
//    -     -     +     +        -  -  [ low*a.high, high,a.low ]
//    -     -     -     +        -  -  [ low*a.high, low*a.low ]
//    -     -     -     -        +  +  [ high*a.high, low * a.low ]
//
template<class _IT> inline interval<_IT>& interval<_IT>::operator*=( const interval<_IT>& a )
   {
   _IT l, h, t;

   if( low >= _IT(0) ) //
      { // both low and high >= 0
      if( a.lower() >= _IT(0) )
         { // a.low >=0, a.high >= 0
         l = low * a.lower();
         h = high * a.upper();
         }
      else
         if( a.upper() >= _IT(0) )
            {  //  a.low < 0, a.high >= 0
            l = high * a.lower();
            h = high * a.upper();
            }
         else
            { // a.low and a.high < 0
            l = high * a.lower();
            h = low * a.upper();
            }
      }
   else
      if( high >= _IT(0) )
         {  // low < 0, high >= 0
         if( a.lower() >= _IT(0) )
            { // a.low >=0, a.high >= 0
            l = low * a.upper();
            h = high * a.upper();
            }
         else
            if( a.upper() >= _IT(0) )
               {  //  a.low < 0, a.high >= 0
               l = low * a.upper(); if ( l > ( t = high * a.lower() ) ) l = t;
               h = high * a.upper(); if ( h < ( t = low * a.lower() ) ) h = t;
               }
            else
               { // a.low and a.high < 0
               l = high * a.lower();
               h = low * a.lower();
               }
         }
      else
         { // low and high are < 0
         if( a.lower() >= _IT(0) )
            { // a.low >=0, a.high >= 0
            l = low * a.upper();
            h = high * a.lower();
            }
         else
            if( a.upper() >= _IT(0) )
               {  //  a.low < 0, a.high >= 0
               l = low * a.upper();
               h = low * a.lower();
               }
            else
               { // a.low and a.high < 0
               l = high * a.upper();
               h = low * a.lower();
               }
        }

   low = l;
   high = h;

   return *this;
   }


// Specialization for float_precision and *= operator
//
template<> inline interval<float_precision>& interval<float_precision>::operator*=( const interval<float_precision>& a )
   {
   float_precision l, h, t;

   l.precision( low.precision() );
   h.precision( low.precision() );
   t.precision( low.precision() );

   l.mode( ROUND_DOWN );
   h.mode( ROUND_UP );

   if( low.sign() > 0 ) //
      { // both low and high >= 0
      if( a.lower().sign() > 0 )
         { // a.low >=0, a.high >= 0
         l = low;  l *= a.lower();
         h = high; h *= a.upper();
         }
      else
         if( a.upper().sign() > 0 )
            {  //  a.low < 0, a.high >= 0
            l = high;  l *= a.lower();
            h = high; h *= a.upper();
            }
         else
            { // a.low and a.high < 0
            l = high; l *= a.lower();
            h = low;  h *= a.upper();
            }
      }
   else
      if( high.sign() > 0 )
         {  // low < 0, high >= 0
         if( a.lower().sign() > 0 )
            { // a.low >=0, a.high >= 0
            l = low;  l *= a.upper();
            h = high; h *= a.upper();
            }
         else
            if( a.upper().sign() > 0 )
               {  //  a.low < 0, a.high >= 0
               t.mode( ROUND_DOWN );
               l = low;  l *= a.upper(); 
			   if( l > ( t = high, t *= a.lower() ) ) 
				   l = t;
               t.mode( ROUND_UP );
               h = high; h *= a.upper(); 
			   if( h < ( t = low, t *= a.lower() ) ) 
				   h = t;
               }
            else
               { // a.low and a.high < 0
               l = high; l *= a.lower();
               h = low;  h *= a.lower();
               }
         }
      else
         { // low and high are < 0
         if( a.lower().sign() > 0 )
            { // a.low >=0, a.high >= 0
            l = low;  l *= a.upper();
            h = high; h *= a.lower();
            }
         else
            if( a.upper().sign() > 0 )
               {  //  a.low < 0, a.high >= 0
               l = low; l *= a.upper();
               h = low; h *= a.lower();
               }
            else
               { // a.low and a.high < 0
               l = high; l *= a.upper();
               h = low; h *= a.lower();
               }
        }

   low = l;
   high = h;
   return *this;
   }

// Specilization for float and *=
// That can work with both managed an unmanaged application.
// HARDWARE_SUPPORT indicate if we use the build in ROUND_DOWN or ROUND_UP Mode when performing the operation
// For now Hardware support we do software emulation to do the interval arithmetic in software
//
template<> inline interval<float>& interval<float>::operator*=( const interval<float>& a )
	{
	float l, h, t;

	if (low >= 0) //
		{ // both low and high >= 0
		if (a.lower() >= 0)
			{ // a.low >=0, a.high >= 0
			l = tofloat(mul_down(low, a.lower()), ROUND_DOWN );
			h = tofloat(mul_up(high, a.upper()), ROUND_UP );
			}
		else
			if (a.upper() >= 0)
				{  //  a.low < 0, a.high >= 0
				l = tofloat(mul_down(high, a.lower()), ROUND_DOWN );
				h = tofloat(mul_up(high, a.upper()), ROUND_UP );
			}
			else
				{ // a.low and a.high < 0
				l = tofloat(mul_down(high, a.lower()), ROUND_DOWN );
				h = tofloat(mul_up(low, a.upper()), ROUND_UP );
				}
		}
	else
		if (high >= 0)
			{  // low < 0, high >= 0
			if (a.lower() >= 0)
				{ // a.low >=0, a.high >= 0
				l = tofloat(mul_down(low, a.upper()), ROUND_DOWN );
				h = tofloat(mul_up(high, a.upper()), ROUND_UP );
				}
			else
				if (a.upper() >= 0)
					{  //  a.low < 0, a.high >= 0
					l = tofloat(mul_down(low, a.upper()), ROUND_DOWN );
					t = tofloat(mul_down(high, a.lower()), ROUND_DOWN);
					if (l > t) 
						l = t;
					h = tofloat(mul_up(high, a.upper()), ROUND_UP ); 
					t = tofloat(mul_up(low, a.lower()), ROUND_UP);
					if (h < t )
						h = t;
				}
				else
					{ // a.low and a.high < 0
					l = tofloat(mul_down(high, a.lower()), ROUND_DOWN );
					h = tofloat(mul_up(low, a.lower()), ROUND_UP );
					}
			}
		else
			{ // low and high are < 0
			if (a.lower() >= 0)
				{ // a.low >=0, a.high >= 0
				l = tofloat(mul_down(low, a.upper()), ROUND_DOWN );
				h = tofloat(mul_up(high, a.lower()), ROUND_DOWN );
				}
			else
				if (a.upper() >= 0)
					{  //  a.low < 0, a.high >= 0
					l = tofloat(mul_down(low, a.upper()), ROUND_DOWN);
					h = tofloat(mul_up(low, a.lower()), ROUND_UP );
					}
				else
					{ // a.low and a.high < 0
					l = tofloat(mul_down(high, a.upper()), ROUND_DOWN );
					h = tofloat(mul_up(low, a.lower()), ROUND_UP );
					}
			}

	low = l;
	high = h;
	return *this;
	}

// Specilization for double and *=
// That can work with both managed an unmanaged application.
// HARDWARE_SUPPORT indicate if we use the build in ROUND_DOWN or ROUND_UP Mode when performing the operation
// For now Hardware support we do software emulation to do the interval arithmetic in software
//
template<> inline interval<double>& interval<double>::operator*=( const interval<double>& a )
	{
	double l, h, t;

	if (low >= 0) //
	{ // both low and high >= 0
		if (a.lower() >= 0)
		{ // a.low >=0, a.high >= 0
			l = mul_down(low, a.lower());
			h = mul_up(high, a.upper());
		}
		else
			if (a.upper() >= 0)
			{  //  a.low < 0, a.high >= 0
				l = mul_down(high, a.lower());
				h = mul_up(high, a.upper());
			}
			else
			{ // a.low and a.high < 0
				l = mul_down(high, a.lower());
				h = mul_up(low, a.upper());
			}
	}
	else
		if (high >= 0)
		{  // low < 0, high >= 0
			if (a.lower() >= 0)
			{ // a.low >=0, a.high >= 0
				l = mul_down(low, a.upper());
				h = mul_up(high, a.upper());
			}
			else
				if (a.upper() >= 0)
				{  //  a.low < 0, a.high >= 0
					l = mul_down(low, a.upper());
					t = mul_down(high, a.lower());
					if (l > t)
						l = t;
					h = mul_up(high, a.upper());
					t = mul_up(low, a.lower());
					if (h < t)
						h = t;
				}
				else
				{ // a.low and a.high < 0
					l = mul_down(high, a.lower());
					h = mul_up(low, a.lower());
				}
		}
		else
		{ // low and high are < 0
			if (a.lower() >= 0)
			{ // a.low >=0, a.high >= 0
				l = mul_down(low, a.upper());
				h = mul_up(high, a.lower());
			}
			else
				if (a.upper() >= 0)
				{  //  a.low < 0, a.high >= 0
					l = mul_down(low, a.upper());
					h = mul_up(low, a.lower());
				}
				else
				{ // a.low and a.high < 0
					l = mul_down(high, a.upper());
					h = mul_up(low, a.lower());
				}
		}

	low = l;
	high = h;
	return *this;
	}

// Works for all other classes
// Please note that this is for all interger classes. interval<int>, interval<long>
// were there is no loss of precision
// Actually there is specialization for both <int> and <int_precision> further down.
//
template<class _IT> inline interval<_IT>& interval<_IT>::operator/=( const interval<_IT>& b )
   {
   interval<_IT> a, c;

   c.low = (_IT)1 / b.upper();
   c.high = (_IT)1 / b.lower();

   a = interval( low, high );
   c *= a;

   low = c.lower();
   high = c.upper();

   return *this;
   }

// Specilization for float_precision and /=
//
template<> inline interval<float_precision>& interval<float_precision>::operator/=( const interval<float_precision>& b )
   {
   float_precision l, h;
   interval<float_precision> c(b);

   l.precision(b.upper().precision());
   l = b.upper();
   l.mode( ROUND_DOWN );
   l = _float_precision_inverse( l );

   h.precision(b.lower().precision());
   h = b.lower();
   h.mode( ROUND_UP );
   h = _float_precision_inverse( h );

   c = interval<float_precision>( l , h );
   *this *= c;

   return *this;
   }

// Specilization for float and /=
// That can work with both managed an unmanaged application.
// HARDWARE_SUPPORT indicate if we use the build in ROUND_DOWN or ROUND_UP Mode when performing the operation
// For now Hardware support we software emulation to do the interval arithmetic in software
//
template<> inline interval<float>& interval<float>::operator/=( const interval<float>& a )
	{
	interval<float> b, c;

	c.low = tofloat(div_down(1, a.upper()), ROUND_DOWN );
	c.high = tofloat(div_up(1, a.lower()), ROUND_UP );
	b = interval(low, high);
	c *= b;

	low = c.lower();
	high = c.upper();
	return *this;
	}

// Specilization for double and /=
// That can work with both managed an unmanaged application.
// HARDWARE_SUPPORT indicate if we use the build in ROUND_DOWN or ROUND_UP Mode when performing the operation
// For now Hardware support we do software emulation to do the interval arithmetic in software
//
template<> inline interval<double>& interval<double>::operator/=( const interval<double>& a )
	{
	interval<double> b, c;

	c.low = div_down(1, a.upper());
	c.high = div_up(1, a.lower());
	b = interval(low, high);
	c *= b;

	low = c.lower();
	high = c.upper();
	return *this;
	}

// Specilization for int_precision and /=
//
template<> inline interval<int_precision>& interval<int_precision>::operator/=( const interval<int_precision>& b )
   {
   float_precision l = b.upper(), h = b.lower();
   interval<float_precision> c(b), a(*this);

   l.mode( ROUND_DOWN );
   l = _float_precision_inverse( l );

   h.mode( ROUND_UP );
   h = _float_precision_inverse( h );

   c = interval<float_precision>( l , h );
   a *= c;

   low = (int_precision)floor(a.lower());
   high = (int_precision)ceil(a.upper()); 

   return *this;
   }


// Specialization for int and /=
//
template<> inline interval<int>& interval<int>::operator/=( const interval<int>& b )
   {
   double tlow, thigh;
   interval<int> a;
   interval<double> c;

   tlow = 1 / (double)b.upper();
   thigh = 1 / (double)b.lower();

   a = interval( low, high );
   c = interval<double>( tlow, thigh );
   c *= a;

   low = (int)floor( c.lower() );
   high = (int)ceil( c.upper() );

   return *this;
   }

// Works on all classes.
// Return the intersection
//
template<class _IT> inline interval<_IT>& interval<_IT>::operator&=(const interval<_IT>& a)
	{
	if (a.lower() > low )
		low = a.lower();
	if (a.upper() < high)
		high = a.upper();
	if (low > high)  // Empty set
		{
		low = 0; high = 0;
		}

	return *this;
	}

// Works on all classes.
// Return the union
//
template<class _IT> inline interval<_IT>& interval<_IT>::operator|=(const interval<_IT>& a)
	{
	if (low > a.upper() || high < a.lower())
		{
		if (a.upper() - a.lower() > high - low)
			{ // return the largest set
			low = a.lower();
			high = a.upper();
			}
		}
	else
		{ // non empty intersection
		if (a.lower() < low)
			low = a.lower();
		if (a.upper() > high)
			high = a.upper();
		}
	}


// Works on all classes.
// Return the set minus
//
template<class _IT> inline interval<_IT>& interval<_IT>::operator^=(const interval<_IT>& a)
	{
	if ( a.lower() < high && a.upper() > low ) // intersection is not empty
		{
		if (a.upper() <= low)
			low = a.upper();
		else
			if (a.high() >= high)
				high = a.lower();
		}

	return *this;
	}

//////////////////////////////////////////////////////////////////////////////////////
///
/// END Essential Operators
///
//////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//    Binary and Unary Operators +,-,*,/
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// Binary + operator
// Specialization for float_precision
//
inline interval<float_precision> operator+(const interval<float_precision>& a, const float_precision& b)
	{
	interval<float_precision> c(a);

	c += interval<float_precision>(b);
	return c;
	}

// Binary + operator
// Works for all classes
//
inline interval<float_precision> operator+(float_precision& a, const interval<float_precision>& b )
{
	interval<float_precision> c(b);

	c += interval<float_precision>(a);
	return c;
}

// Binary + operator
// Works for all classes
//
template<class _IT,class _X> inline interval<_IT> operator+(const interval<_IT>& a, const _X& b)
	{
	interval<_IT> c(a);

	c += interval<_IT>(_IT(b));
	return c;
	}

// Binary + operator
// Works for all classes
//
template<class _IT,class _X> inline interval<_IT> operator+( const _X& a, const interval<_IT>& b)
	{
	interval<_IT> c(b);

	c += interval<_IT>(_IT(a));
	return c;
	}


// Binary + operator
// Works for all classes
//
template<class _IT> inline interval<_IT> operator+(const interval<_IT>& a, const interval<_IT>& b)
	{
	interval<_IT> c(a);

	c += b;
	return c;
	}


// Unary + operator
// Works for all classes
//
template<class _IT> inline interval<_IT> operator+( const interval<_IT>& a )
   {
   return a;
   }


// Binary - operator
// Specialization for float_precision
//
inline interval<float_precision> operator-(const interval<float_precision>& a, const float_precision& b)
	{
	interval<float_precision> c(a);

	c -= interval<float_precision>(b);
	return c;
	}

// Binary - operator
// Works for all classes
//
inline interval<float_precision> operator-(float_precision& a, const interval<float_precision>& b)
	{
	interval<float_precision> c(a);

	c -= b;
	return c;
	}

// Binary - operator
// Works for all classes
//
template<class _IT, class _X> inline interval<_IT> operator-(const interval<_IT>& a, const _X& b)
	{
	interval<_IT> c(a);

	c -= interval<_IT>(_IT(b));
	return c;
	}

// Binary - operator
// Works for all classes
//
template<class _IT, class _X> inline interval<_IT> operator-(const _X& a, const interval<_IT>& b)
	{
	interval<_IT> c(a);

	c -= b;
	return c;
	}

// Binary - operator
// Works for all classes
//
template<class _IT> inline interval<_IT> operator-( const interval<_IT>& a, const interval<_IT>& b )
   {
   interval<_IT> c(a);

   c -= b;
   return c;
   }


// Unary - operator
// Works for all classes
//
template<class _IT> inline interval<_IT> operator-( const interval<_IT>& a )
   {
   interval<_IT> c(0);

   c -= a;
   return c;
   }

// Binary * operator
// Specialization for float_precision
//
inline interval<float_precision> operator*(const interval<float_precision>& a, const float_precision& b)
	{
	interval<float_precision> c(a);

	c *= interval<float_precision>(b);
	return c;
	}

// Binary * operator
// Works for all classes
//
inline interval<float_precision> operator*(float_precision& a, const interval<float_precision>& b)
	{
	interval<float_precision> c(b);

	c *= interval<float_precision>(a);
	return c;
	}

// Binary * operator
// Works for all classes
//
template<class _IT, class _X> inline interval<_IT> operator*(const interval<_IT>& a, const _X& b)
	{
	interval<_IT> c(a);

	c *= interval<_IT>(_IT(b));
	return c;
	}

// Binary * operator
// Works for all classes
//
template<class _IT, class _X> inline interval<_IT> operator*(const _X& a, const interval<_IT>& b)
	{
	interval<_IT> c(b);

	c *= interval<_IT>(_IT(a));
	return c;
	}



// Binary * operator
// Works for all classes
//
template<class _IT> inline interval<_IT> operator*( const interval<_IT>& a, const interval<_IT>& b )
   {
   interval<_IT> c(a);

   c *= b;
   return c;
   }

// Binary / operator
// Specialization for float_precision
//
inline interval<float_precision> operator/(const interval<float_precision>& a, const float_precision& b)
	{
	interval<float_precision> c(a);

	c /= interval<float_precision>(b);
	return c;
	}

// Binary / operator
// Works for all classes
//
inline interval<float_precision> operator/(float_precision& a, const interval<float_precision>& b)
	{
	interval<float_precision> c(a);

	c /= interval<float_precision>(b);
	return c;
	}

// Binary / operator
// Works for all classes
//
template<class _IT, class _X> inline interval<_IT> operator/(const interval<_IT>& a, const _X& b)
	{
	interval<_IT> c(a);

	c /= interval<_IT>(_IT(b));
	return c;
	}

// Binary / operator
// Works for all classes
//
template<class _IT, class _X> inline interval<_IT> operator/(const _X& a, const interval<_IT>& b)
	{
	interval<_IT> c(a);

	c /= b;
	return c;
	}


// Binary / operator
// Works for all classes
//
template<class _IT> inline interval<_IT> operator/( const interval<_IT>& a, const interval<_IT>& b )
   {
   interval<_IT> c(a);

   if ( c == b && b.is_class() != ZERO )
	  c = interval<_IT>(1,1);
   else
      c /= b;

   return c;
   }

// Binary & operator
// Return intersection
// Works for all classes
//
template<class _IT> inline interval<_IT> operator&( const interval<_IT>& a, const interval<_IT>& b )
	{
	interval<_IT> c(a);

	c &= b;
	return c;
	}

// Binary | operator.
// Return union
// Works for all classes
//
template<class _IT> inline interval<_IT> operator|(const interval<_IT>& a, const interval<_IT>& b)
	{
	interval<_IT> c(a);

	c |= b;
	return c;
	}

// Binary ^ operator
// Return set minus
// Works for all classes
//
template<class _IT> inline interval<_IT> operator^(const interval<_IT>& a, const interval<_IT>& b)
	{
	interval<_IT> c(a);

	c ^= b;
	return c;
	}

//////////////////////////////////////////////////////////////////////////////////////
///
/// END Binary and Unary Operators
///
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
///
/// Boolean Interval for == and !=
///
//////////////////////////////////////////////////////////////////////////////////////

template<class _IT> inline bool operator==(const interval<_IT>& a, const interval<_IT>& b)
	{
	return a.lower() == b.lower() && a.upper() == b.upper();
	}

template<class _IT> inline bool operator!=(const interval<_IT>& a, const interval<_IT>& b)
	{
	return a.lower() != b.lower() || a.upper() != b.upper();
	}

//////////////////////////////////////////////////////////////////////////////////////
///
/// END Boolean operators
///
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
///
/// Interval abs(), interior, precedes, intersection, unionsection
///
//////////////////////////////////////////////////////////////////////////////////////

template<class _IT> inline interval<_IT> abs( const interval<_IT>& a )
	{
	if (a.lower() >= _IT(0) )
		return a;
	else
		if (a.upper() <= _IT(0) )
			return -a;

	return interval<_IT>(_IT(0), ( a.upper() > -a.lower() ? a.upper() : -a.lower() ) );
	}

// True if interval a is fully within interval b (interior of b) otherwise false
template<class _IT> inline bool interior(const interval<_IT>& a, const interval<_IT>& b )
	{
	if (b.lower() < a.lower() && a.upper() < b.upper())
		return true;
	else
		return false;
	}

// true if interval a precedes interval b otherwise false
template<class _IT> inline bool precedes(const interval<_IT>& a, const interval<_IT>& b)
	{
	if (a.upper() < b.lower() )
		return true;
	else
		return false;
	}

// Return Intersection of a & b
template<class _IT> inline interval<_IT> intersection(const interval<_IT>& a, const interval<_IT>& b)
	{
	_IT l, u;
	l = std::max(a.lower(), b.lower());
	u = std::min(a.upper(), b.upper());
	return interval<_IT>(l,u);
	}

// Return Union of a & b
template<class _IT> inline interval<_IT> unionsection(const interval<_IT>& a, const interval<_IT>& b)
	{
	_IT l, u;
	l = std::min(a.lower(), b.lower());
	u = std::max(a.upper(), b.upper());
	return interval<_IT>(l, u);
	}
//////////////////////////////////////////////////////////////////////////////////////
///
/// END interval functions
///
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
///
/// Interval sqrt(), log10(), log(), exp() and pow()
///
//////////////////////////////////////////////////////////////////////////////////////

// Support function for correctly converting and float_precision number back to a float
// By default IEE754 round to nearest and that will create incorrect result for interval arithmetic
//
inline float tofloat(const interval<float_precision>& fp, enum round_mode rm )
{
	float f;
	if (rm == ROUND_DOWN) // Conversion to float introduce an error since it is always round to nearest and it should be round up or round_down
		f = fp.lower();
	else
		f = fp.upper();

	if (f != 0)
		{
		float_precision fp1(f, 9);
		if (rm == ROUND_DOWN && fp1 > fp.lower())
			f -= f * 0.5f * FLT_EPSILON;
		if (rm == ROUND_UP && fp1 < fp.upper())
			f += f * 0.5f * FLT_EPSILON;
		}
	return f;
}
/*
inline float tofloat(const interval<double>& fp, enum round_mode rm)
	{
	float f;
	double d;
	if (rm == ROUND_DOWN) // Conversion to float introduce an error since it is always round to nearest and it should be round up or round_down
		f = d = fp.lower();
	else
		f = d = fp.upper();

	if (f != 0)
		{
		if (f != d && rm == ROUND_UP)
			f= nextafterf(f,1.0);
		}
	return f;
	}*/

// Support function for correctly converting and float_precision number back to a float
// By default IEE754 round to nearest and that will create incorrect result for interval arithmetic
//
inline double todouble(const interval<float_precision>& fp, enum round_mode rm)
{
	double  d;
	if (rm == ROUND_DOWN) // Conversion to float introduce an error since it is always round to nearest and it should be round up or round_down
		d = fp.lower();
	else
		d = fp.upper();

	if (d != 0)
		{
		float_precision fp1(d, 17 );
		if (rm == ROUND_DOWN && fp1 > fp.lower())
			d -= d * 0.5 * DBL_EPSILON;
		if (rm == ROUND_UP && fp1 < fp.upper())
			d += d * 0.5 * DBL_EPSILON;
		}
	return d;
}


static interval<double> _intervalexpdouble(const double x)
	{
	int i, k, sign = 1;
	double xn = x;
	const interval<double> c1(1), c2(2);
	interval<double> isum, ix, ixp, ifac, idelta;
	if (xn<0) { xn = -xn; sign = -1; }
	// Argument reduction
	for (k = 0; false && xn > 0.5&& k <= 16; ++k) xn *= 0.5;
	ix = interval<double>(xn); ixp = c1; isum = c1; ifac = interval<double>(1);
	for (i = 1;; ++i)
		{
		ifac *= interval<double>(i);
		ixp *= ix;
		idelta = ixp /interval<double>(ifac);
		if (isum.center() + idelta.center() == isum.center()) 
			break;
		isum += idelta;
		}
	// Reverse reduction
	// Brent enhancement avoid loss of significant digits when x is small.
	if (k>0)
		{
		isum -= c1;
		for (; k > 0; k--)
			isum = (c2 + isum)*isum;
		isum += c1;
		}

	//for (; k>0; --k) 
	//	isum *= isum;
	if (sign<0) 
		isum = c1 / isum;
	return isum;
	}


static interval<double> _intervallogdouble(const double x)
	{
	int i, k, expo = 0;
	interval<double> zn, zsq, sum, delta;
	double xd;
	if (x<0) { throw interval<double>::domain_error(); }
	if (x == 0) { return interval<double>(-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()); }
	if (x == 1) { return interval<double>(0); }
	// Remove the exponent
	expo = (int)trunc(log10(x)); 
	xd = x / pow(10,expo);
	zn = interval<double>(xd);
	// Taylor series of log(x)
	// log(x)=2( z + z^3/3 + z^5/5 ...)
	// where z=(x-1)/(x+1) 
	// Argument reduction
	// In order to get a fast Taylor series result we need to get the fraction closer to 1
	// The fraction part is [1.xxx-9.999] (base 10) 
	// Repeat a series of square root until z < 1.1 
	for (k = 0; false && zn.center()>1.1; ++k) 
		zn = sqrt(zn);
	// Initialize the iteration
	zn = (zn-interval<double>(1))/ (zn+interval<double>(1));
	zsq = zn * zn; 
	sum = zn;
	// Iterate using taylor series log(x) == 2( z + z^3/3 + z^5/5 ... ) 
	for (i = 3;; i += 2)
		{
		zn *= zsq;
		delta = zn/interval<double>(i);
		if (sum.center() + delta.center() == sum.center())
			break;
		sum += delta;
		}
	sum *= interval<double>(pow(2,(k + 1)));
	// Handle expo adjustment
	if (expo != 0)
		{
		sum += interval<double>(expo)*interval<double>(log(10));
		}
	return sum;
	}

// Specilization for sqrt(float_precision)
//
inline interval<float_precision> sqrt( const interval<float_precision>& x )
   {
   float_precision l, u;

   l.assign( x.lower() );  // Assign value, precision and mode
   l.mode( ROUND_DOWN );
   l = sqrt( l );

   u.assign( x.upper() );  // Assign value, precision and mode
   u.mode( ROUND_UP );
   u = sqrt( u );

   return interval<float_precision>( l, u );
   }

// sqrt for float using managed code.
//
inline interval<float> sqrt( const interval<float>& x )
   {
   float lower, upper;
   lower = tofloat(sqrt_down(x.lower()), ROUND_DOWN );
   upper = tofloat(sqrt_up(x.upper()), ROUND_UP );
   return interval<float>( lower, upper );
   }

// sqrt for double using managed code.
//
inline interval<double> sqrt( const interval<double>& x )
   {
   double lower, upper;
   lower = sqrt_down(x.lower());
   upper = sqrt_up(x.upper());
   return interval<double>( lower, upper );
   }



// Specilization for log float_precision
//
inline interval<float_precision> log( const interval<float_precision>& x )
   {
   float_precision l, u;

   l.assign( x.lower() );  // Assign value, precision and mode
   l.mode( ROUND_DOWN );
   l = log( l );

   u.assign( x.upper() );  // Assign value, precision and mode
   u.mode( ROUND_UP );
   u = log( u );

   return interval<float_precision>( l, u );
   }

// log for float using managed code.
//
inline interval<float> log( const interval<float>& x )
	{
	float lower, upper;
	lower = tofloat(_intervallogdouble(x.lower()).lower(), ROUND_DOWN);
	upper = tofloat(_intervallogdouble(x.upper()).upper(), ROUND_UP );
	return interval<float>( lower, upper );
	}

// log for double using managed code.
//
inline interval<double> log( const interval<double>& x )
	{
	double lower, upper;
	lower = _intervallogdouble(x.lower()).lower();
	upper = _intervallogdouble(x.upper()).upper();
	return interval<double>( lower, upper );
	}



// Specilization for log float_precision
//
inline interval<float_precision> log10( const interval<float_precision>& x )
   {
   float_precision l, u;

   l.assign( x.lower() );  // Assign value, precision and mode
   l.mode( ROUND_DOWN );
   l = log10( l );

   u.assign( x.upper() );  // Assign value, precision and mode
   u.mode( ROUND_UP );
   u = log10( u );

   return interval<float_precision>( l, u );
   }



// log10 for float using managed code.
//
inline interval<float> log10( const interval<float>& x )
   {
   float lower, upper;
   interval<float> tmp;
   lower = tofloat(_intervallogdouble(x.lower()).lower(), ROUND_DOWN );
   upper = tofloat(_intervallogdouble(x.upper()).upper(), ROUND_UP );
   tmp = LN10;
   tmp = interval<float>(lower, upper)/tmp;
   lower = tmp.lower();
   upper = tmp.upper();
   return interval<float>( lower, upper );
   }

// log10 for double using managed code.
//
inline interval<double> log10( const interval<double>& x )
	{
	double lower, upper;
	interval<double> tmp;
	lower = _intervallogdouble(x.lower()).lower();
	upper = _intervallogdouble(x.upper()).upper();
	tmp = LN10;
	tmp = interval<double>(lower, upper)/tmp;
	lower = tmp.lower();
	upper = tmp.upper();
	return interval<double>( lower, upper );
	}


// Specilization for exp float_precision
//
inline interval<float_precision> exp( const interval<float_precision>& x )
   {
   float_precision l, u;

   l.assign( x.lower() );  // Assign value, precision and mode
   l.mode( ROUND_DOWN );
   l = exp( l );

   u.assign( x.upper() );  // Assign value, precision and mode
   u.mode( ROUND_UP );
   u = exp( u );

   return interval<float_precision>( l, u );
   }

#ifdef TEST_ONLY_HVE
// exp for float using managed code.
// ONly for test purposed.
//
inline interval<float> exp2( const interval<float>& x )
   {
	interval<float_precision> fx(x);
   float lower, upper;

   fx=exp(fx);
   lower = tofloat(fx, ROUND_DOWN);
   upper = tofloat(fx, ROUND_UP);
   return interval<float>( lower, upper );
   }

// exp for double using managed code.
//
inline interval<double> exp2( const interval<double>& x )
   {
	interval<float_precision> fx(x);
   double lower, upper;

   fx=exp(fx);
   lower = todouble(fx, ROUND_DOWN);
   upper = todouble(fx, ROUND_UP);
   return interval<double>( lower, upper );
   }
#endif



// MSC exp() does not allow rounding control
// So we have to do it manually
// Use a taylor series until their is no more change in the result
// exp(x) == 1 + x + x^2/2!+x^3/3!+....
// Equivalent with the same standard C function call
// use argument reduction via exp(x)=(exp(x/2^k)2^k
// And use Brent enhancement using the double formula:
// expm(x)=exp(x)-1 && expm(2x)=expm(x)(2+expm(x)) on the backend to preseve
// loss of significance digits
//
inline interval<double> exp(const interval<double>& x)
	{
	double lower, upper;
	lower = _intervalexpdouble(x.lower()).lower();
	upper = _intervalexpdouble(x.upper()).upper();
	return interval<double>(lower, upper);
	}

// MSC exp() does not allow rounding control for the exp()
// Since we dont normally do it using float arithmetic (as for sqrt(), log() and log10()) we simply just call the interval<double> version
// of exp() and convert back to float preserving as much accuracy as possible
//
inline interval<float> exp(const interval<float>& x)
	{
	float lower, upper;
	lower = tofloat(_intervalexpdouble(x.lower()).lower(), ROUND_DOWN );
	upper = tofloat(_intervalexpdouble(x.upper()).upper(), ROUND_UP );
	return interval<float>(lower, upper);
	}

// Specilization for pow float_precision
//
inline interval<float_precision> pow( const interval<float_precision>& x, const float_precision& y )
   {
   interval<float_precision> c(x);

   c = log( x );
   c *= interval<float_precision>( y );
   c = exp( c );

   return c;
   }



// MSC pow() does not allow rounding control
// So we have to do it manually
// x^y == exp( y * ln( x ) ) );
// To avoid loss of precision we actually perform the operation using double and then
// convert the result back to float. This is consistent with the pow() that only takes double as an argument.
//
inline interval<float> pow(const interval<float>& x, const float y)
	{
	interval<double> c(x);
	float upper, lower;

	c = log(c);
	c *= interval<double>(y);
	c = exp(c);
	lower = tofloat(c.lower(), ROUND_DOWN );
	upper = tofloat( c.upper(), ROUND_UP );
	/*if (lower > c.lower() )
		lower -= lower * 0.5f * FLT_EPSILON;
	if (upper < c.upper() )
		upper += upper * 0.5f * FLT_EPSILON;
*/
	return interval<float>(lower, upper);
	}

// MSC pow() does not alllow rounding control
// So we have to do it manually
// x^y == exp( y * ln( x ) ) );
//
inline interval<double> pow(const interval<double>& x, const double y)
	{
	interval<double> c(x);

	c = log(c);
	c *= interval<double>(y);
	c = exp(c);

	return c;
	}

//////////////////////////////////////////////////////////////////////////////////////
///
/// END Interval sqrt(), log10(), log(), exp(), pow()
///
//////////////////////////////////////////////////////////////////////////////////////

// Load manifest constant PI for float
//
inline interval<float> int_pifloat()
	{
	interval<float> pif;
	pif.lower(tofloat(PI.lower(), ROUND_DOWN));
	pif.upper(tofloat(PI.upper(), ROUND_UP));
	return pif;;
	}

// Specilization for constant PI float_precision
//
inline interval<float_precision> int_pi(const unsigned int p = float_precision_ctrl.precision() )
	{
	float_precision fx(0,p+1), l(0,p), u(0,p);

	fx = _float_table(_PI, p+1);
	l.mode(ROUND_DOWN);
	l = fx;
	u.mode(ROUND_UP);
	u = fx;

	return interval<float_precision>(l, u);
	}


// Load manifest constant LN2 for float
//
inline interval<float> int_ln2float()
	{
	interval<float> ln2f;

	ln2f.lower(tofloat(LN2.lower(), ROUND_DOWN));
	ln2f.upper(tofloat(LN2.upper(), ROUND_UP));
	return ln2f;;
	}

// Specilization for constant LN2 float_precision
//
inline interval<float_precision> int_ln2(const unsigned int p = float_precision_ctrl.precision())
	{
	float_precision fx(0, p + 1), l(0, p), u(0, p);

	fx = _float_table(_LN2, p + 1);
	l.mode(ROUND_DOWN);
	l = fx;
	u.mode(ROUND_UP);
	u = fx;

	return interval<float_precision>(l, u);
	}

// Load manifest constant LN10 for float
//
inline interval<float> int_ln10float()
	{
	interval<float> ln10f;

	ln10f.lower(tofloat(LN10.lower(), ROUND_DOWN));
	ln10f.upper(tofloat(LN10.upper(), ROUND_UP));
	return ln10f;
	}

// Specilization for constant LN10 float_precision
//
inline interval<float_precision> int_ln10(const unsigned int p = float_precision_ctrl.precision() )
	{
	float_precision fx(0, p + 2), l(0, p), u(0, p);

	fx = _float_table(_LN10, p + 2);
	l.mode(ROUND_DOWN);
	l = fx;
	u.mode(ROUND_UP);
	u = fx;

	return interval<float_precision>(l, u);
	}

// Load manifest constant E for double
//
inline interval<double> int_Edouble()
	{
	const interval<double> e(2.7182818284590451, 2.7182818284590455);
	return e;
	}

// Load manifest constant LN2 for float
//
inline interval<float> int_Efloat()
	{
	interval<double> ed;
	interval<float> ef;

	ed = int_Edouble();
	ef.lower(tofloat(ed.lower(), ROUND_DOWN));
	ef.upper(tofloat(ed.upper(), ROUND_UP));
	return ef;
	}

// Specilization for constant E float_precision
//
inline interval<float_precision> int_E(const unsigned int p = float_precision_ctrl.precision())
	{
	float_precision fx(0, p + 2), l(0, p), u(0, p);

	fx = _float_table(_EXP1, p + 2);
	l.mode(ROUND_DOWN);
	l = fx;
	u.mode(ROUND_UP);
	u = fx;

	return interval<float_precision>(l, u);
	}

//////////////////////////////////////////////////////////////////////////////////////
///
/// END Interval constants
///
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
///
/// Interval sin(), cos(), tan(), asin(), acos(), atan()
///
//////////////////////////////////////////////////////////////////////////////////////

static interval<double> _intervalsindouble(const double x)
	{
	int i, k, sign = 1;
	double xn = x;
	const interval<double> c3(3), c4(4);
	const double PID = 3.1415926535897931;
	interval<double> zn, zsq, isum, idelta;
	if (xn<0) { xn = -xn; sign = -1; }
	// Argument reduction
	//  x has been reduce to between 0..2*PI prior to call
	//if (xn >= 2 * PI)
	//	xn = fmod( xn, 2 * PI);
	//   2) Then reduced further to between 0..PI using sin(x+PI)=-Sin(x) 
	if (xn >= PID)
		{
		xn -= PID; sign *= -1;
		}
	//   3) Finally reduced it to below 0.5/3^reduction factor, using the trisection identity 
	//   The argument reduction is used to reduced the number of taylor iteration and to minimize round off erros and calculation time  
	zn = interval<double>(xn);
	for (k = 0; false && zn.center()>=1; ++k)
		zn /= c3;
	//   4) Then do the taylor.  
	// Taylor series of sin(x)
	// Sin(x) = x - x^3/3! + x^5/5!-...
	// Initialize the iteration
	zsq = zn * zn;  idelta=isum = zn; 
	// Iterate using taylor series sin(x)=x-x^3/3!+x^5/5!...
	for (i = 3;; i += 2)
		{
		idelta = idelta*zsq / interval<double>(i*(i - 1));
		idelta = -idelta; 
		if (isum.center() + idelta.center() == isum.center()) break;
		isum += idelta;
		}
	// Reverse reduction
	for (; k>0; --k) // sin(3x)=sin(x)(3-4*sin(x)^2)
		isum = isum * ( c3 - c4 * isum * isum);
	if (sign<0) isum = -isum;
	return isum;
	}

static interval<double> _intervalcosdouble(const double x)
	{
	int i, k, sign = 1;
	double xn = x;
	const interval<double> c1(1), c3(3), c4(4);
	const double PID = 3.1415926535897931;
	interval<double> zn, zsq, isum, idelta;
	if (xn<0) { xn = -xn; sign = -1; }
	// Argument reduction
	//  x has been reduce to between 0..2*PI prior to call
	//if (xn >= 2 * PI)
	//	xn = fmod(xn, 2 * PI);
	//   2) Then reduced further to between 0..PI using cos(x)=cos(2PI-x) 
	if (xn >= PID)
		xn = 2* PID - xn; 
	//   3) Finally reduced it to below 0.5/3^reduction factor, using the trisection identity 
	//   The argument reduction is used to reduced the number of taylor iteration and to minimize round off erros and calculation time  
	zn = interval<double>(xn);
	for (k = 0; false && zn.center()>=1; ++k)
		zn /= c3;
	//   4) Then do the taylor.  
	// Taylor series of sin(x)
	// Cos(x) = 1 - x^2/2! + x^4/4!-...
	// Initialize the iteration
	zsq = zn * zn;  idelta=isum = c1;
	// Iterate using taylor series  Cos(x) = 1 - x^2/2! + x^4/4!-...
	for (i = 2;; i += 2)
		{
		idelta = idelta * zsq / interval<double>(i*(i-1));
		idelta = -idelta; 
		if (isum.center() + idelta.center() == isum.center()) break;
		isum += idelta;
		}
	// Reverse reduction
	for (; k>0; --k) // cos(3x)=cos(x)(-3+4*cos(x)^2)
		isum = isum * (-c3 - c4 * isum * isum);
	return isum;
	}

//  Use the identity tan(x)=Sin(x)/Sqrt(1-Sin(x)^2)
static interval<double> _intervaltandouble(const double x)
	{
	const interval<double> c1(1);
	const double PID = 3.1415926535897931;
	interval<double> z, zsq;
	//  x has been reduce to between 0..2*PI prior to call
	if (x == PID || x == PID * 3 )
		{
		throw interval<double>::domain_error();
		}
	z = sin(interval<double>(x));
	zsq = sqrt(c1 - z*z);
	if (abs(x) >= PID*0.5&& abs(x) <= PID*1.5)
		zsq = -zsq;
	z /= zsq;
	return z;
	}

//  Use the taylor series. 
//asin(x) == x + x^3/(2*3)+(1*3)x^5/(2*4*5)+(1*3*5)x^7/(2*4*6*7)....
static interval<double> _intervalasindouble(const double x)
	{
	int sign = 1;
	double xn = x;
	interval<double> zn, znsq, sq2, delta, sum, uc, lc, c1(1), c2(2);
	int i, k;
	// uc, lc;
	if (fabs(x)>1)
		throw interval<double>::domain_error(); 
	//	Use argument reduction via the identity 
	//	arcsin(x)=2arcsin(x/(sqrt(2)*sqrt(1+sqrt(1-x*x)))	
	//  The argument reduction is used to reduced the number of taylor iteration
	//	and to minimize round off erros and calculation time
	if( x<0 ) {xn = -xn; sign = -1; }
	zn = interval<double>(xn); sq2 = sqrt(c2);
	for (k = 0; zn.center()>0.5; ++k)
		{
		zn /= sq2 * sqrt( c1 + sqrt( c1 - zn * zn ) );
		}
	sum = zn; znsq = zn * zn; delta = zn;
	for (i = 3;; i += 2)
		{
		uc = interval<double>((i - 2)*(i-2)); lc = interval<double>((i - 1)*i);
		zn = ( uc * znsq ) / lc;
		delta *= zn;
		if (sum.center() + delta.center() == sum.center()) break;
		sum += delta;
		}
	if (k>0) // Reverse reduction
		sum *= interval<double>(pow(2,k));
	if (sign<0) sum = -sum;
	return sum;
	}

// Use the identity. acos(x)==PI/2-asin(x)
static interval<double> _intervalacosdouble(const double x)
	{
	interval<double> z;
	if (fabs(x)>1)
		throw interval<double>::domain_error();
	z = _intervalasindouble(x);
	z = PI * interval<double>(0.5) - z;
	return z;
	}

// Use the Taylor serie. 
// ArcTan(x) = x - x^3/3 + x^5/5 ...
static interval<double> _intervalatandouble(const double x)
	{
	int i, k;
	interval<double> zn, znsq, delta, sum;
	const interval<double> c1(1);
	// First reduce x to abs(x)< 1 to improve taylor series
	// using the identity. ArcTan(x)=2*ArcTan(x/(1+sqrt(1+x^2)))
	zn = interval<double> (x);
	for (k = 0; abs(zn.center())>0.5; ++k)
		zn /= c1+sqrt(c1+zn*zn);
	// Iterate ArcTan(x) = x - x^3/3 + x^5/5 ...
	znsq = zn * zn; sum = zn;
	for (i = 3; ; i += 2)
		{
		zn *= znsq; zn = -zn;
		delta = zn / interval<double>(i);
		if (sum.center() + delta.center() == sum.center()) break;
		sum += delta;
		}
	if (k>0) // Reverse reduction
		sum *= interval<double>(pow(2,k));
	return sum;
	}

static interval<double> _intervalatan2double(const double y, const double x)
	{
	interval<double> u, ipi=PI;
	if (x == 0)
		{ // x == 0
		if (y == 0)
			return interval<double>();
		u = ipi * interval<double>(0.5);
		if (y < 0)
			u = -u;
		}
	else  // x != 0
		{
		u = atan(interval<double>(y) / interval<double>(x));
		if (x < 0)
			if (y < 0)
				u -= ipi;
			else
				u += ipi;
		}
	return u;
	}

// Specilization for sin float_precision
//
inline interval<float_precision> sin(const interval<float_precision>& x)
	{
	float_precision l, u;

	l.assign(x.lower());  // Assign value, precision and mode
	l.mode(ROUND_DOWN);
	l = sin(l);

	u.assign(x.upper());  // Assign value, precision and mode
	u.mode(ROUND_UP);
	u = sin(u);

	return interval<float_precision>(l, u);
	}

// sin for float using managed code.
//
inline interval<float> sin(const interval<float>& x)
	{
	float lower, upper;
	interval<double> t(x.lower(), x.upper());
	t = sin(t);
	lower = tofloat(t.lower(), ROUND_DOWN);
	upper = tofloat(t.upper(), ROUND_UP );
	return interval<float>(lower, upper);
	}

// Sin for double using managed code.
//
inline interval<double> sin(const interval<double>& x)
	{
	double lower, upper;
	interval<double> tc, tb;
	// First reduce x to between 0..2*PI
	tb = interval<double>(2) * PI;
	tc = x / tb;
	tc = interval<double>(tc.lower()<0?ceil(tc.lower()):floor(tc.lower()), tc.upper()<0?ceil(tc.upper()) : floor(tc.upper()));
	tc = x - tc*tb;
	tb = _intervalsindouble(tc.lower());
	lower = tb.lower(); upper = tb.upper();
	tc = _intervalsindouble(tc.upper());
	if (lower > tc.lower()) lower = tc.lower();
	if (upper < tc.upper()) upper = tc.upper();
	return interval<double>(lower, upper);
	}

// Specilization for cos float_precision
//
inline interval<float_precision> cos(const interval<float_precision>& x)
	{
	float_precision l, u;

	l.assign(x.lower());  // Assign value, precision and mode
	l.mode(ROUND_DOWN);
	l = cos(l);

	u.assign(x.upper());  // Assign value, precision and mode
	u.mode(ROUND_UP);
	u = cos(u);

	return interval<float_precision>(l, u);
	}

// cos for float using managed code.
//
inline interval<float> cos(const interval<float>& x)
	{
	float lower, upper;
	interval<double> t(x.lower(), x.upper());
	t = cos(t);
	lower = tofloat(t.lower(), ROUND_DOWN);
	upper = tofloat(t.upper(), ROUND_UP);
	
	return interval<float>(lower, upper);
	}

// Cos for double using managed code.
//
inline interval<double> cos(const interval<double>& x)
	{
	double lower, upper;
	interval<double> tc, tb;
	// First reduce x to between 0..2*PI
	tb = interval<double>(2) * PI;
	tc = x / tb;
	tc = interval<double>(tc.lower()<0 ? ceil(tc.lower()) : floor(tc.lower()), tc.upper()<0 ? ceil(tc.upper()) : floor(tc.upper()));
	tc = x - tc*tb;
	tb = _intervalcosdouble(tc.lower());
	lower = tb.lower(); upper = tb.upper();
	tc = _intervalcosdouble(tc.upper());
	if (lower > tc.lower()) lower = tc.lower();
	if (upper < tc.upper()) upper = tc.upper();
	return interval<double>(lower, upper);
	}

// Specilization for tan float_precision
//
inline interval<float_precision> tan(const interval<float_precision>& x)
	{
	float_precision l, u;

	l.assign(x.lower());  // Assign value, precision and mode
	l.mode(ROUND_DOWN);
	l = tan(l);

	u.assign(x.upper());  // Assign value, precision and mode
	u.mode(ROUND_UP);
	u = tan(u);

	return interval<float_precision>(l, u);
	}

// tan for float using managed code.
//
inline interval<float> tan(const interval<float>& x)
	{
	float lower, upper;
	interval<double> t(x.lower(), x.upper());
	t = tan(t);
	lower = tofloat(t.lower(), ROUND_DOWN);
	upper = tofloat(t.upper(), ROUND_UP);
 
	return interval<float>(lower, upper);
	}

// Tan for double using managed code.
//
inline interval<double> tan(const interval<double>& x)
	{
	double lower, upper;
	interval<double> tb, tc;
	// First reduce x to between 0..2*PI
	tb = interval<double>(2) * PI;
	tc = x / tb;
	tc = interval<double>(tc.lower()<0 ? ceil(tc.lower()) : floor(tc.lower()), tc.upper()<0 ? ceil(tc.upper()) : floor(tc.upper()));
	tc = x - tc*tb;
	tb = _intervaltandouble(tc.lower());
	lower = tb.lower(); upper = tb.upper();
	tc = _intervaltandouble(tc.upper());
	if (lower > tc.lower()) lower = tc.lower();
	if (upper < tc.upper()) upper = tc.upper();
	return interval<double>(lower, upper);
	}

// Specilization for arctan float_precision
//
inline interval<float_precision> atan(const interval<float_precision>& x)
{
	float_precision l, u;

	l.assign(x.lower());  // Assign value, precision and mode
	l.mode(ROUND_DOWN);
	l = atan(l);

	u.assign(x.upper());  // Assign value, precision and mode
	u.mode(ROUND_UP);
	u = atan(u);

	return interval<float_precision>(l, u);
}

// Arctan for float using managed code.
//
inline interval<float> atan(const interval<float>& x)
	{
	float lower, upper;
	lower = tofloat(_intervalatandouble(x.lower()).lower(), ROUND_DOWN);
	upper = tofloat(_intervalatandouble(x.upper()).upper(), ROUND_UP);

	return interval<float>(lower, upper);
	}

// ArcTan for double using managed code.
//
inline interval<double> atan(const interval<double>& x)
	{
	double lower, upper;
	lower = _intervalatandouble(x.lower()).lower();
	upper = _intervalatandouble(x.upper()).upper();
	return interval<double>(lower, upper);
	}

inline interval<double> asin(const interval<double>& x)
	{
	double lower, upper;
	lower = _intervalasindouble(x.lower()).lower();
	upper = _intervalasindouble(x.upper()).upper();
	return interval<double>(lower, upper);
	}

inline interval<double> acos(const interval<double>& x)
	{
	double lower, upper;
	lower = _intervalacosdouble(x.lower()).lower();
	upper = _intervalacosdouble(x.upper()).upper();
	return interval<double>(lower, upper);
	}

inline interval<float> asin(const interval<float>& x)
	{
	float lower, upper;
	lower = tofloat(_intervalasindouble(x.lower()).lower(), ROUND_DOWN);
	upper = tofloat(_intervalasindouble(x.upper()).upper(), ROUND_UP);
	return interval<float>(lower, upper);
	}

// Specilization for arcsin float_precision
//
inline interval<float_precision> asin(const interval<float_precision>& x)
	{
	float_precision l, u;

	l.assign(x.lower());  // Assign value, precision and mode
	l.mode(ROUND_DOWN);
	l = asin(l);

	u.assign(x.upper());  // Assign value, precision and mode
	u.mode(ROUND_UP);
	u = asin(u);

	return interval<float_precision>(l, u);
	}

inline interval<float> acos(const interval<float>& x)
	{
	float lower, upper;
	lower = tofloat(_intervalacosdouble(x.lower()).lower(), ROUND_DOWN);
	upper = tofloat(_intervalacosdouble(x.upper()).upper(), ROUND_UP);
	return interval<float>(lower, upper);
	}

// Specilization for arccos float_precision
//
inline interval<float_precision> acos(const interval<float_precision>& x)
	{
	float_precision l, u;

	l.assign(x.lower());  // Assign value, precision and mode
	l.mode(ROUND_DOWN);
	l = acos(l);

	u.assign(x.upper());  // Assign value, precision and mode
	u.mode(ROUND_UP);
	u = acos(u);

	return interval<float_precision>(l, u);
	}


//////////////////////////////////////////////////////////////////////////////////////
///
/// END Interval sin(), cos(), tan(), asin(), acos(), atan()
///
//////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////
///
/// Interval sinh(), cosh(), tanh(), asinh(), acosh(), atanh()
///
//////////////////////////////////////////////////////////////////////////////////////

// Specilization for sinh float_precision
//
// Use the identity. sinh(x)=0.5*(exp(x)-1/exp(x))
inline interval<float_precision> sinh(const interval<float_precision>& x)
	{
	float_precision l, u;

	l.assign(x.lower());  // Assign value, precision and mode
	l.mode(ROUND_DOWN);
	l = sinh(l);

	u.assign(x.upper());  // Assign value, precision and mode
	u.mode(ROUND_UP);
	u = sinh(u);

	return interval<float_precision>(l, u);
	}

// sinh for float using managed code.
//
// Use the identity. sinh(x)=0.5*(exp(x)-1/exp(x))
inline interval<float> sinh(const interval<float>& x)
	{// do the calculation in double and convert to float at the end.
	interval<double> d(x);
	float lower, upper;
	d = sinh(d);
	lower = tofloat(d.lower(), ROUND_DOWN);
	upper = tofloat(d.upper(), ROUND_UP);
	return interval<float>(lower,upper);
	}

// Sinh for double using managed code.
//
// Use the identity. sinh(x)=0.5*(exp(x)-1/exp(x))
inline interval<double> sinh(const interval<double>& x)
	{
	const interval<double> c1(1), c05(0.5);
	interval<double> e = exp(x);
	return c05 * (e - c1 / e);
	}

// Specilization for cosh float_precision
//
// Use the identity. cosh(x)=0.5*(exp(x)+1/exp(x))
inline interval<float_precision> cosh(const interval<float_precision>& x)
	{
	float_precision l, u;

	l.assign(x.lower());  // Assign value, precision and mode
	l.mode(ROUND_DOWN);
	l = cosh(l);

	u.assign(x.upper());  // Assign value, precision and mode
	u.mode(ROUND_UP);
	u = cosh(u);

	return interval<float_precision>(l, u);
	}

// sinh for float using managed code.
//
// Use the identity. cosh(x)=0.5*(exp(x)+1/exp(x))
inline interval<float> cosh(const interval<float>& x)
	{// do the calculation in double and convert to float at the end.
	interval<double> d(x);
	float lower, upper;
	d = cosh(d);
	lower = tofloat(d.lower(), ROUND_DOWN);
	upper = tofloat(d.upper(), ROUND_UP);
	return interval<float>(lower, upper);
	}

// Sinh for double using managed code.
//
// Use the identity. cosh(x)=0.5*(exp(x)+1/exp(x))
inline interval<double> cosh(const interval<double>& x)
	{
	const interval<double> c1(1), c05(0.5);
	interval<double> e = exp(x);
	return c05 * (e + c1 / e);
	}


// Specilization for tanh float_precision
//
// Use the identity. tanh(x)=(exp(x)-1/exp(x))/(exp(x)+1/exp(x))=(exp(x)^2-1)/(exp(x)^2+1)
inline interval<float_precision> tanh(const interval<float_precision>& x)
	{
	float_precision l, u;

	l.assign(x.lower());  // Assign value, precision and mode
	l.mode(ROUND_DOWN);
	l = tanh(l);

	u.assign(x.upper());  // Assign value, precision and mode
	u.mode(ROUND_UP);
	u = tanh(u);

	return interval<float_precision>(l, u);
	}

// tanh for float using managed code.
//
// Use the identity. tanh(x)=(exp(x)-1/exp(x))/(exp(x)+1/exp(x))=(exp(x)^2-1)/(exp(x)^2+1)
inline interval<float> tanh(const interval<float>& x)
	{// do the calculation in double and convert to float at the end.
	interval<double> d(x);
	float lower, upper;
	d = tanh(d);
	lower = tofloat(d.lower(), ROUND_DOWN);
	upper = tofloat(d.upper(), ROUND_UP);
	return interval<float>(lower, upper);
	}

// tanh for double using managed code.
//
// Use the identity. tanh(x)=(exp(x)-1/exp(x))/(exp(x)+1/exp(x))=(exp(x)^2-1)/(exp(x)^2+1)
inline interval<double> tanh(const interval<double>& x)
	{
	const interval<double> c1(1);
	interval<double> e = exp(x);
	e *= e;
	return ( e - c1 ) / ( e + c1 );
	}


// Specilization for asinh float_precision
//
// Use the identity. asinh(x)=Ln(x)+sqrt(x^2+1)
inline interval<float_precision> asinh(const interval<float_precision>& x)
	{
	float_precision l, u;

	l.assign(x.lower());  // Assign value, precision and mode
	l.mode(ROUND_DOWN);
	l = asinh(l);

	u.assign(x.upper());  // Assign value, precision and mode
	u.mode(ROUND_UP);
	u = asinh(u);

	return interval<float_precision>(l, u);
	}

// asinh for float using managed code.
//
// Use the identity. asinh(x)=Ln(x)+sqrt(x^2+1)
inline interval<float> asinh(const interval<float>& x)
	{// do the calculation in double and convert to float at the end.
	interval<double> d(x);
	float lower, upper;
	d = asinh(d);
	lower = tofloat(d.lower(), ROUND_DOWN);
	upper = tofloat(d.upper(), ROUND_UP);
	return interval<float>(lower, upper);
	}

// aSinh for double using managed code.
//
// Use the identity. asinh(x)=Ln(x+sqrt(x^2+1))
inline interval<double> asinh(const interval<double>& x)
	{
	const interval<double> c1(1);
	return log( x  + sqrt( x*x + c1 ) );
	}

// Specilization for acosh float_precision
//
// Use the identity. acosh(x)=Ln(x)+sqrt(x^2-1)
inline interval<float_precision> acosh(const interval<float_precision>& x)
	{
	float_precision l, u;

	l.assign(x.lower());  // Assign value, precision and mode
	l.mode(ROUND_DOWN);
	l = acosh(l);

	u.assign(x.upper());  // Assign value, precision and mode
	u.mode(ROUND_UP);
	u = acosh(u);

	return interval<float_precision>(l, u);
	}

// acosh for float using managed code.
//
// Use the identity. acosh(x)=Ln(x+sqrt(x^2-1))
inline interval<float> acosh(const interval<float>& x)
	{// do the calculation in double and convert to float at the end.
	interval<double> d(x);
	float lower, upper;
	d = acosh(d);
	lower = tofloat(d.lower(), ROUND_DOWN);
	upper = tofloat(d.upper(), ROUND_UP);
	return interval<float>(lower, upper);
	}

// acosh for double using managed code.
//
// Use the identity. acosh(x)=Ln(x)+sqrt(x^2-1)
inline interval<double> acosh(const interval<double>& x)
	{
	const interval<double> c1(1);
	if (x.lower() < 1 )
		throw interval<double>::domain_error();
	return log(x + sqrt( x*x - c1 ) );
	}


// Specilization for atanh float_precision
//
// Use the identity. atanh(x)=0.5*Ln((1+x)/(1-x))
inline interval<float_precision> atanh(const interval<float_precision>& x)
	{
	float_precision l, u;

	l.assign(x.lower());  // Assign value, precision and mode
	l.mode(ROUND_DOWN);
	l = atanh(l);

	u.assign(x.upper());  // Assign value, precision and mode
	u.mode(ROUND_UP);
	u = atanh(u);

	return interval<float_precision>(l, u);
	}

// atanh for float using managed code.
//
// Use the identity. atanh(x)=0.5*Ln((1+x)/(1-x))
inline interval<float> atanh(const interval<float>& x)
	{// do the calculation in double and convert to float at the end.
	interval<double> d(x);
	float lower, upper;
	d = atanh(d);
	lower = tofloat(d.lower(), ROUND_DOWN);
	upper = tofloat(d.upper(), ROUND_UP);
	return interval<float>(lower, upper);
	}

// atanh for double using managed code.
//
// Use the identity. atanh(x)=0.5*Ln((1+x)/(1-x))
inline interval<double> atanh(const interval<double>& x)
	{
	const interval<double> c1(1), c05(0.5);
	if(x.lower() <=-1 || x.upper() >= 1 )
		throw interval<double>::domain_error();
	return c05*log( (c1+x) / (c1-x) );
	}
//////////////////////////////////////////////////////////////////////////////////////
///
/// END Interval sinh(), cosh(), tanh(), asinh(), acosh(), atanh()
///
//////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
//
// Low Level Interval arithmetic
//
//////////////////////////////////////////////////////////////////////////
static void split(const double &a, double &x, double &y)
	{
	double tmp;
	static const double sigma = ldexp(1.0, 27) + 1;
	tmp = a* sigma;
	x = tmp - (tmp - a);
	y = a - x;
	}

static double succ(const double &x)
	{
	static const double th1 = ldexp(1.0, -969);
	static const double th2 = ldexp(1.0, -1021);
	static const double c1 = ldexp(1.0, -53) + ldexp(1.0, -105);
	static const double c2 = ldexp(1.0, -1074);
	static const double c3 = ldexp(1.0, 53);
	static const double c4 = ldexp(1.0, -53);
	double a, c, e;
	a = fabs(x);
	if (a >= th1) return x + a*c1;
	if (a < th2) return x + c2;
	c = c3*x;
	e = c1*fabs(c);
	return (c + e) * c4;
	}

static double pred(const double &x)
	{
	static const double th1 = ldexp(1.0, -969);
	static const double th2 = ldexp(1.0, -1021);
	static const double c1 = ldexp(1.0, -53) + ldexp(1.0, -105);
	static const double c2 = ldexp(1.0, -1074);
	static const double c3 = ldexp(1.0, 53);
	static const double c4 = ldexp(1.0, -53);
	double a, c, e;
	a = fabs(x);
	if (a >= th1) return x - a*c1;
	if (a < th2) return x - c2;
	c = c3*x;
	e = c1*fabs(c);
	return (c - e) * c4;
	}

static void twosum(const double &a, const double &b, double &x, double &y)
	{
	double tmp;
	x = a + b;
	if (fabs(a) > fabs(b))
		{
		tmp = x - a;
		y = b - tmp;
		}
	else
		{
		tmp = x - b;
		y = a - tmp;
		}
	}

static void twoproduct(const double &a, const double b, double &x, double &y)
	{
	static const double th = ldexp(1.0, 996);
	static const double c1 = ldexp(1.0, -28);
	static const double c2 = ldexp(1.0, 28);
	static const double th2 = ldexp(1.0, 1023);
	double na, nb, a1, a2, b1, b2;
	x = a*b;
	if (fabs(a) > th)
		{
		na = a * c1;
		nb = b * c2;
		}
	else
		if (fabs(b) > th)
			{
			na = a * c2;
			nb = b * c1;
			}
		else
			{
			na = a;
			nb = b;
			}
	split(na, a1, a2);
	split(nb, b1, b2);
	if (fabs(x) > th2)
		y = a2 * b2 - ((((x*0.5) - (a1*0.5)*b1)*2.0 - a2*b1) - a1*b2);
	else
		y = a2 * b2 - (((x - a1*b1) - a2*b1) - a1*b2);
	}


inline float tofloat(const double &d, const enum round_mode rm)
	{
	float f = (float)d; // round to closets

	if (f != d)
		{
		if (rm == ROUND_UP)
			{
			if (f < d  )
				f = nextafterf( f, std::numeric_limits<float>::infinity());
			}
		if (rm == ROUND_DOWN)
			{//
			if ( f > d )
				f = nextafterf( f, -std::numeric_limits<float>::infinity());
			}
		}
	return f;
	/*if (lower > c.lower() )
	lower -= lower * 0.5f * FLT_EPSILON;
	if (upper < c.upper() )
	upper += upper * 0.5f * FLT_EPSILON;
	*/
	}

inline double add_down(const double &x, const double &y)
	{
	double r, r2;
	twosum(x, y, r, r2);
	if (r == std::numeric_limits<double>::infinity())
	{
		if (x == std::numeric_limits<double>::infinity() ||
			y == std::numeric_limits<double>::infinity())
			return r;
		else
			return std::numeric_limits<double>::max();
	}
	else
		if (r == -std::numeric_limits<double>::infinity())
			return r;
	if (r2 < 0)
		return pred(r);
	return r;
	}

inline double add_up(const double &x, const double &y)
	{
	double r, r2;
	twosum(x, y, r, r2);
	if (r == std::numeric_limits<double>::infinity())
		return r;
	else
		if (r == -std::numeric_limits<double>::infinity())
			if (x == -std::numeric_limits<double>::infinity() ||
				y == -std::numeric_limits<double>::infinity())
				return r;
			else
				return -std::numeric_limits<double>::max();
	if (r2 > 0)
		return succ(r);
	return r;
	}

inline double sub_down(const double &x, const double &y)
	{
	double r, r2;
	twosum(x, -y, r, r2);
	if (r == std::numeric_limits<double>::infinity())
		{
		if (x == std::numeric_limits<double>::infinity() ||
			y == std::numeric_limits<double>::infinity())
			return r;
		else
			return std::numeric_limits<double>::max();
		}
	else
		if (r == -std::numeric_limits<double>::infinity())
			return r;
	if (r2 < 0)
		return pred(r);
	return r;
	}

inline double sub_up(const double &x, const double &y)
	{
	double r, r2;
	twosum(x, -y, r, r2);
	if (r == std::numeric_limits<double>::infinity())
		return r;
	else
		if (r == -std::numeric_limits<double>::infinity())
			if (x == -std::numeric_limits<double>::infinity() ||
				y == -std::numeric_limits<double>::infinity())
				return r;
			else
				return -std::numeric_limits<double>::max();
	if (r2 > 0)
		return succ(r);
	return r;
	}

inline double mul_down(const double &x, const double &y)
	{
	double r, r2;
	double s, s2, t;
	static const double th = ldexp(1.0, -969);
	static const double c = ldexp(1.0, 537);
	twoproduct(x, y, r, r2);

	if (r == std::numeric_limits<double>::infinity())
		{
		if (fabs(x) == std::numeric_limits<double>::infinity() ||
			fabs(y) == std::numeric_limits<double>::infinity())
			return r;
		else
			return std::numeric_limits<double>::max();
		}
	else
		if (r == -std::numeric_limits<double>::infinity())
			return r;
	if (fabs(r) >= th)
		{
		if (r2 < 0)
			return pred(r);
		return r;
		}
	else
		{
		twoproduct(x*c, y*c, s, s2);
		t = (r*c)*c;
		if (t > s || (t == s && s2 < 0.0))
			return pred(r);
		}
	return r;
	}

inline double mul_up(const double &x, const double &y)
	{
	double r, r2;
	double s, s2, t;
	static const double th = ldexp(1.0, -969);
	static const double c = ldexp(1.0, 537);
	twoproduct(x, y, r, r2);

	if (r == std::numeric_limits<double>::infinity())
		return r;
	else if (r == -std::numeric_limits<double>::infinity())
		{
		if (fabs(x) == std::numeric_limits<double>::infinity() ||
			fabs(y) == std::numeric_limits<double>::infinity())
			return r;
		else
			return -std::numeric_limits<double>::max();
		}
	if (fabs(r) >= th)
		{
		if (r2 > 0)
			return succ(r);
		return r;
		}
	else
		{
		twoproduct(x*c, y*c, s, s2);
		t = (r*c)*c;
		if (t < s || (t == s && s2 > 0.0))
			return succ(r);
		}
	return r;
	}

inline double div_up(const double &x, const double &y)
	{
	double r, r2;
	double xn, yn, d;
	static const double th1 = ldexp(1.0, -969);
	static const double th2 = ldexp(1.0, 918);
	static const double c1 = ldexp(1.0, 105);
	static const double c2 = ldexp(1.0, -1074);
	if (x == 0 || y == 0 || fabs(x) == std::numeric_limits<double>::infinity() ||
		fabs(y) == std::numeric_limits<double>::infinity() || x != x || y != y)
		return x / y;
	if (y < 0)
		{
		xn = -x;
		yn = -y;
		}
	else
		{
		xn = x;
		yn = y;
		}	
	if (fabs(xn)<th1)
		{
		if (fabs(yn) < th2)
			{
			xn *= c1;
			yn *= c1;
			}
		else
			{
			if (xn < 0) return 0;
			else return c2;
			}
		}
	d = xn / yn;
	if (d == std::numeric_limits<double>::infinity())
		return d;
	else if (d == -std::numeric_limits<double>::infinity())
		return -(std::numeric_limits<double>::max());
	twoproduct(d, yn, r, r2);
	if (r < xn || (r == xn) && r2 < 0)
		return succ(d);
	return d;
	}

inline double div_down(const double &x, const double &y)
	{
	double r, r2;
	double xn, yn, d;
	static const double th1 = ldexp(1.0, -969);
	static const double th2 = ldexp(1.0, 918);
	static const double c1 = ldexp(1.0, 105);
	static const double c2 = ldexp(1.0, -1074);
	if (x == 0 || y == 0 || fabs(x) == std::numeric_limits<double>::infinity() ||
		fabs(y) == std::numeric_limits<double>::infinity() || x != x || y != y)
		return x / y;
	if (y < 0)
		{
		xn = -x;
		yn = -y;
		}
	else
		{
		xn = x;
		yn = y;
		}
	if (fabs(xn)<th1)
		{
		if (fabs(yn) < th2)
			{
			xn *= c1;
			yn *= c1;
			}
		else
			{
			if (xn < 0) return -c2;
			else return 0;
			}
		}
	d = xn / yn;
	if (d == std::numeric_limits<double>::infinity())
		return std::numeric_limits<double>::max();
	else if (d == -std::numeric_limits<double>::infinity())
		return d;
	twoproduct(d, yn, r, r2);
	if (r > xn || (r == xn) && r2 > 0)
		return pred(d);
	return d;
	}

inline double sqrt_up(const double &x)
	{
	double r, r2, d;
	static const double th1 = ldexp(1.0, -969);
	static const double c1 = ldexp(1.0, 106);
	static const double c2 = ldexp(1.0, 53);
	d = sqrt(x);
	if (x < th1)
		{
		double d2, x2;
		x2 = x*c1;
		d2 = d*c2;
		twoproduct(d2, d2, r, r2);
		if (r < x2 || (r == x2) && r2 < 0)
			return succ(d);
		return d;
		}
	twoproduct(d, d, r, r2);
	if (r < x || (r == x) && r2 < 0)
		return succ(d);
	return d;
	}

inline double sqrt_down(const double &x)
	{
	double r, r2, d;
	static const double th1 = ldexp(1.0, -969);
	static const double c1 = ldexp(1.0, 106);
	static const double c2 = ldexp(1.0, 53);
	d = sqrt(x);
	if (x < th1)
		{
		double d2, x2;
		x2 = x*c1;
		d2 = d*c2;
		twoproduct(d2, d2, r, r2);
		if (r > x2 || (r == x2) && r2 > 0)
			return pred(d);
		return d;
		}
	twoproduct(d, d, r, r2);
	if (r > x || (r == x) && r2 > 0)
		return pred(d);
	return d;
	}

//////////////////////////////////////////////////////////////////////////
//
// END Low Level Interval arithmetic
//
//////////////////////////////////////////////////////////////////////////
#endif
