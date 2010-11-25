/** \file real2int.h
  * Source of this convertion routine is posted
  * <a href="http://www.stereopsis.com/FPU.html#convert">here</a>.
  *
  * Description of the floating point storage layout:
  * <a href="http://www.psc.edu/general/software/packages/ieee/ieee.html">IEEE-754 standart</a>.
  *
  * Basically all doubles are stored in memory in form (\b S is sign, \b E is exponent and \b F is fraction):
  *
  <pre>
    S EEEEEEEEEEE FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    0 1        11 12                                                63
  </pre>
  * Number itself is \f$ Sign*2^{Exponent}*1.mantissa \f$ and FFFF..FFF is mantissa.
  *
  * Goal is to round toward zero.
  */

#ifndef MC_REAL2INT_HEADER
#define MC_REAL2INT_HEADER

#undef DEFAULT_CONVERSION

typedef unsigned long uint32;
typedef long int32;

#if BigEndian_
  #define MC_IEXP  0
  #define MC_IMAN  1
#else
  #define MC_IEXP  1
  #define MC_IMAN  0
#endif

// ---------------------------------------------------------------------------
/// double -> int.
// ---------------------------------------------------------------------------
static int32
double2int (double val)
{
#if DEFAULT_CONVERSION
  return (int)(val + 1000) - 1000;
#else
  val = val + 68719476736.0*1.5;
  return ((int32*) &val)[MC_IMAN] >> 16;
#endif
}

#endif
