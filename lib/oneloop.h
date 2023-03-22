/** @file oneloop.h
 *
 *  Definition of functions for one-loop Feynman integrals. */

/*
 *  xloops Copyright (C) 1997,2000-2001 Johannes Gutenberg University Mainz,
 *  Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef __XLOOPS_ONELOOP_H__
#define __XLOOPS_ONELOOP_H__

#include <ginac/ginac.h>

namespace xloops {

using GiNaC::ex;


/*
 *  One-loop one-point function
 *
 *     /\      p
 *     |  D   l
 *     | d l ----
 *     |       t
 *    \/      P
 *
 *  P = l^2 - m^2 + I*rho
 *
 *  p  : tensor degree
 *  m  : mass
 *  t  : power of denominator
 *  rho: small imaginary part of propagator
 */

DECLARE_FUNCTION_4P(OneLoop1Pt)
DECLARE_FUNCTION_2P(Scalar1Pt)


/*
 *  One-loop two-point function
 *
 *     /\      p0 p1
 *     |      l  l
 *     |  D    0  +
 *     | d l --------
 *     |       t1  t2
 *     |      P   P
 *    \/       1   2
 *
 *  P1 = (l0 + q)^2 - l+^2 - m1^2 + I*rho
 *  P2 = l0^2 - l+^2 - m2^2 + I*rho
 *
 *  p0,p1: tensor degrees
 *  q    : external momentum
 *  m1,m2: masses
 *  t1,t2: powers of denominators
 *  rho  : small imaginary part of propagators
 */

DECLARE_FUNCTION_8P(OneLoop2Pt)
DECLARE_FUNCTION_6P(Scalar2Pt)


/*
 *  One-loop three-point function
 *
 *     /\      p0 p1 p2
 *     |      l  l  l
 *     |  D    0  1  +
 *     | d l -------------
 *     |       t1  t2  t3
 *     |      P   P   P
 *    \/       1   2   3
 *
 *  P1 = l0^2 - l1^2 - l+^2 + 2*l0*q10 + q10^2 - m1^2 + I*rho
 *  P2 = l0^2 - l1^2 - l+^2 + 2*l0*q20 - 2*l1*q21 + q20^2 - q21^2 - m2^2 + I*rho
 *  P3 = l0^2 - l1^2 - l+^2 - m3^2 + I*rho
 *
 *  p0,p1,p2   : tensor degrees
 *  q10,q20,q21: external momenta
 *  m1,2,3     : masses
 *  t1,t2,t3   : powers of denominators
 *  rho        : small imaginary part of propagators
 */

DECLARE_FUNCTION_13P(OneLoop3Pt)
DECLARE_FUNCTION_10P(Scalar3Pt)


/*
 *  One-loop two-point function with 2 parallel dimensions (used internally)
 *
 *     /\      p0 p1 p2
 *     |      l  l  l
 *     |  D    0  1  +
 *     | d l -----------
 *     |        t1  t2
 *     |       P   P
 *    \/        1   2
 *
 *  P1 = (l0 + q10)^2 - l1^2 - l+^2 - m1^2 + I*rho
 *  P2 = (l0 + q20)^2 - (l1 + q21)^2 - l+^2 - m2^2 + I*rho
 *
 *  p0,p1,p2   : tensor degrees
 *  q10,q20,q21: external momentum
 *  m1,m2      : masses
 *  t1,t2      : powers of denominators
 *  rho        : small imaginary part of propagators
 */

DECLARE_FUNCTION_11P(S322Pt)


/** Reduce all OneLoopNPt() functions to scalar functions, if possible */
extern ex tensor_reduction(const ex &e, const ex &eps);


/** Return tensor decomposition of one-point function of tensor rank n, and
 *  the list "Cl" of defining equations for the coefficients (whoch correspond
 *  to those of the Passarino-Veltman procedure) in a subs()able format. */
extern ex OneLoopTens1Pt(unsigned n, const ex & m, const ex & t, const ex & rho, const ex & eps, GiNaC::lst & Cl);

/** Return tensor decomposition of two-point function of tensor rank n, and
 *  the list "Cl" of defining equations for the coefficients (whoch correspond
 *  to those of the Passarino-Veltman procedure) in a subs()able format. */
extern ex OneLoopTens2Pt(unsigned n, const ex & q, const ex & m1, const ex & m2, const ex & t1, const ex & t2, const ex & rho, const ex & eps, GiNaC::lst & Cl);

/** Return tensor decomposition of three-point function of tensor rank n, and
 *  the list "Cl" of defining equations for the coefficients (whoch correspond
 *  to those of the Passarino-Veltman procedure) in a subs()able format. */
extern ex OneLoopTens3Pt(unsigned n, const ex & q10, const ex & q20, const ex & q21, const ex & m1, const ex & m2, const ex & m3, const ex & t1, const ex & t2, const ex & t3, const ex & rho, const ex & eps, GiNaC::lst & Cl);

 
} // namespace xloops

#endif // __XLOOPS_ONELOOP_H__
