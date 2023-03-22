/** @file twoloop.h
 *
 *  Definition of functions for two-loop Feynman integrals. */

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

#ifndef __XLOOPS_TWOLOOP_H__
#define __XLOOPS_TWOLOOP_H__

#include <ginac/ginac.h>
#include "oneloop.h"

namespace xloops {

using GiNaC::ex;
	
/*
 *  Two-loop two-point function
 *
 *     /\    /\      p0  p1  r0  r1  t
 *     |     |      l   l   k   k   z
 *     |  D  |  D    0   +   0   +
 *     | d l | d k --------------------
 *     |     |       t1  t2  t3  t4  t5    
 *     |     |      P   P   P   P   P    
 *    \/    \/       1   2   3   4   5  
 *
 *  P1 = (l + q)^2 - m1^2 + I*rho
 *  P2 = l^2 - m2^2 + I*rho
 *  P3 = (l + k)^2 - m3^2 +I*rho
 *  P4 = (k - q)^2 - m4^2 +I*rho
 *  P5 = k^2 - m5^2 + I*rho
 *
 *  p0, p1, r0, r1, t: tensor degrees
 *  q    : external momentum
 *  m1, m2, m3, m4, m5 : masses
 *  t1, t2, t3, t4, t5 : powers of denominators
 *  rho  : small imaginary part of propagators
 *
 *  SYNOPSIS:
 *  TwoLoop2Pt(const ex &p0, const ex &p1, const ex &r0, const ex &r1,
 *             const ex &t, dynamic dyna, int topo, const ex & rho)
 */

DECLARE_FUNCTION_10P(TwoLoop2Pt)

	
} // namespace xloops

#endif // __XLOOPS_TWOLOOP_H__
