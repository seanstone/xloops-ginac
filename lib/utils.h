/** @file utils.h
 *
 *  Definition of utility functions that are not supposed to be visible
 *  outside the library. */

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

#ifndef __XLOOPS_UTILS_H__
#define __XLOOPS_UTILS_H__

#include <ginac/ginac.h>
using namespace GiNaC;

namespace xloops {


/** Return the coefficient of P1^m * P2^n * P3^o * P4^r in A. */
extern ex Pcollect(const ex &A, int m, int n, int o, int r, const symbol &P1, const symbol &P2, const symbol &P3, const symbol &P4);

/** Return coefficient of the transformation
 *    ln^i * lm^j * lr^k * g^2(l) -> CoeffDim(i,j,k,eps) * l^(i+j+k) * g^2(l) */
extern ex CoeffDim(int i, int j, int k, const ex &eps);

/** Reduce general one-loop one-point function to scalar functions. */
extern ex reduce_OneLoop1Pt(int p, const ex &m, int t, const ex &rho, const ex &eps);

/** Reduce general one-loop two-point function to scalar functions. */
extern ex reduce_OneLoop2Pt(int p0, int p1, const ex &q, const ex &m1, const ex &m2, int t1, int t2f, const ex &rho, const ex &eps);

/** Reduce general one-loop three-point function to scalar functions. */
extern ex reduce_OneLoop3Pt(int p0, int p1, int p2, const ex &q10, const ex &q20, const ex &q21, const ex &m1, const ex &m2, const ex &m3, int t1, int t2, int t3, const ex &rho, const ex &eps);

/** Reduce one-loop two-point function with two parallel dimensions, if possible. */
extern ex reduce_S322Pt(int p0, int p1, int p2, const ex &q10, const ex &q20, const ex &q21, const ex &m1, const ex &m2, int t1, int t2, const ex &rho, const ex &eps);

} // namespace xloops

#endif // ndef __XLOOPS_UTILS_H__
