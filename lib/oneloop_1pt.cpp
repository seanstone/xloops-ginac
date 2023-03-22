/** @file oneloop_1pt.cpp
 *
 *  Implementation of functions for one-loop one-point Feynman integrals. */

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

#include <stdexcept>

#include "oneloop.h"
#include "utils.h"

namespace xloops {


/*
 *  Reduce general one-loop one-point function to scalar functions
 */

ex reduce_OneLoop1Pt(int p,
                     const ex &m,
                     int t,
                     const ex &rho, const ex &eps)
{
	// ...
	return OneLoop1Pt(p, m, t, rho);
}


/*
 *  General one-loop one-point function
 */

static ex OneLoop1Pt_eval(const ex &p, const ex &m,
                          const ex &t, const ex &rho)
{
	// The 1-point function vanishes for odd p or zero mass m
	if (m.is_zero() || p.info(info_flags::odd))
		return 0;
	else if (p.is_zero() && t.is_equal(1))
		return Scalar1Pt(m, rho);
	else
		return OneLoop1Pt(p, m, t, rho).hold();
}

static ex OneLoop1Pt_series(const ex &p, const ex &m,
                            const ex &t, const ex &rho,
                            const relational &r, int order, unsigned options)
{
	// Check parameters for validity
	if (!p.info(info_flags::nonnegint) || !t.info(info_flags::posint))
	    throw(std::logic_error("OneLoop1Pt_series: p<0 or t<=0"));
        
	// Check whether we can do the expansion and get parameters
	if (!r.rhs().is_zero())
		throw(std::logic_error("OneLoop1Pt_series: don't know the expansion at eps!=0"));
	const ex eps = r.lhs();
	int t_int = ex_to<numeric>(t).to_int();
	
	// Assemble factors of integral
	ex p_half = p / 2;
	ex M = pow(m, 2) - I*rho;
	ex c = pow(Pi,2)*pow(2*Pi, -2*eps)*I*pow(-1, p_half - t_int)*pow(4*Pi,eps)*pow(M, -eps)*pow(M, 2+p_half-t_int)
		*tgamma(2+p_half-eps)*tgamma(t_int-p_half-2+eps)/tgamma(2-eps)/tgamma(t_int);
	
	/*
	  for (int i=0; i<t_int-1; i++)
	  c *= 1 + p_half - eps - i;
	  c *= I * pow(Pi, 2) * pow(Pi, -eps) * tgamma(eps) / (1 - eps) * pow(M, -eps)
	  * pow(m, 4 + p - t * 2) / factorial(t_int - 1);
	  */
	// Expand into series in eps
	return c.series(r, order);
}

REGISTER_FUNCTION(OneLoop1Pt, eval_func(OneLoop1Pt_eval).
                              series_func(OneLoop1Pt_series));


/*
 *  Scalar one-loop one-point function (p = t = 0)
 */

static ex Scalar1Pt_eval(const ex &m, const ex &rho)
{
	// The scalar 1-point function vanishes zero mass m
	if (m.is_zero())
		return 0;
	else
		return Scalar1Pt(m, rho).hold();
}

static ex Scalar1Pt_series(const ex &m, const ex &rho,
                           const relational &r, int order, unsigned options)
{
	// Check whether we can do the expansion and get parameters
	if (!r.rhs().is_zero())
		throw(std::logic_error("Scalar1Pt_series: don't know the expansion at eps!=0"));
	const ex eps = r.lhs();
	
	// Assemble factors of integral
	ex M = pow(m, 2) - I*rho;
	ex c = I * pow(Pi, 2) * pow(Pi, -eps) * tgamma(eps) / (1 - eps) * pow(M, -eps) * pow(m, 2);
	
	// Expand into series in eps
	return c.series(r, order);
}

REGISTER_FUNCTION(Scalar1Pt, eval_func(Scalar1Pt_eval).
                             series_func(Scalar1Pt_series));


class reduce_map_function : public map_function {
	ex eps;
public:
	reduce_map_function(const ex &eps_) : eps(eps_) {}

	ex operator()(const ex & e)
	{
		if (is_ex_the_function(e, OneLoop1Pt))
			return reduce_OneLoop1Pt(
				ex_to<numeric>(e.op(0)).to_int(),
				e.op(1),
				ex_to<numeric>(e.op(2)).to_int(),
				e.op(3),
				eps
			);
		else if (is_ex_the_function(e, OneLoop2Pt))
			return reduce_OneLoop2Pt(
				ex_to<numeric>(e.op(0)).to_int(), ex_to<numeric>(e.op(1)).to_int(),
				e.op(2),
				e.op(3), e.op(4),
				ex_to<numeric>(e.op(5)).to_int(), ex_to<numeric>(e.op(6)).to_int(),
				e.op(7),
				eps
			);
		else if (is_ex_the_function(e, OneLoop3Pt))
			return reduce_OneLoop3Pt(
				ex_to<numeric>(e.op(0)).to_int(), ex_to<numeric>(e.op(1)).to_int(), ex_to<numeric>(e.op(2)).to_int(),
				e.op(3), e.op(4), e.op(5),
				e.op(6), e.op(7), e.op(8),
				ex_to<numeric>(e.op(9)).to_int(), ex_to<numeric>(e.op(10)).to_int(), ex_to<numeric>(e.op(11)).to_int(),
				e.op(12),
				eps
			);
		else if (is_ex_the_function(e, S322Pt))
			return reduce_S322Pt(
				ex_to<numeric>(e.op(0)).to_int(), ex_to<numeric>(e.op(1)).to_int(), ex_to<numeric>(e.op(2)).to_int(),
				e.op(3), e.op(4), e.op(5),
				e.op(6), e.op(7),
				ex_to<numeric>(e.op(8)).to_int(), ex_to<numeric>(e.op(9)).to_int(),
				e.op(10),
				eps
			);
		else
			return e.map(*this);
	}
};

// Reduce all OneLoopNPt() functions to scalar functions
ex tensor_reduction(const ex &e, const ex &eps)
{
	// Perform reduction
	reduce_map_function map_reduce(eps);
	ex r = map_reduce(e).expand();

	// Collect by ScalarNPt() functions
	ex X1 = wild(1), X2 = wild(2), X3 = wild(3), X4 = wild(4), X5 = wild(5), X6 = wild(6),
	   X7 = wild(7), X8 = wild(8), X9 = wild(9), X10 = wild(10), X11 = wild(11);
	lst fcns;
	r.find(Scalar1Pt(X1, X2), fcns);
	r.find(Scalar2Pt(X1, X2, X3, X4, X5, X6), fcns);
	r.find(Scalar3Pt(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10), fcns);
	r.find(S322Pt(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11), fcns);
	return r.collect(fcns);
}


} // namespace xloops
