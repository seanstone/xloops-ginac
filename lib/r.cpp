/** @file r.cpp
 *
 *  Implementation of necessary R functions. */

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

#include "r.h"

using namespace GiNaC;

namespace xloops {


/*
 *  Helper function: extract numeric term from sum
 */

static ex NumericTerm(const ex &x)
{
	if (is_exactly_a<numeric>(x))
		return x;
	else if (is_a<add>(x)) {
		for (int i=x.nops()-1; i>=0; i--)
			if (is_exactly_a<numeric>(x.op(i)))
				return x.op(i);
	}
	return 0;
}


/*
 *  Two-dimensional R function
 *
 *  t:    index
 *  b1,2: parameters
 *  z1,2: arguments
 */

static ex R2_eval(const ex &t,
                  const ex &b1, const ex &b2,
                  const ex &z1, const ex &z2,
                  const ex &rho )
{
	// R function is zero when index is zero
	if (t.is_zero())
		return 0;

	// Equal arguments or vanishing parameters reduce the R function to the power function
	if (z1.is_equal(z2) || b2.is_zero())
		return pow(z1, t);
	else if (b1.is_zero())
		return pow(z2, t);

	// Sum of parameters
	ex beta = b1 + b2;

	// Feynman parameter formula
	if (t.is_equal(-beta))
		return pow(z1, -b1) * pow(z2, -b2);

	// Directly calculable cases
	if (b1 == 1 && b2 == 1) {
		if (t == -1)
			return (log(z1) - log(z2)) / (z1 - z2);
		else
			return (pow(z1, 1+t) - pow(z2, 1+t)) / ((1+t) * (z1 - z2));
	}

	// Numeric terms in index and parameters
	ex t_num = NumericTerm(t);
	ex b1_num = NumericTerm(b1);

	// Decrease second parameter
	if (b2 > 1 && !beta.is_zero())
		return (beta - 1) / ((b2 - 1) * (z1 - z2)) * (z2 * R2(t, b1, b2-1, z1, z2, rho) - R2(t+1, b1, b2-1, z1, z2, rho));

	// Decrease first parameter
	if (b1_num > 0 && b2.is_equal(1) && !(t.is_equal(-1)))
		return b1 / ((t+1) * (z2 - z1)) * (R2(t+1, b1-1, 1, z1, z2, rho) - pow(z1, t+1));

	// Increase first parameter
	if (b1_num <= -1 && !beta.is_zero())
		return ((beta + t) * R2(t, b1+1, b2, z1, z2, rho) - t * z1 * R2(t-1, b1+1, b2, z1, z2, rho)) / beta;

	// Decrease index
	if (t_num >= 1 && b2.is_equal(1))
		return (b1 * pow(z1, t) + t * z2 * R2(t-1, b1, b2, z1, z2, rho)) / (b1 + t);
	if (t_num > 1)
		return -((t - 1) * z1 * z2 * R2(t-2, b1, b2, z1, z2, rho) + ((1 - t - b1) * z1 + (1 - t - b2) * z2) * R2(t-1, b1, b2, z1, z2, rho)) / (t + b1 + b2 - 1);

	// Increase index
	if (t_num < 0 && !(t.is_equal(-1))){
		return -1 / ((t+1) * z1 * z2) * ((beta+t+1) * R2(t+2, b1, b2, z1, z2, rho) - ((t+b1+1)*z1 + (t+b2+1)*z2) * R2(t+1, b1, b2, z1, z2, rho));
	}
	return R2(t, b1, b2, z1, z2, rho).hold();
}

static ex R2_deriv(const ex &t,
                   const ex &b1, const ex &b2,
                   const ex &z1, const ex &z2,
                   const ex &rho,
                   unsigned deriv_param)
{
	if (deriv_param == 3)		// d/dz1
		return t * b1 / (b1 + b2) * R2(t-1, b1+1, b2, z1, z2, rho);
	else if (deriv_param == 4)	// d/dz2
		return t * b2 / (b1 + b2) * R2(t-1, b1, b2+1, z1, z2, rho);
	else
		throw(std::logic_error("don't know the derivative of this R2 function"));
}

static ex R2_series(const ex &t,
                    const ex &b1, const ex &b2,
                    const ex &z1, const ex &z2,
                    const ex &rho,
                    const relational &r,
                    int order,
                    unsigned options)
{
	const symbol &eps = ex_to<symbol>(r.lhs());
	if (!r.rhs().is_zero())
		throw(std::logic_error("don't know the series expansion of R2 function at this point"));

	// We only handle some very few special cases
	ex e;
	if (t.is_equal(-eps) && b1.is_equal(eps - numeric(1,2)) && b2.is_equal(1)) {
		if (z1 != I*rho)
			e = 1 + eps * R2ex1(z1, z2) + pow(eps, 2) * R2ex2(z1, z2) + Order(pow(eps, 3));
		else
			e = 1 + (eps + pow(eps, 2))*(-I*Pi - 2*log(2) - log(z2)) + Order(pow(eps, 3));
	} else if (t.is_equal(-eps) && b1.is_equal(eps) && b2.is_equal(1)) {
		ex lambda = log(z2) * (eta(z1-z2, 1/(1-z2)) - eta(z1-z2, -1/z2)) + log(1-z1/z2) * eta(z1, 1/z2);
		e = 1 - eps * log(z2) + pow(eps, 2) * (pow(log(z2), 2) / 2 + Li2(1 - z1/z2) + lambda) + Order(pow(eps, 3));
	} else if (t.info(info_flags::negint) && b1.is_equal(eps) && b2.is_equal(1 - t - eps)) {
		e = pow(z2, eps + t) * R2(-eps, -t, 1, z1, z2, rho);
	} else{
		// cout << "R2(" << t << "," << b1 << "," << b2 << "," << z1 << "," << z2 << ")" <<endl;
		throw(std::logic_error("don't know the series expansion of this particular R2 function"));
	}
	return e.series(r, order);
}

REGISTER_FUNCTION(R2, eval_func(R2_eval).
					  derivative_func(R2_deriv).
					  series_func(R2_series));


/*
 *  Expansion coefficients of the two-dimensional R function
 */

static ex R2ex1_eval(const ex &z1, const ex &z2)
{
	return R2ex1(z1, z2).hold();
}

static ex R2ex1_is_really(const ex &z1, const ex &z2)
{
    ex c = sqrt(1 - z1/z2);
    ex logp = log(1 + c);
    ex logm = log(1 - c);
    return c * (logm - logp + I*Pi) - log(-z1) - I*Pi;
}

static ex R2ex1_evalf(const ex &z1, const ex &z2)
{
	return R2ex1_is_really(z1, z2).evalf();
}

static ex R2ex1_deriv(const ex &z1, const ex &z2,
                      unsigned deriv_param)
{
	symbol t1, t2;
	ex e = R2ex1_is_really(t1, t2);
	switch (deriv_param) {
		case 0:
			e = e.diff(t1);
			break;
		case 1:
			e = e.diff(t2);
			break;
	}
	return e.subs(t1 == z1).subs(t2 == z2);
}

static ex R2ex2_eval(const ex &z1, const ex &z2)
{
	return R2ex2(z1, z2).hold();
}

static ex R2ex2_is_really(const ex &z1, const ex &z2)
{
	ex c = sqrt(1 - z1/z2);
	ex logp = log(1 + c);
	ex logm = log(1 - c);
	return   (1+c) * Li2(1-(1-c)/(1+c)) + (1+c) * pow(logm, 2)
			+(1-c) * Li2(1-(1+c)/(1-c)) + (1-c) * pow(logp, 2)
			+ pow(log(z2), 2)/2 + 2*pow(log(z1), 2) - 2*log(z2)*log(z1)
			+ (log(z2) - 2*log(z1)) * ((1+c)*logm + (1-c)*logp)
			- I*Pi/sqrt(-z2) * sqrt(z1-z2) * (2*log(2)+log(z1-z2));
}

static ex R2ex2_evalf(const ex &z1, const ex &z2)
{
	return R2ex2_is_really(z1, z2).evalf();
}

static ex R2ex2_deriv(const ex &z1, const ex &z2,
                      unsigned deriv_param)
{
	symbol t1, t2;
	ex e = R2ex1_is_really(t1, t2);
	switch (deriv_param) {
		case 0:
			e = e.diff(t1);
			break;
		case 1:
			e = e.diff(t2);
			break;
	}
	return e.subs(t1 == z1).subs(t2 == z2);
}

REGISTER_FUNCTION(R2ex1, eval_func(R2ex1_eval).
						 evalf_func(R2ex1_evalf).
						 derivative_func(R2ex1_deriv));
REGISTER_FUNCTION(R2ex2, eval_func(R2ex2_eval).
						 evalf_func(R2ex2_evalf).
						 derivative_func(R2ex2_deriv));


/*
 *  Three-dimensional R function
 *
 *  t:     index
 *  b1..3: parameters
 *  z1..3: arguments
 */

static ex R3_eval(const ex &t,
                  const ex &b1, const ex &b2, const ex &b3,
                  const ex &z1, const ex &z2, const ex &z3,
                  const ex &rho )
{
	// R function is zero when index is zero
	if (t.is_zero())
		return 0;

	// Vanishing parameters recude the R3 function to an R2 function
	if (b1.is_zero())
		return R2(t, b2, b3, z2, z3, rho);
	if (b2.is_zero())
		return R2(t, b1, b3, z1, z3, rho);
	if (b3.is_zero())
		return R2(t, b1, b2, z1, z2, rho);

	// Sum of parameters
	ex beta = b1 + b2 + b3;

	// Feynman parameter formula
	if (t.is_equal(-beta))
		return pow(z1, -b1) * pow(z2, -b2) * pow(z3, -b3);

	// Numeric terms in index and parameters
	ex t_num = NumericTerm(t);
	ex b1_num = NumericTerm(b1);
	ex b2_num = NumericTerm(b2);

	// Increase first parameter
	if (b1_num <= -1 && !beta.is_zero())
      return ((beta + t) * R3(t, b1+1, b2, b3, z1, z2, z3,rho) - t * z1 * R3(t-1, b1+1, b2, b3, z1, z2, z3,rho)) / beta;
  
	// Increase second parameter
	if (b2_num <= -1 && !beta.is_zero())
		return ((beta + t) * R3(t, b1, b2+1, b3, z1, z2, z3,rho) - t * z2 * R3(t-1, b1, b2+1, b3, z1, z2, z3,rho)) / beta;

	return R3(t, b1, b2, b3, z1, z2, z3, rho).hold();
}

static ex R3_deriv(const ex &t,
                   const ex &b1, const ex &b2, const ex &b3,
                   const ex &z1, const ex &z2, const ex &z3,
                   const ex &rho,
                   unsigned deriv_param)
{
	if (deriv_param == 4)		// d/dz1
		return t * b1 / (b1 + b2 + b3) * R3(t-1, b1+1, b2, b3, z1, z2, z3, rho);
	else if (deriv_param == 5)	// d/dz2
		return t * b2 / (b1 + b2 + b3) * R3(t-1, b1, b2+1, b3, z1, z2, z3, rho);
	else if (deriv_param == 6)	// d/dz3
		return t * b3 / (b1 + b2 + b3) * R3(t-1, b1, b2, b3+1, z1, z2, z3, rho);
	else
		throw(std::logic_error("don't know the derivative of this R3 function"));
}

static ex R3_series(const ex &t,
                    const ex &b1, const ex &b2, const ex &b3,
                    const ex &z1, const ex &z2, const ex &z3,
                    const ex &rho,
                    const relational &r,
                    int order,
                    unsigned options)
{
    const symbol &s = ex_to<symbol>(r.lhs());
    
	if (t.is_equal(-2*s) && b1.is_equal(s) && b2.is_equal(s) && b3.is_equal(1) && r.rhs().is_zero()) {
		ex e = 1 - 2 * s * log(z3) + pow(s, 2) * R3ex2(z1, z2, z3) + pow(s, 3) * R3ex3(z1, z2, z3) + Order(pow(s, 4));
		return e.series(r, order);
	} else
		throw(std::logic_error("don't know the series expansion of this particular R3 function"));
}

REGISTER_FUNCTION(R3, eval_func(R3_eval).
					  derivative_func(R3_deriv).
					  series_func(R3_series));


/*
 *  Expansion coefficients of the three-dimensional R function
 */

static ex R3ex2_eval(const ex &z1, const ex &z2, const ex &z3)
{
	return R3ex2(z1, z2, z3).hold();
}

static ex R3ex2_is_really(const ex &z1, const ex &z2, const ex &z3)
{
	ex a = 1 - z1/z3;
	ex b = 1 - z2/z3;
    return 2 * (log(a) * eta(z1, z3) + log(b) * eta(z2, z3) + Li2(a) + Li2(b) + pow(log(z3), 2));
}

static ex R3ex2_evalf(const ex &z1, const ex &z2, const ex &z3)
{
    return R3ex2_is_really(z1, z2, z3).evalf();
}

static ex R3ex3_eval(const ex &z1, const ex &z2, const ex &z3)
{
	return R3ex3(z1, z2, z3).hold();
}

static ex R3ex3_evalf(const ex &z1, const ex &z2, const ex &z3)
{
	//!! this requires the Spence function S12
	return R3ex3(z1, z2, z3).hold();
}

static ex R3ex_deriv(const ex &z1, const ex &z2, const ex &z3,
                     unsigned diff_param)
{
	return 0;
}

REGISTER_FUNCTION(R3ex2, eval_func(R3ex2_eval).
						 evalf_func(R3ex2_evalf).
						 derivative_func(R3ex_deriv));
REGISTER_FUNCTION(R3ex3, eval_func(R3ex3_eval).
						 evalf_func(R3ex3_evalf).
						 derivative_func(R3ex_deriv));


/** Substitute any RXexY() functions by the real expressions. */
ex subs_Rex(const ex &e)
{
	ex Z1 = wild(1), Z2 = wild(2), Z3 = wild(3);
	return e.subs(lst(
		R2ex1(Z1, Z2) == R2ex1_is_really(Z1, Z2),
		R2ex2(Z1, Z2) == R2ex2_is_really(Z1, Z2),
		R3ex2(Z1, Z2, Z3) == R3ex2_is_really(Z1, Z2, Z3)
	));
}


} // namespace xloops
