/** @file oneloop_2pt.cpp
 *
 *  Implementation of functions for one-loop two-point Feynman integrals. */

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
#include "r.h"
#include "utils.h"

namespace xloops {


/*           i-t1  j-t2
 *  Reduce P1    P2     to Scalar2Pt() and Scalar1Pt()
 */

static ex SP2Pt(const ex &q,
                const ex &m1, const ex &m2,
                int n1, int n2,
                const ex &rho, const ex &eps)
{
	ex c;

	if (n1 < 0 && n2 < 0) {

		c = Scalar2Pt(q, m1, m2, -n1, -n2, rho);

	} else {

		symbol l0("l0"), P1("P1"), P2("P2");

		if (n1 < 0) {

			// The case P2^n2 / P1^n1
			ex b = pow(P1 - 2*q*(l0-q) - pow(q, 2) + pow(m1, 2) - pow(m2, 2), n2) * pow(P1, n1);
			b = b.expand();
			for (int i=0; i<=n2; i+=2) {
				for (int j=1; j<=-n1; j++) {
					ex co = b.coeff(l0, i).coeff(P1, -j);
					if (!co.is_zero())
						c += co * OneLoop1Pt(i, m1, j, rho);
				}
				c *= CoeffDim(i, 0, -1, eps);
			}

		} else if (n2 < 0) {

			// The case P1^n1 / P2^n2
			ex b = pow(2*q*l0 + P2 + pow(q, 2) - pow(m1, 2) + pow(m2, 2), n1) * pow(P2, n2);
			b = b.expand();
			for (int i=0; i<=n1; i+=2) {
				for (int j=1; j<=-n2; j++) {
					ex co = b.coeff(l0, i).coeff(P2, -j);
					if (!co.is_zero())
						c += co * OneLoop1Pt(i, m2, j, rho);
				}
				c *= CoeffDim(i, 0, -1, eps);
			}
		}
	}

	return c;
}


/*
 *  Reduce general one-loop two-point function to scalar functions
 */

ex reduce_OneLoop2Pt(int p0, int p1,
                     const ex &q,
                     const ex &m1, const ex &m2,
                     int t1, int t2,
                     const ex &rho, const ex &eps)
{
	ex c;

	if (q.is_zero()) {

		// The special case q=0
		if (pow(m1,2).is_equal(pow(m2,2))) {

			// The special cases where q=0, m1^2-m2^2=0
			c = OneLoop1Pt(p0+p1, m1, t1+t2, rho) * CoeffDim(p0, p1, -1, eps);

		} else {

			// The special case where q=0, m1^2-m2^2!=0
			symbol M1("M1"), M2("M2");
			ex b1 = diff(1/(M1-M2), M2, t2-1);
			ex b2 = diff(1/(M1-M2), M1, t1-1);
			ex b3, b;
			for (int i=0; i<=t1-1; i++)
				b += binomial(t1-1, i) * factorial(i) * OneLoop1Pt(p0+p1, m1, i+1, rho) * diff(b1, M1, t1-1-i);
			for (int j=0; j<=t2-1; j++)
				b3 += binomial(t2-1, j) * factorial(j) * OneLoop1Pt(p0+p1, m2, j+1, rho) * diff(b2, M2, t2-1-j);
			c = subs(b-b3, lst(M1 == pow(m1,2), M2 == pow(m2,2))) / (factorial(t1-1) * factorial(t2-1));
			c *= CoeffDim(p0, p1, -1, eps);
		}

	} else {
		
		// The general case q!=0: decompose into Scalar2Pt() and Scalar1Pt() functions
		symbol P1("P1"), P2("P2"), P3("P3"), P4("P4");
		ex l0 = (P1 - P2 - pow(q, 2) + pow(m1, 2) - pow(m2, 2)) / (2*q);
		ex l1q = pow(l0, 2) - P2 - pow(m2, 2) + I*rho;
		ex b = pow(l0, p0) * pow(l1q, numeric(p1,2));
		b = b.expand();
		for (int i=0; i<=p0+p1; i++)
			for (int j=0; j<=p0+p1; j++){
				c += Pcollect(b, i, j, 0, 0, P1, P2, P3, P4) * SP2Pt(q, m1, m2, i-t1, j-t2, rho, eps);
			}
	}

	// Reduce the remaining OneLoop1Pt() functions
	return tensor_reduction(c, eps);
}


/*
 *  General one-loop two-point function
 */

static ex OneLoop2Pt_eval(const ex &p0, const ex &p1,
                          const ex &q,
                          const ex &m1, const ex &m2,
                          const ex &t1, const ex &t2,
                          const ex &rho)
{
	// The 2-point function vanishes for odd p1 and for odd p0 (when q is zero)
	if (p1.info(info_flags::odd))
		return 0;
	else if (q.is_zero() && p0.info(info_flags::odd))
		return 0;
	else if (p0.is_zero() && p1.is_zero())
		return Scalar2Pt(q, m1, m2, t1, t2, rho);
	else
		return OneLoop2Pt(p0, p1, q, m1, m2, t1, t2, rho).hold();
}

static ex OneLoop2Pt_series(const ex &p0, const ex &p1,
                            const ex &q,
                            const ex &m1, const ex &m2,
                            const ex &t1, const ex &t2,
                            const ex &rho,
                            const relational &r, int order, unsigned options)
{
	// Check parameters for validity
	if (!p0.info(info_flags::nonnegint) || !p1.info(info_flags::nonnegint)
	 || !t1.info(info_flags::posint) || !t2.info(info_flags::posint))
		throw(std::logic_error("OneLoop2Pt_series: p0<0 or p1<0 or t1<=0 or t2<=0"));

	// Check whether we can do the expansion and get parameters
	if (!r.rhs().is_zero())
		throw(std::logic_error("OneLoop2Pt_series: don't know the expansion at eps!=0"));
	const ex eps = r.lhs();

	// Reduce to scalar functions
	ex c = reduce_OneLoop2Pt(
		ex_to<numeric>(p0).to_int(), ex_to<numeric>(p1).to_int(),
		q,
		m1, m2,
		ex_to<numeric>(t1).to_int(), ex_to<numeric>(t2).to_int(),
		rho,
		eps
	);

	// Expand into series in eps
	return c.series(r, order);
}

REGISTER_FUNCTION(OneLoop2Pt, eval_func(OneLoop2Pt_eval).
                              series_func(OneLoop2Pt_series));


/*
 *  Scalar one-loop two-point function (p0 = p1 = 0)
 */

static ex Scalar2Pt_eval(const ex &q,
                         const ex &m1, const ex &m2,
                         const ex &t1, const ex &t2,
                         const ex &rho)
{
	if (q.is_zero() && m1.is_zero() && m2.is_zero())
		return 0;
	else
		return Scalar2Pt(q, m1, m2, t1, t2, rho).hold();
}

static ex Scalar2Pt_series(const ex &q,
                           const ex &m1, const ex &m2,
                           const ex &t1f, const ex &t2f,
                           const ex &rho,
                           const relational &r, int order, unsigned options)
{
	// Check parameters for validity
	if (!t1f.info(info_flags::posint) || !t2f.info(info_flags::posint))
		throw(std::logic_error("Scalar2Pt_series: t1<=0 or t2<=0"));

	// Check whether we can do the expansion and get parameters
	if (!r.rhs().is_zero())
		throw(std::logic_error("Scalar2Pt_series: don't know the expansion at eps!=0"));
	const ex eps = r.lhs();
	int t1 = ex_to<numeric>(t1f).to_int();
	int t2 = ex_to<numeric>(t2f).to_int();

	// Check for special cases
	symbol M1("M1"), M2("M2");
	ex c;
	if (q.is_zero()) {
		if (m1.is_equal(m2)) {

			// q == 0, m1 == m2 is a trivial 1-point function
			c = OneLoop1Pt(0, m1, t1f + t2f, rho);
			
		} else {

			// q == 0, m1 != m2
			ex x1 = (m1.is_zero() ? 0 : M1 * pow(M1-I*rho, -eps));
			ex x2 = (m2.is_zero() ? 0 : M2 * pow(M2-I*rho, -eps));
			ex CC = I * pow(Pi, 2) * tgamma(eps) * pow(Pi, -eps) / (1-eps);
			ex b = (x1 - x2) / (M1 - M2);
			c = b.diff(M1, t1-1).diff(M2, t2-1);
			c = subs(c, lst(M1 == pow(m1,2), M2 == pow(m2,2)));
			c *= CC / (factorial(t1-1) * factorial(t2-1));
		}

	} else {
	
		// General case, reduce to R and Gamma functions
		ex Mdp = (pow(q,2) + M1 - M2).expand();
		ex Mdm = (pow(q,2) - M1 + M2).expand();
		ex CC = I * pow(Pi, 2) * exp(I*Pi*eps) * pow(Pi, -eps) * tgamma(eps) / ((2-4*eps) * pow(q,2));
		ex b = Mdp * R2(-eps, eps-numeric(1,2), 1, -M1+I*rho, -pow(Mdp,2) / (4*pow(q,2)), rho)
		     + Mdm * R2(-eps, eps-numeric(1,2), 1, -M2+I*rho, -pow(Mdm,2) / (4*pow(q,2)), rho);
	 	c = b.diff(M1, t1-1).diff(M2, t2-1);
	 	c = subs(c, lst(M1 == pow(m1,2), M2 == pow(m2,2)));
	 	c *= CC / (factorial(t1-1) * factorial(t2-1));
	}

	// Expand into series in eps
	return c.series(r, order);
}

REGISTER_FUNCTION(Scalar2Pt, eval_func(Scalar2Pt_eval).
                             series_func(Scalar2Pt_series));

} // namespace xloops
