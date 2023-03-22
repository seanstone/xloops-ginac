/** @file oneloop_3pt.cpp
 *
 *  Implementation of functions for one-loop three-point Feynman integrals. */

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


/*
 *  This function S322Pt() is used to shift the one-loop two-point functions
 *  with 2 parallel dimensions to OneLoop2Pt() with 1 parallel dimension
 *  (using a Lorentz shift).
 *
 *  P1 = (l0+q10)^2 - l1^2 - l+^2 - m1^2
 *  P2 = (l0+q20)^2 - (l1+q21)^2 - l+^2 - m2^2
 */

ex reduce_S322Pt(int p0, int p1, int p2,
                 const ex &q10, const ex &q20, const ex &q21,
                 const ex &m1, const ex &m2,
                 int t1, int t2,
                 const ex &rho, const ex &eps)
{
	// Check whether the function can be reduced to Scalar2Pt() functions
	ex Q = q20 - q10;
	ex c;
	if (q21.is_zero()) {

		// The case q21 = 0 is involved because the general three-point case
		// also calls it with q21 = 0.
		symbol l0("l0");
		ex d = pow(l0 - q10, p0).expand();
		for (int n=0; n<=p0; n++) {
			ex co = d.coeff(l0, n);
			if (!co.is_zero())
				c += co * OneLoop2Pt(n, p1+p2, Q, m2, m1, t1, t2, rho);
		}
		c *= CoeffDim(-1, p1, p2, eps);

	} else if (!Q.is_zero() && !(pow(Q, 2) - pow(q21, 2)).expand().is_zero()) {

		// Use Lorentz boost, It doesn't matter whether Q^2-q21^2 is greater or smaller than zero.
		symbol l0("l0"), l1("l1");
		ex P = sqrt(power(Q, 2) - power(q21, 2));
		ex lgamma = Q / P, lbeta = q21 / Q;
		ex d = (power(lgamma * (l0+lbeta*l1) - q10, p0) * power(lgamma * (lbeta*l0+l1), p1)).expand();
		for (int n=0; n<=p0+p1; n++) {
			for (int m=0; m<=p0+p1; m+=2) {
				ex co = d.coeff(l0, n).coeff(l1, m);
				if (!co.is_zero())
					c += co * OneLoop2Pt(n, m+p2, P, m2, m1, t2, t1, rho) * CoeffDim(-1, m, p2, eps);
			} 
		}

	} else {

		// Cannot reduce
		return S322Pt(p0, p1, p2, q10, q20, q21, m1, m2, t1, t2, rho);
	}

	// Reduce the remaining OneLoop2Pt() functions
	return tensor_reduction(c, eps);
}

static ex S322Pt_eval(const ex &p0, const ex &p1, const ex &p2,
                      const ex &q10, const ex &q20, const ex &q21,
                      const ex &m1, const ex &m2,
                      const ex &t1, const ex &t2,
                      const ex &rho)
{
	// This function vanishes for odd p2
	if (p2.info(info_flags::odd))
		return 0;
	else
		return S322Pt(p0, p1, p2, q10, q20, q21, m1, m2, t1, t2, rho).hold();
}

static ex S322Pt_series(const ex &p0f, const ex &p1f, const ex &p2f,
                        const ex &q10, const ex &q20, const ex &q21,
                        const ex &m1, const ex &m2,
                        const ex &t1f, const ex &t2f,
                        const ex &rho,
                        const relational &r, int order, unsigned options)
{
	// Check parameters for validity
	if (!p0f.info(info_flags::nonnegint) || !p1f.info(info_flags::nonnegint) || !p2f.info(info_flags::nonnegint)
	 || !t1f.info(info_flags::posint) || !t2f.info(info_flags::posint))
		throw(std::logic_error("S322Pt_series: p0<0 or p1<0 or p2<0 or t1<=0 or t2<=0"));

	// Check whether we can do the expansion and get parameters
	if (!r.rhs().is_zero())
		throw(std::logic_error("S322Pt_series: don't know the expansion at eps!=0"));
	const ex eps = r.lhs();

	// Reduce to scalar functions
	ex red = reduce_S322Pt(
		ex_to<numeric>(p0f).to_int(), ex_to<numeric>(p1f).to_int(), ex_to<numeric>(p2f).to_int(),
		q10, q20, q21,
		m1, m2,
		ex_to<numeric>(t1f).to_int(), ex_to<numeric>(t2f).to_int(),
		rho,
		eps);

	// If the reduction worked, we need can series expand the result
	if (!is_ex_the_function(red, S322Pt))
		return red.series(eps, order);

	// Reduction not possible:
	// (Q==0 and pow(Q,2)-pow(q21,2)!=0) or (Q!=0 and pow(Q,2)-pow(q21,2)==0)
	// Solve by residue theorem for Q^2-q21^2==0. Don't boost.

	// OK, as it is, I took the results below from Lars Bruecher,
	// but I still doubt this result. I will check it again some day.
	// Don't blame me about this part if there is any problem.

	symbol M1("M1"), M2("M2");
	int p0 = ex_to<numeric>(p0f).to_int();
	int p1 = ex_to<numeric>(p1f).to_int();
	int p2 = ex_to<numeric>(p2f).to_int();
	int t1 = ex_to<numeric>(t1f).to_int();
	int t2 = ex_to<numeric>(t2f).to_int();

	ex Q = q20 -q10;
	ex s = (Q.is_equal(q21) ? -1 : 1);

	ex r0[2] = {q10, 0};
	ex r1[2] = {q20, q21};
	ex a[2]  = {2*Q, -2*Q};
	ex b[2]  = {-2*q21, 2*q21};
	ex c1[2] = {pow(Q,2)-pow(q21,2)-M2+M1, pow(Q,2)-pow(q21,2)-M1+M2};
	ex MM[2] = {M1, M2};
	
	ex c;
	if (Q.is_zero()) { // Q^2-q21^2 != 0

		ex y0[2], yp[2],ym[2];
		for (int k=0; k<=1; k++) {
			y0[k] = c1[k]/b[k]*s + r0[k];
			yp[k] = c1[k]/b[k]*s + sqrt(pow(c1[k], 2) / pow(b[k], 2) + MM[k] - I*rho);
			ym[k] = c1[k]/b[k]*s - sqrt(pow(c1[k], 2) / pow(b[k], 2) + MM[k] - I*rho);
		}
		for (int k=0; k<=1; k++) {
			ex res1;
			for (int i=0; i<=p0; i++) {
				res1 += binomial(p0, i) * pow(-y0[k], p0-i) * beta(2*eps-p2-i-1, i+1)
				      * (R2(1+i+p2-2*eps, eps-numeric(p2,2), eps-numeric(p2,2), -yp[k], -ym[k], rho)
				       - pow(-1, i+p2-2*eps) * R2(1+i+p2-2*eps, eps-numeric(p2,2), eps-numeric(p2,2), yp[k], ym[k], rho));
			}
			c += res1 * pow(-c1[k]/b[k] - r1[k], p1);
		}

	} else { // Q^2-q21^2 == 0

		for (int k=0; k<=1; k++) {
			ex res1;
			for (int i=0; i<=p0; i++) {
				res1 += binomial(p0, i) * pow(-(c1[k]/a[k]+2*r0[k]), p0-i) * beta(-p1-numeric(p2,2)+eps-i-1, i+1)
				      * (R2(p1+numeric(p2,2)-eps+i+1, -p1, -numeric(p2,2)+eps, c1[k]+2*r1[k]*s, 2*a[k]*(MM[k]-I*rho)/c1[k], rho)
				       - pow(-1, i+p1+numeric(p2,2)-eps) * R2(p1+numeric(p2,2)-eps+i+1, -p1, -numeric(p2,2)+eps, -(c1[k]+2*r1[k]*s), -2*a[k]*(MM[k]-I*rho)/c1[k], rho));
			}
			c += res1 * pow(-c1[k]/2/a[k], numeric(p2,2)-eps);
		}

		c *= pow(numeric(1,2), p0) * pow(-s/2, p1);
	}


	c *= -I * pow(Pi, 2-eps) * tgamma(eps-numeric(p2,2)) * tgamma(1-eps+numeric(p2,2)) * pow(-1, numeric(p2,2)-eps) / tgamma(1-eps);
	
	c = c.diff(M1, t1-1).diff(M2, t2-1);
	c *= 1 / (factorial(t1-1) * factorial(t2-1));
	c = c.subs(lst(M1 == pow(m1, 2), M2 == pow(m2, 2)));

	// Expand into series in eps
	return c.series(eps, order);
}

REGISTER_FUNCTION(S322Pt, eval_func(S322Pt_eval).
						  series_func(S322Pt_series));


/*
 *  SP3Pt() is used to reduce the tensor integrals of the form
 *
 *                     ni          ni    nj
 *                    P           P     P 
 *                     i           i     j
 *                ---------- and -----------
 *                 P     P           P
 *                  j     k           k
 *
 *  to OneLoop1Pt() and two-point functions with 2 parallel dimensions.
 *  This operation is valid for any values of q10, q20, q21 and any masses.
 */

static ex SP3Pt(const ex &q10, const ex &q20, const ex &q21,
                const ex &m1, const ex &m2, const ex &m3,
                int n1, int n2, int n3,
                const ex &rho, const ex &eps)
{
	symbol l0("l0"), l1("l1"), P1("P1"), P2("P2"), P3("P3");
	ex c00 = pow(q10, 2) - pow(m1, 2) + pow(m3, 2);
	ex c10 = 2*q20;
	ex c11 = pow(m3, 2) - pow(m2, 2) + pow(q20, 2) - pow(q21, 2);
	ex c20 = pow(m3, 2) - I*rho;
	ex c30 = pow(q10,2) - pow(q20, 2) + pow(q21, 2) + pow(m2, 2) - pow(m1, 2);

	int n11 = (n1 == 0) ? 1 : n1,
	    n22 = (n2 == 0) ? 1 : n2,
	    n33 = (n3 == 0) ? 1 : n3;
	
	ex c;
	
	if (n11*n22*n33 < 0) {
		
		if (n1 < 0 && n2 < 0 && n3 < 0) {
			
			c = Scalar3Pt(q10, q20, q21, m1, m2, m3, -n1, -n2, -n3, rho);
			
		} else {

			if (n1 < 0) {

				// The case (P2^n2 * P3^n3) / P1^(-n1)
				ex e = P1 - c00 - 2*q10*(l0-q10);
				ex d = e + c10*(l0-q10) - 2*q21*l1 + c11;                               
				ex b = pow(d, n2) * pow(e, n3) * pow(P1, n1);
				b = b.expand();
				for (int i=0; i<=n2+n3; i+=2) {
					for (int j=0; j<=n2+n3; j+=2) {
						for (int k=1; k<=-n1; k++) {
							ex co = b.coeff(l0, i).coeff(l1, j).coeff(P1, -k);
							if (!co.is_zero())
								c += co * OneLoop1Pt(i+j, m1, k, rho);
						}
						c *= CoeffDim(i, j, -1, eps);
					}
				}

			} else if (n2 < 0) {

				// The case (P1^n1 * P3^n3) / P2^(-n2)
				ex e = 2*q21*(l1-q21) - c10*(l0-q20) + P2 - c11;
				ex d = 2*q10*(l0-q20) + e + c00;
				ex b = pow(d, n1) * pow(e, n3) * pow(P2, n2);
				b = b.expand();
				for (int i=0; i<=n1+n3; i+=2) {
					for (int j=0; j<=n1+n3; j+=2) {
						for (int k=1; k<=-n2; k++) {
							ex co = b.coeff(l0, i).coeff(l1, j).coeff(P2, -k);
							if (!co.is_zero())
								c += co * OneLoop1Pt(i+j, m2, k, rho);
						}
						c *= CoeffDim(i, j, -1, eps);
					}
				}

			} else {

				// The case (P1^n1 * P2^n2) / P3^(-n3)
				ex e = c10*l0 + P3 + c11 - 2*q21*l1;
				ex d = 2*q10*l0 + P3 + c00;
				ex b = pow(d, n1) * pow(e, n2) * pow(P3, n3);
				b = b.expand();
				for (int i=0; i<=n1+n2; i+=2) {
					for (int j=0; j<=n1+n2; j+=2) {
						for (int k=1; k<=-n3; k++) {
							ex co = b.coeff(l0, i).coeff(l1, j).coeff(P3, -k);
							if (!co.is_zero())
								c += co * OneLoop1Pt(i+j, m3, k, rho);
						}
						c *= CoeffDim(i, j, -1, eps);
					}
				}
			}
		}

	} else if (n1 < 0 || n2 < 0 || n3 < 0) {

		ex c1, c2;

		if (n1 >= 0) {
		  
		  // The case P1^n1 / (P2^(-n2) * P3^(-n3))
		  ex d = 2*q10*l0 + P3 + c00;
		  ex b = pow(d, n1) * pow(P2, n2) * pow(P3, n3);
		  b = b.expand();
		  for (int i=0; i<=n1; i++) {
			for (int j= 0; j <= n1+n3; j++) {
			  ex e = Pcollect(b, i, j, n2, 0, l0, P3, P2, l1) * pow(l0-q20, i) * pow(2*q21*(l1-q21)+P2-c10*(l0-q20)-c11, j) * pow(P2, n2);
			  e = e.expand();
			  for (int l=0; l<=i+j; l+=2) {
				for (int k=0; k<=i+j; k+=2) {
				  for (int h=1; h<=-n2; h++) {
					ex co = Pcollect(e, l, k, -h, 0, l0, l1, P2, P3);
					if (!co.is_zero())
					  c1 += co * OneLoop1Pt(l+k, m2, h, rho);
				  }
				  c1 *= CoeffDim(l, k, -1, eps);
				}
			  }
			}
			for (int s=1; s<=-n3; s++){
			  ex co = Pcollect(b, i, n2, -s, 0, l0, P2, P3, l1);
			  if (!co.is_zero())
				c2 += co * S322Pt(i, 0, 0, 0, q20, q21, m3, m2, s, -n2, rho);
			}
		  }
		  
		} else if (n2 >= 0) {
		  
		  // The case P2^n2 / (P1^(-n1) * P3^(-n3))
		  ex d = P3 + c10*l0 - 2*q21*l1 + c11;
		  ex b = pow(d, n2) * pow(P1, n1) * pow(P3, n3);
		  b = b.expand();
		  for (int i=0; i<=n2; i++) {
			for (int j=0; j<=n2; j++) {
			  for (int k=0; k<=n2+n3; k++) {
				ex e = Pcollect(b, i, j, k, n1, l0, l1, P3, P1) * pow(l0-q10, i) * pow(l1, j) * pow(P1-c00-2*q10*(l0-q10), k) * pow(P1, n1);
				e = e.expand();
				for (int l=0; l<=i+j+k; l+=2) {
				  for (int o=0; o<=i+j+k; o+=2) {
					for (int h=1; h<=-n1; h++) {
					  ex co = Pcollect(e, l, o, 0, -h, l0, l1, P3, P1);
					  if (!co.is_zero())
						c1 += co * OneLoop1Pt(l+o, m1, h, rho);
					}
					c1 *= CoeffDim(l, o, -1, eps);
				  }
				}
			  }
			  for (int s=1; s<=-n3; s++) {
				ex co = Pcollect(b, i, j, n1, -s, l0, l1, P1, P3);
				if (!co.is_zero())
				  c2 += co * OneLoop2Pt(i, j, q10, m1, m3, -n1, s, rho);
			  }
			  c2 *= CoeffDim(-1, j, 0, eps);
			}
		  }
		  
		} else {
		  
		  // The case P3^n3 / (P1^(-n1) * P2^(-n2))
		  ex d = P1 - c00 - 2*q10*l0;
		  ex b = pow(d, n3) * pow(P1, n1) * pow(P2, n2);
		  b = b.expand();
		  for (int i=0; i<=n3; i++) {
			for (int j=0; j<=n3+n1; j++) {
			  ex e = Pcollect(b, i, 0, j, n2, l0, l1, P1, P2) * pow(l0-q20, i) * pow(2*(q10-q20)*(l0-q20)+2*q21*(l1-q21)+P2+c30, j) * pow(P2, n2);
			  e = e.expand();
			  for (int k=0; k<=i+j; k+=2) {
				for (int l=0; l<=i+j; l+=2) {
				  for (int h=1; h<=-n2; h++) {
					ex co = Pcollect(e, k, l, 0, -h, l0, l1, P1, P2);
					if (!co.is_zero())
					  c1 += co * OneLoop1Pt(k+l, m2, h, rho);
				  }
				  c1 *= CoeffDim(k, l, -1, eps);
				}
			  }
			}
			for (int s=1; s<=-n1; s++) {
			  ex co = Pcollect(b, i, 0, -s, n2, l0, l1, P1, P2);
			  if (!co.is_zero())
				c2 += co * S322Pt(i, 0, 0, q10, q20, q21, m1, m2, s, -n2, rho);
			}
		  }
		}

		c = c1 + c2;

	} else
		return 0;

	return c;
}


/*
 *  Reduce general one-loop three-point function to scalar functions
 */

ex reduce_OneLoop3Pt(int p0, int p1, int p2,
                     const ex &q10, const ex &q20, const ex &q21,
                     const ex &m1, const ex &m2, const ex &m3,
                     int t1, int t2, int t3,
                     const ex &rho, const ex &eps)
{
	symbol P1("P1"), P2("P2"), P3("P3"), P4("P4"), M1("M1"), M3("M3");
	ex c;

	if (!(q10*q21).is_zero()) {
		// The case q10*q21 != 0
		ex c00 = pow(q10, 2) - pow(m1, 2) + pow(m3, 2);
		ex c10 = 2*q20;
		ex c11 = pow(m3,2) - pow(m2, 2) + pow(q20, 2) - pow(q21, 2);
		ex c20 = pow(m3,2) - I*rho;

		ex l0 = (P1-P3-c00) / (2*q10);
		ex l1 = (c10*l0+P3-P2+c11) / (2*q21);
		ex l2 = pow(l0, 2) - pow(l1, 2) - P3 - c20;

		ex A = pow(l0, p0) * pow(l1, p1) * pow(l2, numeric(p2,2));
		A = A.expand();

		for (int i=0; i<=p0+p1+p2; i++) {
			for (int j=0; j<=p0+p1+p2; j++) {
				for ( int k=0; k<=p0+p1+p2; k++) {
					ex co = Pcollect(A, i, j, k, 0, P1, P2, P3, P4);
					if (!co.is_zero())
						c += co * SP3Pt(q10, q20, q21, m1, m2, m3, i-t1, j-t2, k-t3, rho, eps);
				}
			}
		}

	} else {

		// The case q10*q21 = 0
		if (q10.is_zero()) {

			// q10 = 0, q21 != 0
			if (!(pow(m1, 2) - pow(m3, 2)).is_zero()) {
				ex A1 = diff(1/(M1-M3), M3, t3-1);
				ex A3 = diff(1/(M1-M3), M1, t1-1);
				ex d, e;
				for (int i=0; i<=t1-1; i++)
					d += binomial(t1-1, i) * factorial(i) * S322Pt(p0, p1, p2, 0, q20, q21, m1, m2, i+1, t2, rho) * diff(A1, M1, t1-1-i);
				for (int i=0; i<=t3-1; i++)
					e += binomial(t3-1, i) * factorial(i) * S322Pt(p0, p1, p2, 0, q20, q21, m3, m2, i+1, t2, rho) * diff(A3, M3, t3-1-i);
				c = subs(d-e, lst(M1 == pow(m1, 2), M3 == pow(m3,2))) / (factorial(t1-1) * factorial(t3-1));
			} else { // P1=P3
				c = S322Pt(p0, p1, p2, 0, q20, q21, m1, m2, t1+t3, t2, rho);
			}

		} else {

			// q10 != 0, q21 = 0
			if  ((p1 & 1) == 0) {
				
				ex l0 = (P1 - P3 + pow(m1, 2) - pow(m3, 2) - pow(q10, 2)) / (2*q10);
				ex l1 = pow(l0, 2) - P3 - pow(m3, 2) + I*rho;
				ex A = pow(l0, p0) * pow(l1, numeric(p1+p2,2));
				A = A.expand();
				for (int i=0; i<=p0+p2+p1; i++) {
					for (int j=0; j<=p0+p2+p1; j++) {
						ex co = Pcollect(A, i, 0, j, 0, P1, P2, P3, P4);
						if (!co.is_zero())
							c += co * SP3Pt(q10, q20, 0, m1, m2, m3, i-t1, -t2, j-t3, rho, eps);
					}
				}
				c *= CoeffDim(-1, p1, p2, eps);
			}
		}
	}
	// Reduce the remaining S322Pt() functions
	return tensor_reduction(c, eps);
}


/*
 *  General one-loop three-point function
 */

static ex OneLoop3Pt_eval(const ex &p0, const ex &p1, const ex &p2,
                          const ex &q10, const ex &q20, const ex &q21,
                          const ex &m1, const ex &m2, const ex &m3,
                          const ex &t1, const ex &t2, const ex &t3,
                          const ex &rho)
{
	// The 3-point function vanishes for odd p2
	if (p2.info(info_flags::odd))
		return 0;
	else if (p0.is_zero() && p1.is_zero() && p2.is_zero())
		return Scalar3Pt(q10, q20, q21, m1, m2, m3, t1, t2, t3, rho);
	else
		return OneLoop3Pt(p0, p1, p2, q10, q20, q21, m1, m2, m3, t1, t2, t3, rho).hold();
}

static ex OneLoop3Pt_series(const ex &p0, const ex &p1, const ex &p2,
                            const ex &q10, const ex &q20, const ex &q21,
                            const ex &m1, const ex &m2, const ex &m3,
                            const ex &t1, const ex &t2, const ex &t3,
                            const ex &rho,
                            const relational &r, int order, unsigned options)
{
	// Check parameters for validity
	if (!p0.info(info_flags::nonnegint) || !p1.info(info_flags::nonnegint) || !p2.info(info_flags::nonnegint)
	 || !t1.info(info_flags::posint) || !t2.info(info_flags::posint) || !t3.info(info_flags::posint))
		throw(std::logic_error("OneLoop3Pt_series: p0<0 or p1<0 or p2<0 or t1<=0 or t2<=0 or t3<=0"));

	// Check whether we can do the expansion and get parameters
	if (!r.rhs().is_zero())
		throw(std::logic_error("OneLoop3Pt_series: don't know the expansion at eps!=0"));
	const ex eps = r.lhs();

	// Reduce to scalar functions
	ex c = reduce_OneLoop3Pt(
		ex_to<numeric>(p0).to_int(), ex_to<numeric>(p1).to_int(), ex_to<numeric>(p2).to_int(),
		q10, q20, q21,
		m1, m2, m3,
		ex_to<numeric>(t1).to_int(), ex_to<numeric>(t2).to_int(), ex_to<numeric>(t3).to_int(),
		rho,
		eps
	);

	// Expand into series in eps
	return c.series(r, order);
}

REGISTER_FUNCTION(OneLoop3Pt, eval_func(OneLoop3Pt_eval).
                              series_func(OneLoop3Pt_series));


/*
 *  Scalar one-loop three-point function (p0 = p1 = p2 = 0)
 */

static ex Scalar3Pt_eval(const ex &q10, const ex &q20, const ex &q21,
                         const ex &m1, const ex &m2, const ex &m3,
                         const ex &t1, const ex &t2, const ex &t3,
                         const ex &rho)
{
	return Scalar3Pt(q10, q20, q21, m1, m2, m3, t1, t2, t3, rho).hold();
}

static ex Scalar3Pt_series(const ex &q10, const ex &q20, const ex &q21,
                           const ex &m1, const ex &m2, const ex &m3,
                           const ex &t1f, const ex &t2f, const ex &t3f,
                           const ex &rho,
                           const relational &r, int order, unsigned options)
{
	// Check parameters for validity
	if (!t1f.info(info_flags::posint) || !t2f.info(info_flags::posint) || !t3f.info(info_flags::posint))
		throw(std::logic_error("Scalar3Pt_series: t1<=0 or t2<=0 or t3<=0"));

	// Check whether we can do the expansion and get parameters
	if (!r.rhs().is_zero())
		throw(std::logic_error("Scalar2Pt_series: don't know the expansion at eps!=0"));
	const ex eps = r.lhs();
	int t1 = ex_to<numeric>(t1f).to_int();
	int t2 = ex_to<numeric>(t2f).to_int();
	int t3 = ex_to<numeric>(t3f).to_int();

	// Squares of masses and momenta
	symbol M1("M1"), M2("M2"), M3("M3");
	ex m_2[3] = {M1, M2, M3};
	ex q10_2 = pow(q10, 2);
	ex q20_2 = pow(q20, 2);
	ex q21_2 = pow(q21, 2);

	// Abbreviations (here, we change the order of the index a[k,l]-> a_l_[k])
	ex a1[3] = {2*(q20-q10), -2*(q20-q10), 2*q10};
	ex a2[3] = {-2*q10, -2*q20, 2*q20};
	ex b1[3] = {-2*q21, 2*q21, 0};
	ex b2[3] = {0, 2*q21, -2*q21};
	ex c1[3] = {pow(q20-q10,2)-q21_2-m_2[1]+m_2[0], pow(q20-q10,2)-q21_2-m_2[0]+m_2[1], q10_2-m_2[0]+m_2[2]};
	ex c2[3] = {q10_2-m_2[2]+m_2[0], q20_2-q21_2-m_2[2]+m_2[1], q20_2-q21_2-m_2[1]+m_2[2]};

	// First, check infrared singularity
	ex IR[3] = {1, 1, 1};
	if (m3.is_zero() && expand(q10_2 - pow(m1, 2)).is_zero()
		&& expand(q20_2 - q21_2 - pow(m2, 2)).is_zero()) {
		IR[2] = 0;
	} else if (m1.is_zero() && expand(q10_2 - pow(m3, 2)).is_zero()
			   && expand(pow(q20-q10, 2) - q21_2 - pow(m2, 2)).is_zero()) {
		IR[0] = 0;
	} else if (m2.is_zero() && expand(q20_2 - q21_2 - pow(m3, 2)).is_zero()
			   && expand(pow(q20-q10,2)-q21_2-pow(m1,2)).is_zero()) {
		IR[1] = 0;
	} // otherwise IR convergent
	/*
	 *  There are 3 groups of the special cases:
	 *   1) q10 = 0, with any dynamical values, then the three-point functions reduce
	 *      to two-point functions, trivial.
	 *   2) q21 = 0 and no other constraints, we treat it separately. If additionally
	 *      q10 = q20 or q20 = 0 then it will be reduced to two-point functions.
	 *   3) Collinear divergencies (special cases 11,12, 13 [Lars]), for me the known
	 *      results are not clear then I will not do anything till I have free time.
	 */

	ex c;
	if (evalf(q10).is_zero()) {

		// Reduce to two-point functions
		if (m1.is_equal(m3))
			c = S322Pt(0, 0, 0, 0, q20, q21, m1, m2, t3+t1, t2, rho);
		else {
			ex A1 = diff(1/(M1-M3), M3, t3-1);
			ex A3 = diff(1/(M1-M3), M1, t1-1);
			ex d, e;
			for (int i=0; i<=t1-1; i++)
				d += binomial(t1-1, i) * factorial(i) * S322Pt(0,0,0, 0, q20, q21, m1, m2, i+1, t2, rho) * diff(A1, M1, t1-1-i);
			for (int i=0; i<=t3-1; i++)
				e += binomial(t3-1, i) * factorial(i) * S322Pt(0,0,0, 0, q20, q21, m3, m2, i+1, t2, rho) * diff(A3, M3, t3-1-i);
			c = subs(d-e, lst(M1 == pow(m1, 2), M3 == pow(m3,2))) / (factorial(t1-1) * factorial(t3-1));
		}

	} else if (!evalf(q21).is_zero()) {
		// Should consider to the special cases 11 and 12.
		ex res;
		for (int k=0; k<3 ; k++) {
			if ( !IR[k].is_zero()) {
				ex y31 = ( -c1[k]*(a2[k] + b2[k]) + c2[k]*(a1[k] + b1[k])) / (a2[k]*b1[k]-a1[k]*b2[k]);
				ex y32 = y31;
				ex root1 = sqrt(pow(c1[k], 2) + (pow(b1[k], 2) - pow(a1[k], 2)) * (m_2[k] - I*rho));
				ex root2 = sqrt(pow(c2[k], 2) + (pow(b2[k], 2) - pow(a2[k], 2)) * (m_2[k] - I*rho));
				ex lsum;
				
				if ( (q20-q21).normal() == 0) {
					// The special  case 11, the terms (k,l) = {(1,2), (2,2)} are vanished.
					ex S1 = (a1[k] - b1[k]) / (a1[k] + b1[k]);
					ex y11 = (c1[k] + root1) / (a1[k] - b1[k]);
					ex y21 = (c1[k] - root1) / (a1[k] - b1[k]);
					ex y12 = (c2[k] + root2) / (a2[k] - b2[k]);
					ex y22 = (c2[k] - root2) / (a2[k] - b2[k]);
					
					if (k==0){
						ex S2 = (a2[k] - b2[k]) / (a2[k] + b2[k]);
						lsum = - pow(S1, -eps)
							* (R3(-2 * eps, eps, eps, 1, y11, y21, y31, rho) + R3(-2 * eps, eps, eps, 1, -y11, -y21, -y31, rho))
							+ pow(S2, -eps)
							* (R3(-2 * eps, eps, eps, 1, y12, y22, y32, rho) + R3(-2 * eps, eps, eps, 1, -y12, -y22, -y32, rho));	
					} else {
						lsum = - pow(S1, -eps)
							* (R3(-2 * eps, eps, eps, 1, y11, y21, y31, rho) + R3(-2 * eps, eps, eps, 1, -y11, -y21, -y31, rho));	
						
					}
					res += lsum / (a1[k]*b2[k] - a2[k]*b1[k]);

					//====================================================			
					
				} else if( (q20-q21-q10).normal() == 0 ) {
					// The special  case 11, the terms (k,l) = {(0,1), (1,1)} are vanished.
					ex S2 = (a2[k] - b2[k]) / (a2[k] + b2[k]);
					ex y11 = (c1[k] + root1) / (a1[k] - b1[k]);
					ex y21 = (c1[k] - root1) / (a1[k] - b1[k]);
					ex y12 = (c2[k] + root2) / (a2[k] - b2[k]);
					ex y22 = (c2[k] - root2) / (a2[k] - b2[k]);
					if (k==2){
						ex S1 = (a1[k] - b1[k]) / (a1[k] + b1[k]);	
						lsum = - pow(S1, -eps)
							* (R3(-2 * eps, eps, eps, 1, y11, y21, y31, rho) + R3(-2 * eps, eps, eps, 1, -y11, -y21, -y31, rho))
							+ pow(S2, -eps)
							* (R3(-2 * eps, eps, eps, 1, y12, y22, y32, rho) + R3(-2 * eps, eps, eps, 1, -y12, -y22, -y32, rho));
					} else{
						lsum = pow(S2, -eps)
							* (R3(-2 * eps, eps, eps, 1, y12, y22, y32, rho) + R3(-2 * eps, eps, eps, 1, -y12, -y22, -y32, rho));
				
					}
					res += lsum / (a1[k]*b2[k] - a2[k]*b1[k]);
					//======================================================
					
				} else if( (q20+q21).normal() == 0) {

					// The special  case 12-a, the terms (k,l) = {(1,2), (2,2)} are reduced to R2().
					
					ex y42 = (m_2[k]-I*rho)*(a2[k]+b2[k])/(2*c2[k]);
					ex S1 = (a1[k] - b1[k]) / (a1[k] + b1[k]);
					ex S3 = 2*c2[k]/(a2[k]+b2[k]);
					if (k==0){
						ex S2 = (a2[k] - b2[k]) / (a2[k] + b2[k]);
						ex y11 = (c1[k] + root1) / (a1[k] - b1[k]);
						ex y21 = (c1[k] - root1) / (a1[k] - b1[k]);
						ex y12 = (c2[k] + root2) / (a2[k] - b2[k]);
						ex y22 = (c2[k] - root2) / (a2[k] - b2[k]);
						lsum = - pow(S1, -eps)
							* (R3(-2 * eps, eps, eps, 1, y11, y21, y31, rho) + R3(-2 * eps, eps, eps, 1, -y11, -y21, -y31, rho))
							+ pow(S2, -eps)
							* (R3(-2 * eps, eps, eps, 1, y12, y22, y32, rho) + R3(-2 * eps, eps, eps, 1, -y12, -y22, -y32, rho));
					}
					else{
						ex y11 = (c1[k] + root1) / (a1[k] - b1[k]);
						ex y21 = (c1[k] - root1) / (a1[k] - b1[k]);
						lsum = - pow(S1, -eps)
							* (R3(-2 * eps, eps, eps, 1, y11, y21, y31, rho) + R3(-2 * eps, eps, eps, 1, -y11, -y21, -y31, rho))
							+pow(S3,-eps)*2
							*(R2(-eps, eps, 1, y42, y32, rho) + pow(-1,eps) * R2(-eps, eps, 1, -y42, -y32, rho));
						// Note that we have to multiply with 2 in the second term due to
						// the factorization of beta(2*eps,1) out.
					}
					res += lsum / (a1[k]*b2[k] - a2[k]*b1[k]);

					//===============================================
					
				} else if ((q20+q21-q10).normal() == 0) {
					// The special  case 12, the terms (k,l) = {(0,1), (1,1)} are vanished.
					ex y41 = (m_2[k]-I*rho)*(a1[k]+b1[k])/(2*c1[k]);
					ex S2 = (a2[k] - b2[k]) / (a2[k] + b2[k]);
					ex S3 = 2*c1[k]/(a1[k]+b1[k]);
					if (k==2){
						ex S1 = (a1[k] - b1[k]) / (a1[k] + b1[k]);
						ex y11 = (c1[k] + root1) / (a1[k] - b1[k]);
						ex y21 = (c1[k] - root1) / (a1[k] - b1[k]);
						ex y12 = (c2[k] + root2) / (a2[k] - b2[k]);
						ex y22 = (c2[k] - root2) / (a2[k] - b2[k]);
						lsum = - pow(S1, -eps)
							* (R3(-2 * eps, eps, eps, 1, y11, y21, y31, rho) + R3(-2 * eps, eps, eps, 1, -y11, -y21, -y31, rho))
							+ pow(S2, -eps)
							* (R3(-2 * eps, eps, eps, 1, y12, y22, y32, rho) + R3(-2 * eps, eps, eps, 1, -y12, -y22, -y32, rho));
					} else{
						ex y12 = (c2[k] + root2) / (a2[k] - b2[k]);
						ex y22 = (c2[k] - root2) / (a2[k] - b2[k]);
						lsum = -2*pow(S3, -eps)*
							(R2(-eps, eps, 1, y41, y31, rho) + pow(-1, eps) * R2(-eps, eps, 1, -y41, -y31, rho))
							+pow(S2, -eps)
							* (R3(-2 * eps, eps, eps, 1, y12, y22, y32, rho) + R3(-2 * eps, eps, eps, 1, -y12, -y22, -y32, rho));
					}
					
					res += lsum / (a1[k]*b2[k] - a2[k]*b1[k]);
					
				} else if ( y31 == 0) {
					// Do not know how to do yet
					if (c1[k].is_zero() && c2[k].is_zero()){
						res += 0;
					} else {
						throw(std::logic_error("Sorry! This special case is not implemented yet"));
					}
				} else {
					// General case
					// y[k,l](i)-> y_i_l
					ex y11 = (c1[k] + root1) / (a1[k] - b1[k]);
					ex y21 = (c1[k] - root1) / (a1[k] - b1[k]);
					ex y12 = (c2[k] + root2) / (a2[k] - b2[k]);
					ex y22 = (c2[k] - root2) / (a2[k] - b2[k]);
					ex S1 = (a1[k] - b1[k]) / (a1[k] + b1[k]);
					ex S2 = (a2[k] - b2[k]) / (a2[k] + b2[k]);
					ex lsum = - pow(S1, -eps)
						* (R3(-2 * eps, eps, eps, 1, y11, y21, y31, rho) + R3(-2 * eps, eps, eps, 1, -y11, -y21, -y31, rho))
						+ pow(S2, -eps)
						* (R3(-2 * eps, eps, eps, 1, y12, y22, y32, rho) + R3(-2 * eps, eps, eps, 1, -y12, -y22, -y32, rho));
					res += lsum / (a1[k]*b2[k] - a2[k]*b1[k]);

				}
			} else {
				res += 0;
			}
		}
		res *= -I * pow(Pi, 2) * pow(Pi, -eps) * tgamma(eps) * beta(2*eps,numeric(1));
		res = res.diff(M1, t1-1).diff(M2, t2-1).diff(M3, t3-1);
		res = subs(res, lst(M1 == pow(m1, 2), M2 == pow(m2, 2), M3 == pow(m3, 2)));
		c = res / (factorial(t1-1) * factorial(t2-1) * factorial(t3-1));
		
	} else {

		// This is the case where (q21==0 && q10 !=0).
		// I don't know whether the case q21==0 can appear in SM, but to keep
		// the calculation as general as possible, I try to consider all possibilities.

		if (q20.is_zero()) {

			// Reduce to two-point functions
			if (m2.is_equal(m3))
				c = OneLoop2Pt(0, 0, q10, m1, m2, t1, t2+t3, rho);
			else {
				ex A1 = diff(1/(M2-M3), M3, t3-1);
				ex A3 = diff(1/(M2-M3), M2, t2-1);
				ex d, e;
				for (int i=0; i<=t2-1; i++)
					d += binomial(t2-1, i) * factorial(i) * OneLoop2Pt(0, 0, q10, m1, m2, t1, i+1, rho) * diff(A1, M2, t2-1-i);
				for (int i=0; i<=t3-1; i++)
					e += binomial(t3-1, i) * factorial(i) * OneLoop2Pt(0, 0, q10, m1, m3, t1, i+1, rho) * diff(A3, M3, t3-1-i);
				c = subs(d-e, lst(M2 == pow(m2, 2), M3 == pow(m3,2))) / (factorial(t2-1) * factorial(t3-1));
			}

		} else if (q10.is_equal(q20)) {

			// Reduce to two-point functions
			if (m2.is_equal(m1))
				c = OneLoop2Pt(0, 0, -q10, m3, m1, t3, t2+t1, rho);
			else {
				ex A1 = diff(1/(M1-M2), M2, t2-1);
				ex A3 = diff(1/(M1-M2), M1, t1-1);
				ex d, e;
				for (int i=0; i<=t1-1; i++)
					d += binomial(t1-1, i) * factorial(i) * OneLoop2Pt(0, 0, -q10, m3, m1, t3, i+1, rho) * diff(A1, M1, t1-1-i);
				for (int i=0; i<=t2-1; i++)
					e += binomial(t2-1, i) * factorial(i) * OneLoop2Pt(0, 0, -q10, m3, m2, t3, i+1, rho) * diff(A3, M2, t2-1-i);
				c = subs(d-e, lst(M2 == pow(m2, 2), M1 == pow(m1,2))) / (factorial(t2-1) * factorial(t1-1));
			}

		} else { // q10 != q20, No more crazy singularity here
              
			ex res;
			for (int k=0; k<3 ; k++) {
				if ( !IR[k].is_zero() ){
					ex y11 = sqrt(m_2[k] - I*rho);
					ex y21 = -y11;
					ex y31 = c2[k]/a2[k] + I*rho;
					ex y32 = c1[k]/a1[k] + I*rho;
					ex F1 = 1 / (c1[k]*a2[k] - c2[k]*a1[k]);
					ex lsum = R3(-2*eps + 1, -numeric(1,2)+eps, -numeric(1,2)+ eps, 1, y11, y21, y31, rho)
						+ R3(-2*eps + 1, -numeric(1,2)+eps,-numeric(1,2)+ eps, 1, -y11, -y21, -y31, rho)
						- (R3(-2*eps + 1, -numeric(1,2)+eps, -numeric(1,2)+ eps, 1, y11, y21, y32, rho)
						   + R3(-2*eps + 1, -numeric(1,2)+eps,-numeric(1,2)+ eps, 1, -y11, -y21, -y32, rho));
					
					res += lsum*F1*beta(2*eps-1, 1);
				} else {
					res += 0;
				}
			}
			res *= pow(Pi, numeric(3,2)) * pow(Pi, -eps) * exp(I*Pi*eps) * tgamma(numeric(1,2)+eps) / (numeric(1,2)-eps);
			res = res.diff(M1, t1-1).diff(M2, t2-1).diff(M3, t3-1);
			res = subs(res, lst(M1 == pow(m1, 2), M2 == pow(m2, 2), M3 == pow(m3, 2)));
			c = res / (factorial(t1-1) * factorial(t2-1) * factorial(t3-1));
		}
	}
	// Expand into series in eps
	return c.series(r, order);
}

REGISTER_FUNCTION(Scalar3Pt, eval_func(Scalar3Pt_eval).
                             series_func(Scalar3Pt_series));

} // namespace xloops
