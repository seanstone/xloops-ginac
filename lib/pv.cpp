/** @file pv.cpp
 *
 *  Implementation of PV tensor decomposition. */

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
#include <sstream>

#include "oneloop.h"
using namespace GiNaC;

namespace xloops {


/** This main function of this class is to store all the indices and
 *  coefficient symbols that are used in the algorithms (via the factory
 *  methods get_index() and get_c?pt()). */
class pv {
public:
	pv() {}
	~pv() {}

	// Symmetrize removing constant factor
	ex pv_symmetrize(const ex & e)
	{
		ex r = symmetrize(e);
		if (is_a<add>(r))
			return r * r.nops();
		else
			return r;
	}

	// Return general Lorentz tensor of rank n with no external momenta
	ex general_tensor_1pt(unsigned n)
	{
		if (n & 1)
			throw std::invalid_argument("general_tensor_1pt: parameter must be even");
		if (n == 0)
			return get_c1pt(0);

		ex res = 1;
		for (int i=0; i<n/2; i++)
			res *= lorentz_g(get_index(2*i), get_index(2*i+1));
		if (n > 2)
			res = pv_symmetrize(res);
		return res * get_c1pt(n);
	}

	// Return general Lorentz tensor of rank n with one external momentum
	ex general_tensor_2pt(unsigned n)
	{
		ex res;
		for (int i=0; i<=n/2; i++) {
			ex fac = 1;
			for (int j=0; j<i; j++)
				fac *= lorentz_g(get_index(2*j), get_index(2*j+1));
			for (int j=2*i; j<n; j++)
				fac *= indexed(get_q(1), get_index(j));
			if (n > 2)
				fac = pv_symmetrize(fac);
			res += fac * get_c2pt(n, i);
		}
		return res;
	}

	// Return general Lorentz tensor of rank n with two external momenta
	ex general_tensor_3pt(unsigned n)
	{
		ex res;
		for (int i=0; i<=n/2; i++) {
			ex fac = 1;
			for (int j=0; j<i; j++)
				fac *= lorentz_g(get_index(2*j), get_index(2*j+1));
			for (int k=0; k<=n-2*i; k++) {
				ex fac2 = fac;
				for (int l=2*i; l<n-k; l++)
					fac2 *= indexed(get_q(2), get_index(l));
				for (int l=n-k; l<n; l++)
					fac2 *= indexed(get_q(1), get_index(l));
				if (n > 1)
					fac2 = pv_symmetrize(fac2);
				res += fac2 * get_c3pt(n, i, k);
			}
		}
		return res;
	}

	// Get normalizing constant for PV tensor given the rank and dimension
	// of orthogonal space
	ex norm_const(unsigned n_orth, const ex & D_orth)
	{
		if (n_orth == 0)
			return 1;
		else if (n_orth & 1)
			return 0;

		int num = 1;
		ex den = 1;
		for (int i=0; i<=(n_orth-2)/2; i++) {
			num *= n_orth - 2*i - 1;
			den *= D_orth + 2*i;
		}
		return num / den.expand();
	}

	// Construct linear combination of OneLoop1Pt() functions for tensor of
	// given (PO-decomposed) rank
	ex pv_tensor_1pt(unsigned n_orth, const ex & m, const ex & t, const ex & rho, const ex & eps)
	{
		return OneLoop1Pt(n_orth, m, t, rho) * norm_const(n_orth, 4 - 2*eps);
	}

	// Construct linear combination of OneLoop2Pt() functions for tensor of
	// given (PO-decomposed) rank
	ex pv_tensor_2pt(unsigned n_0, unsigned n_orth, const ex & q, const ex & m1, const ex & m2, const ex & t1, const ex & t2, const ex & rho, const ex & eps)
	{
		return OneLoop2Pt(n_0, n_orth, q, m1, m2, t1, t2, rho) * norm_const(n_orth, 3 - 2*eps);
	}

	// Construct linear combination of OneLoop3Pt() functions for tensor of
	// given (PO-decomposed) rank
	ex pv_tensor_3pt(unsigned n_0, unsigned n_1, unsigned n_orth, const ex & q10, const ex & q20, const ex & q21, const ex & m1, const ex & m2, const ex & m3, const ex & t1, const ex & t2, const ex & t3, const ex & rho, const ex & eps)
	{
		return OneLoop3Pt(n_0, n_1, n_orth, q10, q20, q21, m1, m2, m3, t1, t2, t3, rho) * norm_const(n_orth, 2 - 2*eps);
	}

	// Return index number "n" (either create a new one for this number,
	// or return a previously created one)
	ex get_index(unsigned n)
	{
		std::map<unsigned, ex>::iterator it = indices.find(n);
		if (it != indices.end())
			return it->second;
		std::ostringstream s;
		s << "mu" << n << std::ends;
		return indices.insert(std::map<unsigned, ex>::value_type(n, varidx(symbol(s.str()), 4))).first->second;
	}

	// Return momentum number "n" (either create a new one for this number,
	// or return a previously created one)
	ex get_q(unsigned n)
	{
		std::map<unsigned, ex>::iterator it = momenta.find(n);
		if (it != momenta.end())
			return it->second;
		std::ostringstream s;
		s << "q" << n << std::ends;
		return momenta.insert(std::map<unsigned, ex>::value_type(n, symbol(s.str()))).first->second;
	}

	// Return 1-point coefficient "n" (either create a new one for this number,
	// or return a previously created one)
	ex get_c1pt(unsigned n)
	{
		std::map<unsigned, ex>::iterator it = coeffs.find(n);
		if (it != coeffs.end())
			return it->second;
		std::ostringstream s;
		s << "C" << n << std::ends;
		return coeffs.insert(std::map<unsigned, ex>::value_type(n, symbol(s.str()))).first->second;
	}

	// Return 2-point coefficient "n,m"
	ex get_c2pt(unsigned n, unsigned m)
	{
		return get_c1pt(n*10 + m);
	}

	// Return 3-point coefficient "n,m,l"
	ex get_c3pt(unsigned n, unsigned m, unsigned l)
	{
		return get_c1pt(n*100 + m*10 + l);
	}

private:
	std::map<unsigned, ex> indices;
	std::map<unsigned, ex> momenta;
	std::map<unsigned, ex> coeffs;
};


// One-loop one-point tensor decomposition
ex OneLoopTens1Pt(unsigned n, const ex & m, const ex & t, const ex & rho, const ex & eps, lst & Cl)
{
	Cl = lst();

	if (n & 1)
		return 0;

	pv P;

	// General tensor structure
	ex g = P.general_tensor_1pt(n);

	// Extract one tensor component (replace all eta~mu~nu by eta~0~0 (= 1))
	ex c = g.subs(lorentz_g(varidx(wild(1), 4), varidx(wild(2), 4)) == 1);

	// Decompose into OneLoop1Pt() functions, solve for Cx
	ex a = c.coeff(P.get_c1pt(n), 1);
	ex f = P.pv_tensor_1pt(n, m, t, rho, eps) / a;

	Cl.append(P.get_c1pt(n) == f);
	return g;
}


// One-loop two-point tensor decomposition
ex OneLoopTens2Pt(unsigned n, const ex & q0, const ex & m1, const ex & m2, const ex & t1, const ex & t2, const ex & rho, const ex & eps, lst & Cl)
{
	Cl = lst();
	pv P;

	// General tensor structure
	ex g = P.general_tensor_2pt(n);

	// Momentum components
	lst ql;
	ql.append(indexed(P.get_q(1), varidx(0, 4)) == q0);
	ql.append(indexed(P.get_q(1), varidx(1, 4)) == 0);

	// Compute the coefficients
	lst E;
	for (int i1=n/2; i1>=0; i1--) {

		// Extract one tensor component
		lst sl;
		int i2;
		for (i2=0; i2<2*i1; i2++)
			sl.append(P.get_index(i2) == 1);
		for (; i2<n; i2++)
			sl.append(P.get_index(i2) == 0);
		ex c = g.subs(sl).subs(ql);

		// Decompose into OneLoop2Pt functions, solve for Cxy
		// (Cxy * a + b = pv_tensor_2pt -> Cxy = (pv_tensor_2pt - b) / a)
		ex a = c.coeff(P.get_c2pt(n, i1), 1);
		ex b = c.coeff(P.get_c2pt(n, i1), 0);
		ex f = (P.pv_tensor_2pt(n-2*i1, 2*i1, q0, m1, m2, t1, t2, rho, eps) - b) / a;
		f = f.subs(E);
		E.append(P.get_c2pt(n, i1) == f);

		// Collect by OneLoop2Pt() functions
		lst fcns;
		f.find(OneLoop2Pt(wild(0), wild(1), q0, m1, m2, t1, t2, rho), fcns);
		Cl.append(P.get_c2pt(n, i1) == f.expand().collect(fcns));
	}

	return g;
}


// One-loop three-point tensor decomposition
ex OneLoopTens3Pt(unsigned n, const ex & q10, const ex & q20, const ex & q21, const ex & m1, const ex & m2, const ex & m3, const ex & t1, const ex & t2, const ex & t3, const ex & rho, const ex & eps, lst & Cl)
{
	Cl = lst();
	pv P;

	// General tensor structure
	ex g = P.general_tensor_3pt(n);

	// Momentum components
	lst ql;
	ql.append(indexed(P.get_q(1), varidx(0, 4)) == q10);
	ql.append(indexed(P.get_q(1), varidx(1, 4)) == 0);
	ql.append(indexed(P.get_q(1), varidx(2, 4)) == 0);
	ql.append(indexed(P.get_q(2), varidx(0, 4)) == q20);
	ql.append(indexed(P.get_q(2), varidx(1, 4)) == q21);
	ql.append(indexed(P.get_q(2), varidx(2, 4)) == 0);

	// Compute the coefficients
	lst E;
	for (int i1=n/2; i1>=0; i1--) {
		for (int i2=n; i2>=2*i1; i2--) {

			// Extract one tensor component
			lst sl;
			int i3;
			for (i3=0; i3<2*i1; i3++)
				sl.append(P.get_index(i3) == 2);
			for (; i3<i2; i3++)
				sl.append(P.get_index(i3) == 1);
			for (; i3<n; i3++)
				sl.append(P.get_index(i3) == 0);
			ex c = g.subs(sl).subs(ql);

			// Decompose into OneLoop3Pt functions, solve for Cxyz
			// (Cxyz * a + b = pv_tensor_3pt -> Cxyz = (pv_tensor_3pt - b) / a)
			ex a = c.coeff(P.get_c3pt(n, i1, n-i2), 1);
			ex b = c.coeff(P.get_c3pt(n, i1, n-i2), 0);
			ex f = (P.pv_tensor_3pt(n-i2, i2-2*i1, 2*i1, q10, q20, q21, m1, m2, m3, t1, t2, t3, rho, eps) - b) / a;
			f = f.subs(E);
			E.append(P.get_c3pt(n, i1, n-i2) == f);

			// Collect by OneLoop3Pt() functions
			lst fcns;
			f.find(OneLoop3Pt(wild(0), wild(1), wild(2), q10, q20, q21, m1, m2, m3, t1, t2, t3, rho), fcns);
			Cl.append(P.get_c3pt(n, i1, n-i2) == f.expand().collect(fcns));
		}
	}

	return g;
}

} // namespace xloops
