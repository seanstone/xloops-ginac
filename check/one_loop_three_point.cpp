/** @file one_loop_three_point.cpp
 *
 *  Check consistency of OneLoop3Pt() function. */

/*
 *  xloops Copyright (C) 1997,2000 Johannes Gutenberg University Mainz, Germany
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

#include <iostream>
#include "xloops.h"

using namespace std;
using namespace GiNaC;
using namespace xloops;

static unsigned one_loop_three_point1(void)
{
	unsigned result = 0;
	symbol eps("eps"), rho("rho"), q10("q10"), q20("q20"), q21("q21");
	symbol m1("m1"), m2("m2"), m3("m3");
	
	// Cancellation of q-less propagator reduces to a 2-point function with 2-dimensional parallel space
	ex a = OneLoop3Pt(2, 0, 0, q10, q20, q21, m1, m2, m3, 1, 1, 1, rho)
		 - OneLoop3Pt(0, 2, 0, q10, q20, q21, m1, m2, m3, 1, 1, 1, rho)
		 - OneLoop3Pt(0, 0, 2, q10, q20, q21, m1, m2, m3, 1, 1, 1, rho)
		 - pow(m3, 2) * OneLoop3Pt(0, 0, 0, q10, q20, q21, m1, m2, m3, 1, 1, 1, rho);
	a -= S322Pt(0, 0, 0, q10, q20, q21, m1, m2, 1, 1, rho);
	a = tensor_reduction(a, eps);

	if (!a.match(rho * wild())) {
		clog << "Function should be linear in rho, but is " << a << endl;
		++result;
	}

	return result;
}

static unsigned one_loop_three_point2(void)
{
	unsigned result = 0;
	symbol eps("eps"), rho("rho"), q10("q10"), q20("q20"), q21("q21");
	symbol m1("m1"), m2("m2"), m3("m3");
	
	// Scalar three-point function is finite
	ex a = OneLoop3Pt(0, 0, 0, q10, q20, q21, m1, m2, m3, 1, 1, 1, rho);
	a = a.series(eps == 0, 2);
	
	ex a1 = a.coeff(eps, -2).normal();
	ex a2 = a.coeff(eps, -1).normal();
	
	if (!a1.is_zero()) {
		clog << "Order eps^-2 is " << a1 << ", but should be zero\n";
		++result;
	}
	if (!a2.is_zero()) {
		clog << "Order eps^-1 is " << a2 << ", but should be zero\n";
		++result;
	}
	
	return result;
}

static unsigned one_loop_three_point3(void)
{
	unsigned result = 0;
	symbol eps("eps"), rho("rho"), q("q"), p0("p0"), p1("p1"), m("m"), e("e"), g0("g0"), g1("g1"), gorth("gorth");

	// Relation of vertex function to difference of self energies (Ward identity)
	ex T000 = OneLoop3Pt(0, 0, 0, q, -p0, -p1, m, 0, m, 1, 1, 1, rho);
	ex T200 = OneLoop3Pt(2, 0, 0, q, -p0, -p1, m, 0, m, 1, 1, 1, rho);
	ex T020 = OneLoop3Pt(0, 2, 0, q, -p0, -p1, m, 0, m, 1, 1, 1, rho);
	ex T002 = OneLoop3Pt(0, 0, 2, q, -p0, -p1, m, 0, m, 1, 1, 1, rho);
	ex T100 = OneLoop3Pt(1, 0, 0, q, -p0, -p1, m, 0, m, 1, 1, 1, rho);
	ex T010 = OneLoop3Pt(0, 1, 0, q, -p0, -p1, m, 0, m, 1, 1, 1, rho);
	ex T001 = OneLoop3Pt(0, 0, 1, q, -p0, -p1, m, 0, m, 1, 1, 1, rho);
	ex T110 = OneLoop3Pt(1, 1, 0, q, -p0, -p1, m, 0, m, 1, 1, 1, rho);
	ex T101 = OneLoop3Pt(1, 0, 1, q, -p0, -p1, m, 0, m, 1, 1, 1, rho);

	ex D = 4 - 2*eps;
	ex G = -pow(e,3) * (D*m*pow(q,2)*T000
		 + q*g0*(pow(m,2)*(2-D)*T000-(2-D)*(T200-T020-T002))
		 + 2*D*m*q*T100
		 + pow(q,2)*(2-D)*(g0*T100-g1*T010-gorth*T001)
		 + 2*(2-D)*q*(g0*T200-g1*T110-gorth*T101));
	ex g = tensor_reduction(G, eps);

	ex Tx1 = S322Pt(0, 0, 0, 0, p0, p1, 0, m, 1, 1, rho);
	ex Tx2 = S322Pt(1, 0, 0, 0, p0, p1, 0, m, 1, 1, rho);
	ex Tx3 = S322Pt(0, 1, 0, 0, p0, p1, 0, m, 1, 1, rho);
	ex Tx4 = S322Pt(0, 0, 1, 0, p0, p1, 0, m, 1, 1, rho);
	ex Ty1 = S322Pt(0, 0, 0, 0, p0+q, p1, 0, m, 1, 1, rho);
	ex Ty2 = S322Pt(1, 0, 0, 0, p0+q, p1, 0, m, 1, 1, rho);
	ex Ty3 = S322Pt(0, 1, 0, 0, p0+q, p1, 0, m, 1, 1, rho);
	ex Ty4 = S322Pt(0, 0, 1, 0, p0+q, p1, 0, m, 1, 1, rho);

	ex Sx = -pow(e,2) * (D*m*Tx1 + (2-D)*(p0*g0-p1*g1)*Tx1
		  + (2-D)*(g0*Tx2-g1*Tx3-gorth*Tx4));
	ex Sy = -pow(e,2) * (D*m*Ty1 + (2-D)*((p0+q)*g0-p1*g1)*Ty1
		  + (2-D)*(g0*Ty2-g1*Ty3-gorth*Ty4));
	ex s = (e * (Sx - Sy));

	ex diff = tensor_reduction(g - s, eps);

	// We're only interested in the numerator
	diff = diff.expand(expand_options::expand_function_args).numer();

	// Re-collect by scalar functions to make the series expansion easier
	diff = tensor_reduction(diff, eps);

	// Expand into Laurent series
	diff = diff.series(eps == 0, 4).subs(rho == 0);

	ex c1 = expand(coeff(diff, eps, -1));
	ex c2 = expand(coeff(diff, eps, 0));
	ex c3 = expand(coeff(diff, eps, 1));

	if (!c1.is_zero()) {
		clog << "Order eps^-1 is " << c1 << ", but should be zero\n";
		++result;
	}
	if (!c2.is_zero()) {
		clog << "Order eps^0 is " << c2 << ", but should be zero\n";
		++result;
	}
	if (!c3.is_zero()) {
		clog << "Order eps^1 is " << c3 << ", but should be zero\n";
		++result;
	}

	return result;
}

unsigned one_loop_three_point(void)
{
	unsigned result = 0;
	
	cout << "checking consistency of OneLoop3Pt" << flush;
	clog << "---------consistency of OneLoop3Pt:" << endl;
	
	result += one_loop_three_point1();  cout << '.' << flush;
	result += one_loop_three_point2();  cout << '.' << flush;
	result += one_loop_three_point3();  cout << '.' << flush;
	
	if (! result) {
		cout << " passed ";
		clog << "(no output)" << endl;
	} else {
		cout << " failed ";
	}
	
	return result;
}
