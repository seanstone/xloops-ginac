/** @file one_loop_two_point.cpp
 *
 *  Check consistency of OneLoop2Pt() function. */

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

static unsigned one_loop_two_point1(void)
{
	unsigned result = 0;
	symbol eps("eps"), rho("rho"), q("q"), m1("m1"), m2("m2");

	// Cancellation of q-less propagator reduces to a 1-point function
	ex a = OneLoop2Pt(2, 0, q, m1, m2, 1, 1, rho) - OneLoop2Pt(0, 2, q, m1, m2, 1, 1, rho) - pow(m2, 2) * OneLoop2Pt(0, 0, q, m1, m2, 1, 1, rho);
	a -= OneLoop1Pt(0, m1, 1, rho);
	a = tensor_reduction(a, eps);

	if (!a.match(rho * wild())) {
		clog << "Function should be linear in rho, but is " << a << endl;
		++result;
	}

	return result;
}

static unsigned one_loop_two_point2(void)
{
	unsigned result = 0;
	symbol eps("eps"), rho("rho"), q("q"), m("m");

	// Transversality of vacuum polarization
	ex a = OneLoop2Pt(2, 0, q, m, m, 1, 1, rho) + OneLoop2Pt(0, 2, q, m, m, 1, 1, rho)
		 + q * OneLoop2Pt(1, 0, q, m, m, 1, 1, rho) + pow(m, 2) * OneLoop2Pt(0, 0, q, m, m, 1, 1, rho);
	a = tensor_reduction(a, eps);

	if (!a.match(rho * wild())) {
		clog << "Function should be linear in rho, but is " << a << endl;
		++result;
	}

	return result;
}

unsigned one_loop_two_point(void)
{
	unsigned result = 0;

	cout << "checking consistency of OneLoop2Pt" << flush;
	clog << "---------consistency of OneLoop2Pt:" << endl;

	result += one_loop_two_point1();  cout << '.' << flush;
	result += one_loop_two_point2();  cout << '.' << flush;

	if (! result) {
		cout << " passed ";
		clog << "(no output)" << endl;
	} else {
		cout << " failed ";
	}

	return result;
}
