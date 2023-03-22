/** @file utils.cpp
 *
 *  Implementation of utility functions that are not supposed to be visible
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

#include "utils.h"
#include "po_redux.h"
#include "assertion.h"
namespace xloops {


ex Pcollect(const ex &A, int m, int n, int o, int r, const symbol &P1, const symbol &P2, const symbol &P3, const symbol &P4)
{
	return A.expand().coeff(P1, m).coeff(P2, n).coeff(P3, o).coeff(P4, r);
}

static po_redux & get_redux_object(void)
{
	static po_redux * redux = new po_redux;
	return *redux;
}

ex CoeffDim(int i, int j, int k, const ex &eps)
{
	po_redux_powers_vector v;
	ex Dim = 4 - 2*eps;
	ex Dim0 = 3 - 2*eps;
    
	if (k == -1) {
		v = get_redux_object().transform_1_to_0(i, j/2, 0, 0, 0, Dim);
		ASSERT(v.size() == 1);
	} else if (i == -1) {
		v = get_redux_object().transform_n_plus_1_to_n(j, k/2, 0, 0, 0, Dim0, 1);
		ASSERT(v.size() == 1);
	} else if (j == -1) {
		v = get_redux_object().transform_1_to_0(i, k/2, 0, 0, 0, Dim0);
		ASSERT(v.size() == 1);
	} else {
		v = get_redux_object().transform_2_to_0(i, j, k/2, 0, 0, 0, 0, Dim);
		ASSERT(v.size() == 1);
	}
	return v[0].coeff;
}

} // namespace xloops
