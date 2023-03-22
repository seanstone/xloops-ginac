/** @file po_redux.h
 *
 *  Definition of parallel space dimensional reduction procedure. */

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

#ifndef __XLOOPS_PO_REDUX_H__
#define __XLOOPS_PO_REDUX_H__

#include <vector>
#include <list>
#include <ginac/ginac.h>

namespace xloops {

using GiNaC::ex;
using GiNaC::symbol;
using GiNaC::exvector;
using GiNaC::indexed;
using GiNaC::scalar_products;


// Helper types definitions

struct po_redux_powers
{
	ex coeff;
	int pow_kq;
	int pow_lq;
	int pow_kl;
};

typedef std::vector<po_redux_powers> po_redux_powers_vector;

struct tensor_entry
{
	// store the generic orthogonal space tensor for
	// the integral k^\mu1...k^\mu(kpot) l^\nu1...l^\nu(lpot)
	// together with the indices with which it was generated
	// (indices will be replaced when needed)
	tensor_entry(int kp, int lp, ex t, const exvector & iv) :
	  kpot(kp), lpot(lp), tensor(t), idxv(iv) { }

	int kpot;
	int lpot;
	ex tensor;
	exvector idxv;
};

typedef std::vector<int> intvector;

class po_redux
{
public:
	po_redux();
	po_redux_powers_vector transform_1_to_0(int pow_k0_old, int pow_koq_old,
	                                        int pow_l0_old, int pow_loq_old,
	                                        int pow_koloz_old, ex Dim_new);
	po_redux_powers_vector transform_2_to_0(int pow_k0_old,
	                                        int pow_k1_old,
	                                        int pow_koq_old,
	                                        int pow_l0_old,
	                                        int pow_l1_old,
	                                        int pow_loq_old,
	                                        int pow_koloz_old,
	                                        ex Dim_new);
	po_redux_powers_vector transform_n_plus_1_to_n(int pow_k1_old,
	                                               int pow_koq_old,
	                                               int pow_l1_old,
	                                               int pow_loq_old,
	                                               int pow_koloz_old,
	                                               ex DimO_new,
	                                               int DimP_new);
	ex get_tensor(int kpot, int lpot, const exvector & idvx);

protected:
	bool is_perm_ok(const intvector & iv) const {
		// inlined, because it is called very often...
		int len=iv.size();
		for (int i=1; i<=len-3; i+=2) {
			if (iv[i-1]>iv[i]) return false;
			if (iv[i-1]>iv[i+1]) return false;
		}
		return iv[len-2]<iv[len-1];
	}
	int number_of_mixed_g(intvector iv, int kpot, int lpot) const;
	exvector generate_g_terms(const exvector & idxv,
	                          int kpot, int lpot) const;
	tensor_entry generate_tensor(int kpot, int lpot) const;
	ex evaluate_tensor_0(ex num, int pow_k, int pow_l);
	ex evaluate_tensor_1(ex num, int pow_k, int pow_l);
	po_redux_powers_vector collect_powers(const ex & res,
	                                      const symbol & kq_sym,
	                                      const symbol & lq_sym,
	                                      const symbol & kl_sym);
	void collect_powers_single_term(po_redux_powers_vector & v,
	                                const ex & res, const symbol & kq_sym,
	                                const symbol & lq_sym,
	                                const symbol & kl_sym);

protected:
	std::list<tensor_entry> tensors;

	symbol Dim;

	symbol k;
	symbol l;
	symbol kq;
	symbol lq;
	symbol kdotl;

	symbol k1_old, koq_old, l1_old, loq_old, koloz_old;
	symbol k1_new, koq_new, l1_new, loq_new, koloz_new;
	symbol k0_old, l0_old, k0_new, l0_new, kq_new, lq_new, kdotl_new;
	
	scalar_products sp;
};

} // namespace xloops

#endif // ndef __XLOOPS_PO_REDUX_H__
