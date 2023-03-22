/** @file po_redux.cpp
 *
 *  Implementation of parallel space dimensional reduction procedure. */

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

#include <iostream>
#include <algorithm>

#include "po_redux.h"
#include "assertion.h"

using namespace std;
using namespace GiNaC;

namespace xloops {

/** Substitute multiple indices in an expression.
 *
 *  @param e Expression to substitute in
 *  @param idxv_subs Vector of indices being substituted
 *  @param idxv_repl Vector of indices to replace by (1:1 correspondence to idxv_subs)
 *  @return expression with substituted indices */
static ex subs_indices(const ex & e, const exvector & idxv_subs, const exvector & idxv_repl)
{
	GINAC_ASSERT(idxv_subs.size() == idxv_repl.size());
	ex res=e;
	for (unsigned i=0; i<idxv_subs.size(); ++i)
		res = res.subs(idxv_subs[i] == idxv_repl[i]);
	return res;
}


// default constructor

po_redux::po_redux() :
	Dim("Dim"),
	k("k"), l("l"),
	kq("kq"), lq("lq"), kdotl("kdotl"),
	k1_old("k1'"), koq_old("koq'"), l1_old("l1'"), loq_old("loq'"), koloz_old("koloz'"),
	k1_new("k1"), koq_new("koq"), l1_new("l1"), loq_new("loq"), koloz_new("koloz"),
	k0_old("k0'"), l0_old("l0'"), k0_new("k0"), l0_new("l0"),
	kq_new("kq"), lq_new("lq"), kdotl_new("kdotl")
{
	// register scalar products kv.kv=kq, lv.lv=lq, kv.lv=kdotl
	sp.add(k, k, kq);
	sp.add(l, l, lq);
	sp.add(k, l, kdotl);

	// add the kpot=0, lpot=0 tensor as a special case
	tensors.push_back(tensor_entry(0,
								   0,
								   1,
								   exvector()));
}    

int po_redux::number_of_mixed_g(vector<int> iv, int kpot, int lpot) const
{
	int degree=kpot+lpot;
	int nmg=0;
	for (int i=0; i<=degree-2; ++i) {
		if ((iv[i]<kpot)&&(iv[i+1]>=kpot)) nmg++;
	}
	return nmg;
}

exvector po_redux::generate_g_terms(const exvector & idxv,
									int kpot, int lpot) const
{
	// generate all g^{\mu \nu} terms
	// input indices are in idxv
	// idxv[0]..idxv[kpot-1]:         indices for k^\mu
	// idxv[kpot]..idxv[kpot+lpot.1]: indices for l^\nu
	// the n-th entry in the return list is a sum of all g^{\mu\nu} terms
	// which mix n k/l indices
	
	int degree=kpot+lpot;
	
	ASSERT(degree % 2 == 0);
	ASSERT(degree<=(int)idxv.size());

	vector<int> iv;
	unsigned count=0;
	exvector res;
	res.resize(degree/2+1);

	// number of generated terms is (degree-1)!!
	unsigned refsize=1;
	for (int i=1; i<=degree-1; i+=2) {
		refsize *= i;
	}

	// init permutation vector
	iv.resize(degree);
	for (int i=0; i<degree; ++i) {
		iv[i]=i;
	}

	// cycle through all (kpot+lpot)! permutations
	do {
		if (is_perm_ok(iv)) {
			// consider only "sorted" terms (g^{\nu\mu} -> g^{\mu\nu})
			count ++;
			ex term=1;
			for (int i=0; i<=degree-2; i+=2) {
				term *= lorentz_g(idxv[iv[i]],idxv[iv[i+1]]);
			}
			int nmg=number_of_mixed_g(iv,kpot,lpot);
			res[nmg] += term;
		}
	} while (next_permutation(iv.begin(), iv.end()));
		
	// test number of terms
	ASSERT(count==refsize);

	return res;
}

tensor_entry po_redux::generate_tensor(int kpot, int lpot) const
{
	// Creates an ansatz for the generic orthogonal space tensor
	// for the integral k^\mu1...k^\mu(kpot) l^\nu1...l^\nu(lpot).
	// The PO-space dimension D=g^\mu_\mu is replaced later
	// (in evaluate_tensor0/1) with the correct dimension for the
	// considered case.
	// The g^{\mu\nu} terms from generate_g_terms are multiplied
	// with an undetermined coefficient, contracted and the resulting
	// set of linear equations is solved for the coefficients.
	
	int degree=kpot+lpot;

	ASSERT(degree>0);
	ASSERT(degree % 2 == 0);
	
	exvector idxv_contra;
	exvector idxv_co;
	idxv_contra.reserve(degree);
	idxv_co.reserve(degree);
	for (int i=0; i<degree; ++i) {
		symbol s;
		varidx ix(s, Dim);
		idxv_contra.push_back(ix);
		idxv_co.push_back(ix.toggle_variance());
	}

	exvector g_terms=generate_g_terms(idxv_contra,kpot,lpot);

	ex num=1;
	for (int i=0; i<kpot; ++i) {
		num *= indexed(k,idxv_contra[i]);
	}
	for (int i=kpot; i<degree; ++i) {
		num *= indexed(l,idxv_contra[i]);
	}

	symbol * varvec=new symbol[g_terms.size()];

	lst eqns;
	lst vars;
	ex ansatz;
	for (int i=0; i<(int)g_terms.size(); ++i) {
		//cout << "g_terms[" << i << "]=" << g_terms[i] << endl; 
		ansatz += (*(varvec+i))*g_terms[i];
	}
	// cout << "ansatz=" << ansatz << endl;
	for (int i=0; i<(int)g_terms.size(); ++i) {
		ex contract_element=g_terms[i];
		if (is_exactly_a<add>(contract_element)) {
			contract_element=contract_element.op(0);
		}
		contract_element=subs_indices(contract_element,idxv_contra,idxv_co);
		ex lhs=simplify_indexed(num*contract_element,sp);
		//cout << "lhs1: " << num*g_terms[i] << endl;
		//cout << "lhs2: " << lhs << endl;
		ex rhs=simplify_indexed(ansatz*contract_element,sp);
		//cout << "rhs1: " << ansatz*g_terms[i] << endl;
		//cout << "rhs2: " << rhs << endl;
		eqns.append(lhs==rhs);
		vars.append(*(varvec+i));
	}
	// cout << "eqns=" << eqns << endl;
	// cout << "vars=" << vars << endl;
	ex sol=lsolve(eqns,vars);
	// cout << "sol=" << sol << endl;
	ex sol2=ansatz.subs(sol).normal();
	// cout << "tensor(" << kpot <<" ," << lpot << ")=" << sol2 << endl;
	
	delete [] varvec;
	return tensor_entry(kpot,lpot,sol2,idxv_contra);
}

ex po_redux::get_tensor(int kpot, int lpot, const exvector & idxv)
{
	// search in list
	for (list<tensor_entry>::const_iterator cit=tensors.begin();
		 cit!=tensors.end(); ++cit) {
		if ((cit->kpot==kpot)&&(cit->lpot==lpot)) {
			return subs_indices(cit->tensor,cit->idxv,idxv);
		}
	}

	// not found: generate it and append to list
	tensors.push_back(generate_tensor(kpot,lpot));
	return subs_indices(tensors.rbegin()->tensor,tensors.rbegin()->idxv,idxv);
}

po_redux_powers_vector po_redux::transform_1_to_0(int pow_k0_old,
												  int pow_koq_old,
												  int pow_l0_old,
												  int pow_loq_old,
												  int pow_koloz_old,
												  ex Dim_new)
{
	// transform
	// k_0^pow_k0_old*(k_\perp^2)^pow_koq_old
	// *l_0^pow_l0_old*(l_\perp^2)^pow_loq_old
	// *(k_\perp l_\perp z)^pow_koloz_old
	// to \sum (k^2)^r (l^2)^s (k.l)^t

	ex num=power(k0_old,pow_k0_old)*power(koq_old,pow_koq_old)*
		   power(l0_old,pow_l0_old)*power(loq_old,pow_loq_old)*
		   power(koloz_old,pow_koloz_old);

	ex num2=num.subs(lst(k0_old==k0_new,
						 koq_old==k0_new*k0_new-kq_new,
						 l0_old==l0_new,
						 loq_old==l0_new*l0_new-lq_new,
						 koloz_old==k0_new*l0_new-kdotl_new)).expand();

	int pow_k=pow_k0_old+2*pow_koq_old+pow_koloz_old;
	int pow_l=pow_l0_old+2*pow_loq_old+pow_koloz_old;
	// cout << "num2=" << num2 << endl;
	ex res;
	if (is_exactly_a<add>(num2)) {
		for (int i=0; i<num2.nops(); ++i) {
			res += evaluate_tensor_0(num2.op(i),pow_k,pow_l);
		}
	} else {
		res=evaluate_tensor_0(num2,pow_k,pow_l);
	}
	res = res.subs(Dim == Dim_new).expand();
	return collect_powers(res,kq_new,lq_new,kdotl_new);
}

po_redux_powers_vector po_redux::transform_2_to_0(int pow_k0_old,
												  int pow_k1_old,
												  int pow_koq_old,
												  int pow_l0_old,
												  int pow_l1_old,
												  int pow_loq_old,
												  int pow_koloz_old,
												  ex Dim_new)
{
	// transform
	// k_0^pow_k0_old*k_1^pow_k1_old*(k_\perp^2)^pow_koq_old
	// *l_0^pow_l0_old*l_1^pow_l1_old*(l_\perp^2)^pow_loq_old
	// *(k_\perp l_\perp z)^pow_koloz_old
	// to \sum (k^2)^r (l^2)^s (k.l)^t

	ex num=power(k0_old,pow_k0_old)*power(k1_old,pow_k1_old)*power(koq_old,pow_koq_old)*
		   power(l0_old,pow_l0_old)*power(l1_old,pow_l1_old)*power(loq_old,pow_loq_old)*
		   power(koloz_old,pow_koloz_old);

	ex num2=num.subs(lst(k0_old==k0_new,
						 k1_old==k1_new,
						 koq_old==k0_new*k0_new-k1_new*k1_new-kq_new,
						 l0_old==l0_new,
						 l1_old==l1_new,
						 loq_old==l0_new*l0_new-l1_new*l1_new-lq_new,
						 koloz_old==k0_new*l0_new-k1_new*l1_new-kdotl_new)).expand();

	int pow_k=pow_k0_old+pow_k1_old+2*pow_koq_old+pow_koloz_old;
	int pow_l=pow_l0_old+pow_l1_old+2*pow_loq_old+pow_koloz_old;
	// cout << "num2=" << num2 << endl;
	ex res;
	if (is_exactly_a<add>(num2)) {
		for (int i=0; i<num2.nops(); ++i) {
			res += evaluate_tensor_0(num2.op(i),pow_k,pow_l);
		}
	} else {
		res=evaluate_tensor_0(num2,pow_k,pow_l);
	}
	res = res.subs(Dim == Dim_new).expand();
	return collect_powers(res,kq_new,lq_new,kdotl_new);
}

po_redux_powers_vector po_redux::transform_n_plus_1_to_n(int pow_k1_old,
														 int pow_koq_old,
														 int pow_l1_old,
														 int pow_loq_old,
														 int pow_koloz_old,
														 ex DimO_new,
														 int DimP_new)
{
	// transform
	// k_1^pow_k1_old*(k_\perp^2)^pow_koq_old
	// *l_1^pow_l1_old*(l_\perp^2)^pow_loq_old
	// *(k_\perp l_\perp z)^pow_koloz_old
	// to \sum (\tilde{k}_\perp^2)^r (\tilde{l}_\perp^2)^s
	//         (\tilde{k}_\perp \tilde{l}_\perp \tilde{z})^t

	ex num=power(k1_old,pow_k1_old)*power(koq_old,pow_koq_old)*
		   power(l1_old,pow_l1_old)*power(loq_old,pow_loq_old)*
		   power(koloz_old,pow_koloz_old);
	
	ex num2=num.subs(lst(k1_old==k1_new,
						 koq_old==koq_new-k1_new*k1_new,
						 l1_old==l1_new,
						 loq_old==loq_new-l1_new*l1_new,
						 koloz_old==koloz_new-k1_new*l1_new)).expand();
	
	int pow_k=pow_k1_old+2*pow_koq_old+pow_koloz_old;
	int pow_l=pow_l1_old+2*pow_loq_old+pow_koloz_old;
	// cout << "num2=" << num2 << endl;
	ex res;
	if (is_exactly_a<add>(num2)) {
		for (int i=0; i<num2.nops(); ++i) {
			res += evaluate_tensor_1(num2.op(i),pow_k,pow_l);
		}
	} else {
		res=evaluate_tensor_1(num2,pow_k,pow_l);
	}
	res = res.subs(Dim == DimO_new+DimP_new).expand();
	return collect_powers(res,koq_new,loq_new,koloz_new);
}

ex po_redux::evaluate_tensor_0(ex num, int pow_k, int pow_l)
{
	exvector idxv;
	idxv.resize(pow_k+pow_l);

	int pow_k0_new=num.degree(k0_new); 
	int pow_k1_new=num.degree(k1_new); 
	int pow_kq_new=num.degree(kq_new); 
	int pow_l0_new=num.degree(l0_new); 
	int pow_l1_new=num.degree(l1_new); 
	int pow_lq_new=num.degree(lq_new); 
	int pow_kdotl_new=num.degree(kdotl_new);

	ex num2=power(k0_new,pow_k0_new)*power(k1_new,pow_k1_new)*power(kq_new,pow_kq_new)*
			power(l0_new,pow_l0_new)*power(l1_new,pow_l1_new)*power(lq_new,pow_lq_new)*
			power(kdotl_new,pow_kdotl_new);
	ex coeff=normal(num/num2);

	ASSERT(is_exactly_a<numeric>(coeff));
	ASSERT(pow_k==pow_k0_new+pow_k1_new+2*pow_kq_new+pow_kdotl_new);
	ASSERT(pow_l==pow_l0_new+pow_l1_new+2*pow_lq_new+pow_kdotl_new);

	for (int i=0; i<pow_k0_new; ++i) {
		idxv[i]=varidx(0,Dim);
	}
	for (int i=0; i<pow_k1_new; ++i) {
		idxv[pow_k0_new+i]=varidx(1,Dim);
	}
	for (int i=0; i<pow_kq_new; ++i) {
		symbol s_mu, s_nu;
		varidx mu(s_mu,Dim);
		varidx nu(s_nu,Dim);
		idxv[pow_k0_new+pow_k1_new+2*i]=mu;
		idxv[pow_k0_new+pow_k1_new+2*i+1]=nu;
		coeff *= lorentz_g(mu.toggle_variance(),nu.toggle_variance());
	}
	for (int i=0; i<pow_l0_new; ++i) {
		idxv[pow_k+i]=varidx(0,Dim);
	}
	for (int i=0; i<pow_l1_new; ++i) {
		idxv[pow_k+pow_l0_new+i]=varidx(1,Dim);
	}
	for (int i=0; i<pow_lq_new; ++i) {
		symbol s_mu, s_nu;
		varidx mu(s_mu,Dim);
		varidx nu(s_nu,Dim);
		idxv[pow_k+pow_l0_new+pow_l1_new+2*i]=mu;
		idxv[pow_k+pow_l0_new+pow_l1_new+2*i+1]=nu;
		coeff *= lorentz_g(mu.toggle_variance(),nu.toggle_variance());
	}
	for (int i=0; i<pow_kdotl_new; ++i) {
		symbol s_mu, s_nu;
		varidx mu(s_mu,Dim);
		varidx nu(s_nu,Dim);
		idxv[pow_k0_new+pow_k1_new+2*pow_kq_new+i]=mu;
		idxv[pow_k+pow_l0_new+pow_l1_new+2*pow_lq_new+i]=nu;
		coeff *= lorentz_g(mu.toggle_variance(),nu.toggle_variance());
	}
	ex t=get_tensor(pow_k,pow_l,idxv).subs(lst(kq==kq_new,
											   lq==lq_new,
											   kdotl==kdotl_new));
	//cout << "t=" << t << endl;
	ex t2=simplify_indexed(coeff*t,sp);
	//cout << "t2=" << t2 << endl;
	return t2;
}       

ex po_redux::evaluate_tensor_1(ex num, int pow_k, int pow_l)
{
	exvector idxv;
	idxv.resize(pow_k+pow_l);

	int pow_k1_new=num.degree(k1_new); 
	int pow_koq_new=num.degree(koq_new); 
	int pow_l1_new=num.degree(l1_new); 
	int pow_loq_new=num.degree(loq_new); 
	int pow_koloz_new=num.degree(koloz_new);

	ex num2=power(k1_new,pow_k1_new)*power(koq_new,pow_koq_new)*
			power(l1_new,pow_l1_new)*power(loq_new,pow_loq_new)*
			power(koloz_new,pow_koloz_new);
	ex coeff=normal(num/num2);

	ASSERT(is_exactly_a<numeric>(coeff));
	ASSERT(pow_k==pow_k1_new+2*pow_koq_new+pow_koloz_new);
	ASSERT(pow_l==pow_l1_new+2*pow_loq_new+pow_koloz_new);

	for (int i=0; i<pow_k1_new; ++i) {
		idxv[i]=varidx(1,Dim);
	}
	for (int i=0; i<pow_koq_new; ++i) {
		symbol s_mu, s_nu;
		varidx mu(s_mu,Dim-1);
		varidx nu(s_nu,Dim-1);
		idxv[pow_k1_new+2*i]=mu;
		idxv[pow_k1_new+2*i+1]=nu;
		coeff *= -lorentz_g(mu.toggle_variance(),nu.toggle_variance());
	}
	for (int i=0; i<pow_l1_new; ++i) {
		idxv[pow_k+i]=varidx(1,Dim);
	}
	for (int i=0; i<pow_loq_new; ++i) {
		symbol s_mu, s_nu;
		varidx mu(s_mu,Dim-1);
		varidx nu(s_nu,Dim-1);
		idxv[pow_k+pow_l1_new+2*i]=mu;
		idxv[pow_k+pow_l1_new+2*i+1]=nu;
		coeff *= -lorentz_g(mu.toggle_variance(),nu.toggle_variance());
	}
	for (int i=0; i<pow_koloz_new; ++i) {
		symbol s_mu, s_nu;
		varidx mu(s_mu,Dim-1);
		varidx nu(s_nu,Dim-1);
		idxv[pow_k1_new+2*pow_koq_new+i]=mu;
		idxv[pow_k+pow_l1_new+2*pow_loq_new+i]=nu;
		coeff *= -lorentz_g(mu.toggle_variance(),nu.toggle_variance());
	}
	ex t=get_tensor(pow_k,pow_l,idxv).subs(lst(kq==-koq_new,
											   lq==-loq_new,
											   kdotl==-koloz_new,
											   Dim == Dim-1));
	// cout << "t=" << t << endl;
	ex t2=simplify_indexed(coeff*t,sp);
	// cout << "t2=" << t2 << endl;
	return t2;
}       
	
po_redux_powers_vector po_redux::collect_powers(const ex & res,
												const symbol & kq_sym,
												const symbol & lq_sym,
												const symbol & kl_sym)
{
	po_redux_powers_vector v;
	if (is_exactly_a<add>(res)) {
		for (int i=0; i<res.nops(); ++i) {
			collect_powers_single_term(v,res.op(i), kq_sym, lq_sym, kl_sym);
		}
	} else {
		collect_powers_single_term(v,res, kq_sym, lq_sym, kl_sym);
	}

	for (po_redux_powers_vector::iterator it=v.begin(); it!=v.end(); ++it) {
		it->coeff=(it->coeff).normal();
	}

	return v;
}

void po_redux::collect_powers_single_term(po_redux_powers_vector & v,
										  const ex & res,
										  const symbol & kq_sym,
										  const symbol & lq_sym,
										  const symbol & kl_sym)
{
	ASSERT(!is_exactly_a<add>(res));
	int pow_kq=res.degree(kq_sym);
	int pow_lq=res.degree(lq_sym);
	int pow_kl=res.degree(kl_sym);
	ex coeff=res/(power(kq_sym,pow_kq)*power(lq_sym,pow_lq)*power(kl_sym,pow_kl));
	coeff=coeff.normal();

	ASSERT(!coeff.has(kq_sym));
	ASSERT(!coeff.has(lq_sym));
	ASSERT(!coeff.has(kl_sym));

	// search this power combination in v, add coeff if found, append if not
	bool found=false;
	po_redux_powers_vector::iterator it=v.begin();
	while ((!found)&&(it!=v.end())) {
		if ((it->pow_kq==pow_kq)&&(it->pow_lq==pow_lq)&&(it->pow_kl==pow_kl)) {
			found=true;
			it->coeff += coeff;
		}
		++it;
	}
	if (!found) {
		po_redux_powers p;
		p.coeff=coeff;
		p.pow_kq=pow_kq;
		p.pow_lq=pow_lq;
		p.pow_kl=pow_kl;
		v.push_back(p);
	}
}
	
} // namespace xloops
