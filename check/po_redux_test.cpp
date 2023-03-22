/** @file po_redux_test.cpp
 *
 *  Check po_redux functions. */

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
#include "po_redux.h"

using namespace std;
using namespace GiNaC;
using namespace xloops;

static unsigned check(const ex & s, const ex & r)
{
	if (normal(s-r)!=0) {
		clog << "transformed numerator is " << s.normal() 
			 << " (should be " << r << ")" << endl;
		return 1;
	}
	// clog << s << "==" << r << ", ok" << endl;
	return 0;
}

static ex transform_1_0(po_redux & redux, int pk0, int pkoq, int pl0, int ploq,
						int pkoloz, const ex & kq, const ex & lq,
						const ex & kl, const ex & dim)
{
	po_redux_powers_vector v;
	ex res;
		
	v=redux.transform_1_to_0(pk0,pkoq,pl0,ploq,pkoloz,dim);
	for (po_redux_powers_vector::const_iterator cit=v.begin(); cit!=v.end(); ++cit) {
		res += cit->coeff*power(kq,cit->pow_kq)*power(lq,cit->pow_lq)*power(kl,cit->pow_kl);
	}

	return res;
}

static unsigned test_1_0(void)
{
	unsigned result=0;

	po_redux redux;
	symbol kq("(k^2)");
	symbol lq("(l^2)");
	symbol kdotl("k.l");
	symbol dim("dim");

	ex s;

	// check for tensors 0th degree (scalar case)
	
	s  = transform_1_0(redux,0,0,0,0,0,kq,lq,kdotl,dim);
	result += check(s,1);

	// check for tensors 2nd degree
	
	// k^2 = k_0^2 - k_\perp^2
	s  = transform_1_0(redux,2,0,0,0,0,kq,lq,kdotl,dim);
	s -= transform_1_0(redux,0,1,0,0,0,kq,lq,kdotl,dim);
	result += check(s,kq);

	// l^2 = l_0^2 - l_\perp^2
	s  = transform_1_0(redux,0,0,2,0,0,kq,lq,kdotl,dim);
	s -= transform_1_0(redux,0,0,0,1,0,kq,lq,kdotl,dim);
	result += check(s,lq);

	// (k.l) = k_0 l_0 - k_\perp l_\perp z
	s  = transform_1_0(redux,1,0,1,0,0,kq,lq,kdotl,dim);
	s -= transform_1_0(redux,0,0,0,0,1,kq,lq,kdotl,dim);
	result += check(s,kdotl);

	// check for tensors 4th degree

	// (k^2)^2 = k_0^4 + (k_\perp^2)^2 - 2 k_0^2 k_\perp^2
	s  = transform_1_0(redux,4,0,0,0,0,kq,lq,kdotl,dim);
	s += transform_1_0(redux,0,2,0,0,0,kq,lq,kdotl,dim);
	s -= 2*transform_1_0(redux,2,1,0,0,0,kq,lq,kdotl,dim);
	result += check(s,pow(kq,2));

	// (l^2)^2 = l_0^4 + (l_\perp^2)^2 - 2 l_0^2 l_\perp^2
	s  = transform_1_0(redux,0,0,4,0,0,kq,lq,kdotl,dim);
	s += transform_1_0(redux,0,0,0,2,0,kq,lq,kdotl,dim);
	s -= 2*transform_1_0(redux,0,0,2,1,0,kq,lq,kdotl,dim);
	result += check(s,pow(lq,2));

	// (k.l)^2 = k_0^2 l_0^2 + (k_\perp l_\perp z)^2 - 2 k_0 l_0 (k_\perp l_\perp z)
	s  = transform_1_0(redux,2,0,2,0,0,kq,lq,kdotl,dim);
	s += transform_1_0(redux,0,0,0,0,2,kq,lq,kdotl,dim);
	s -= 2*transform_1_0(redux,1,0,1,0,1,kq,lq,kdotl,dim);
	result += check(s,pow(kdotl,2));

	// k^2 (k.l) = k_0^3 l_0 - k_0 k_\perp^2 l_0 - k_0^2 (k_\perp l_\perp z) + k_\perp^2 (k_\perp l_\perp z)
	s  = transform_1_0(redux,3,0,1,0,0,kq,lq,kdotl,dim);
	s -= transform_1_0(redux,1,1,1,0,0,kq,lq,kdotl,dim);
	s -= transform_1_0(redux,2,0,0,0,1,kq,lq,kdotl,dim);
	s += transform_1_0(redux,0,1,0,0,1,kq,lq,kdotl,dim);
	result += check(s,kq*kdotl);

	// check for tensor 6th degree

	// k^2 (k.l) l^2 = k_0^3 l_0^3 - k_0^2 l_0^2 (k_\perp l_\perp z) - k_0^3 l_0 l_\perp^2 + k_0^2 l_\perp^2 (k_\perp l_\perp z) - k_0 k_\perp^2 l_0^3 + k_\perp^2 l_0^2 (k_\perp l_\perp z) + k_0 k_\perp^2 l_0 l_\perp^2 - k_\perp^2 l_\perp^2 (k_\perp l_\perp z)
	s  = transform_1_0(redux,3,0,3,0,0,kq,lq,kdotl,dim);
	s -= transform_1_0(redux,2,0,2,0,1,kq,lq,kdotl,dim);
	s -= transform_1_0(redux,3,0,1,1,0,kq,lq,kdotl,dim);
	s += transform_1_0(redux,2,0,0,1,1,kq,lq,kdotl,dim);
	s -= transform_1_0(redux,1,1,3,0,0,kq,lq,kdotl,dim);
	s += transform_1_0(redux,0,1,2,0,1,kq,lq,kdotl,dim);
	s += transform_1_0(redux,1,1,1,1,0,kq,lq,kdotl,dim);
	s -= transform_1_0(redux,0,1,0,1,1,kq,lq,kdotl,dim);
	result += check(s,kq*kdotl*lq);

	return result;
}

static ex transform_2_0(po_redux & redux, int pk0, int pk1, int pkoq, int pl0,
						int pl1, int ploq, int pkoloz, const ex & kq,
						const ex & lq, const ex & kl, const ex & dim)
{
	po_redux_powers_vector v;
	ex res;
		
	v=redux.transform_2_to_0(pk0,pk1,pkoq,pl0,pl1,ploq,pkoloz,dim);
	for (po_redux_powers_vector::const_iterator cit=v.begin(); cit!=v.end(); ++cit) {
		res += cit->coeff*power(kq,cit->pow_kq)*power(lq,cit->pow_lq)*power(kl,cit->pow_kl);
	}

	return res;
}

static unsigned test_2_0(void)
{
	unsigned result=0;

	po_redux redux;
	symbol kq("(k^2)");
	symbol lq("(l^2)");
	symbol kdotl("k.l");
	symbol dim("dim");

	ex s;

	// check for tensors 0th degree (scalar case)
	
	s  = transform_2_0(redux,0,0,0,0,0,0,0,kq,lq,kdotl,dim);
	result += check(s,1);

	// check for tensors 2nd degree
	
	// k^2 = k_0^2 - k_1^2 - k_\perp^2
	s  = transform_2_0(redux,2,0,0,0,0,0,0,kq,lq,kdotl,dim);
	s -= transform_2_0(redux,0,2,0,0,0,0,0,kq,lq,kdotl,dim);
	s -= transform_2_0(redux,0,0,1,0,0,0,0,kq,lq,kdotl,dim);
	result += check(s,kq);

	// l^2 = l_0^2 - l_1^2 - l_\perp^2
	s  = transform_2_0(redux,0,0,0,2,0,0,0,kq,lq,kdotl,dim);
	s -= transform_2_0(redux,0,0,0,0,2,0,0,kq,lq,kdotl,dim);
	s -= transform_2_0(redux,0,0,0,0,0,1,0,kq,lq,kdotl,dim);
	result += check(s,lq);

	// (k.l) = k_0 l_0 - k_1 l_1 - k_\perp l_\perp z
	s  = transform_2_0(redux,1,0,0,1,0,0,0,kq,lq,kdotl,dim);
	s -= transform_2_0(redux,0,1,0,0,1,0,0,kq,lq,kdotl,dim);
	s -= transform_2_0(redux,0,0,0,0,0,0,1,kq,lq,kdotl,dim);
	result += check(s,kdotl);

	// k^2 (k.l) = k_0^3 l_0 - k_0^2 k_1 l_1 - k_0^2 (k_\perp l_\perp z) + k_0 k_1^2 l0 + k_1^3 l_1 + k_1^2 (k_\perp l_\perp z) - k_0 k_\perp^2 l_0 + k_1 k_\perp^2 l_1 + k_\perp^2 (k_\perp l_\perp z)
	s  = transform_2_0(redux,3,0,0,1,0,0,0,kq,lq,kdotl,dim);
	s -= transform_2_0(redux,2,1,0,0,1,0,0,kq,lq,kdotl,dim);
	s -= transform_2_0(redux,2,0,0,0,0,0,1,kq,lq,kdotl,dim);
	s -= transform_2_0(redux,1,2,0,1,0,0,0,kq,lq,kdotl,dim);
	s += transform_2_0(redux,0,3,0,0,1,0,0,kq,lq,kdotl,dim);
	s += transform_2_0(redux,0,2,0,0,0,0,1,kq,lq,kdotl,dim);
	s -= transform_2_0(redux,1,0,1,1,0,0,0,kq,lq,kdotl,dim);
	s += transform_2_0(redux,0,1,1,0,1,0,0,kq,lq,kdotl,dim);
	s += transform_2_0(redux,0,0,1,0,0,0,1,kq,lq,kdotl,dim);
	result += check(s,kq*kdotl);

	return result;
}

static ex transform_2_1(po_redux & redux, int pk1, int pkoq, int pl1, int ploq,
						int pkoloz, const ex & koq, const ex & loq,
						const ex & koloz, const ex & dimo)
{
	po_redux_powers_vector v;
	ex res;
		
	v=redux.transform_n_plus_1_to_n(pk1,pkoq,pl1,ploq,pkoloz,dimo,1);
	for (po_redux_powers_vector::const_iterator cit=v.begin(); cit!=v.end(); ++cit) {
		res += cit->coeff*power(koq,cit->pow_kq)*power(loq,cit->pow_lq)*power(koloz,cit->pow_kl);
	}

	return res;
}

static unsigned test_2_1(void)
{
	unsigned result=0;

	po_redux redux;
	symbol koq("(k~_\\perp^2)");
	symbol loq("(l~_\\perp^2)");
	symbol koloz("(k~_\\perp l~_\\perp z~)");
	symbol dimo("dimo");

	ex s;

	// check for tensors 0th degree (scalar case)
	
	s  = transform_2_1(redux,0,0,0,0,0,koq,loq,koloz,dimo);
	result += check(s,1);

	// check for tensors 2nd degree
	
	// k~_\perp^2 = k_1^2 + k_\perp^2
	s  = transform_2_1(redux,2,0,0,0,0,koq,loq,koloz,dimo);
	s += transform_2_1(redux,0,1,0,0,0,koq,loq,koloz,dimo);
	result += check(s,koq);

	// l~_\perp^2 = l_1^2 + l_\perp^2
	s  = transform_2_1(redux,0,0,2,0,0,koq,loq,koloz,dimo);
	s += transform_2_1(redux,0,0,0,1,0,koq,loq,koloz,dimo);
	result += check(s,loq);

	// k~_\perp l~_\perp z~ = k_1 l_1 + k_\perp l_\perp z
	s  = transform_2_1(redux,1,0,1,0,0,koq,loq,koloz,dimo);
	s += transform_2_1(redux,0,0,0,0,1,koq,loq,koloz,dimo);
	result += check(s,koloz);

	// check for tensors 4th degree

	// (k~_\perp^2) = k_1^4 + (k_\perp^2)^2 + 2 k_1^2 k_\perp^2
	s  = transform_2_1(redux,4,0,0,0,0,koq,loq,koloz,dimo);
	s += transform_2_1(redux,0,2,0,0,0,koq,loq,koloz,dimo);
	s += 2*transform_2_1(redux,2,1,0,0,0,koq,loq,koloz,dimo);
	result += check(s,pow(koq,2));

	// (l~_\perp^2) = l_1^4 + (l_\perp^2)^2 + 2 l_1^2 l_\perp^2
	s  = transform_2_1(redux,0,0,4,0,0,koq,loq,koloz,dimo);
	s += transform_2_1(redux,0,0,0,2,0,koq,loq,koloz,dimo);
	s += 2*transform_2_1(redux,0,0,2,1,0,koq,loq,koloz,dimo);
	result += check(s,pow(loq,2));

	// (k~_\perp l~_\perp z~)^2 = k_1^2 l_1^2 + (k_\perp l_\perp z)^2 + 2 k_1 l_1 (k_\perp l_\perp z)
	s  = transform_2_1(redux,2,0,2,0,0,koq,loq,koloz,dimo);
	s += transform_2_1(redux,0,0,0,0,2,koq,loq,koloz,dimo);
	s += 2*transform_2_1(redux,1,0,1,0,1,koq,loq,koloz,dimo);
	result += check(s,pow(koloz,2));

	// k~_\perp^2 (k~_\perp l~_\perp z~) = k_1^3 l_1 + k_1 k_\perp^2 l_1 + k_1^2 (k_\perp l_\perp z) + k_\perp^2 (k_\perp l_\perp z)
	s  = transform_2_1(redux,3,0,1,0,0,koq,loq,koloz,dimo);
	s += transform_2_1(redux,1,1,1,0,0,koq,loq,koloz,dimo);
	s += transform_2_1(redux,2,0,0,0,1,koq,loq,koloz,dimo);
	s += transform_2_1(redux,0,1,0,0,1,koq,loq,koloz,dimo);
	result += check(s,koq*koloz);

	// k~_\perp^2 (k~_\perp l~_\perp z~) l~_\perp^2 = k_1^3 l_1^3 + k_1^2 l_1^2 (k_\perp l_\perp z) + k_1^3 l_1 l_\perp^2 + k_1^2 l_\perp^2 (k_\perp l_\perp z) + k_1 k_\perp^2 l_1^3 + k_\perp^2 l_1^2 (k_\perp l_\perp z) + k_1 k_\perp^2 l_1 l_\perp^2 + k_\perp^2 l_\perp^2 (k_\perp l_\perp z)
	s  = transform_2_1(redux,3,0,3,0,0,koq,loq,koloz,dimo);
	s += transform_2_1(redux,2,0,2,0,1,koq,loq,koloz,dimo);
	s += transform_2_1(redux,3,0,1,1,0,koq,loq,koloz,dimo);
	s += transform_2_1(redux,2,0,0,1,1,koq,loq,koloz,dimo);
	s += transform_2_1(redux,1,1,3,0,0,koq,loq,koloz,dimo);
	s += transform_2_1(redux,0,1,2,0,1,koq,loq,koloz,dimo);
	s += transform_2_1(redux,1,1,1,1,0,koq,loq,koloz,dimo);
	s += transform_2_1(redux,0,1,0,1,1,koq,loq,koloz,dimo);
	result += check(s,koq*koloz*loq);

	return result;
}

unsigned po_redux_test(void)
{
	unsigned result=0;

	cout << "checking parallel space reduction" << flush;
	clog << "---------parallel space reduction:" << endl;

	result += test_1_0();  cout << '.' << flush;
	result += test_2_0();  cout << '.' << flush;
	result += test_2_1();  cout << '.' << flush;

	if (!result) {
		cout << " passed ";
		clog << "(no output)" << endl;
	} else {
		cout << " failed ";
	}

	return result;
}
