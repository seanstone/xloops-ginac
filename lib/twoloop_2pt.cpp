/** @file Twoloop_2pt.cpp
 *
 *  Implementation of functions for Two-loop two-point Feynman integrals. */

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
#include "twoloop.h"
#include "r.h"
#include "utils.h"

// using namespace GiNaC;

namespace xloops {

static ex TwoLoop2Pt_eval(const ex &p0, const ex &p1,const ex &r0, const ex &r1, const ex &s,
						  const ex &q,
						  const ex &m ,
						  const ex &t,
						  const ex &topo,
						  const ex &rho)
{
	// The 2-point function vanishes for odd p1, r1 and for odd p0, r0 when q is zero
	if (p1.info(info_flags::odd) || r1.info(info_flags::odd) ){
		return 0;
	}
	else if (q.is_zero() && (p0.info(info_flags::odd) || r0.info(info_flags::odd)  )){
		return 0;
	}
	else
		return TwoLoop2Pt(p0, p1, r0, r1, s, q, m.op(0), t.op(0), topo, rho).hold();
}
	
static ex TwoLoop2Pt_series(const ex &p0, const ex &p1, const ex &r0, const ex &r1, const ex &s,
							const ex &q,
							const ex &m,
							const ex &t,
							const ex &topo,
                            const ex &rho,
							const relational &r, int order, unsigned options)
{
	// Check parameters for validity

	int tz = ex_to<numeric>(s).to_int();
	int p1_half =  ex_to<numeric>(p1/2).to_int();
	int p1_int = ex_to<numeric>(p1).to_int();
	int p0_int = ex_to<numeric>(p0).to_int();
	int r1_int = ex_to<numeric>(r1).to_int();
	int r0_int = ex_to<numeric>(r0).to_int();


	
	if (!p0.info(info_flags::nonnegint) || !p1.info(info_flags::nonnegint)
		|| !r0.info(info_flags::nonnegint) || !r1.info(info_flags::nonnegint)
		|| !s.info(info_flags::nonnegint)
		|| !t.op(0).info(info_flags::posint) || !t.op(1).info(info_flags::posint)
		|| !t.op(2).info(info_flags::posint) || !t.op(3).info(info_flags::posint)
		|| !t.op(4).info(info_flags::posint)
		)
		throw(std::logic_error("TwoLoop2Pt_series: p0<0 or p1<0 or r0<0 or r1<0 or t[i]<=0"));
	
	
	// Check whether we can do the expansion and get parameters
	
	if (!r.rhs().is_zero())
		throw(std::logic_error("OneLoop2Pt_series: don't know the expansion at eps!=0"));
	const ex eps = r.lhs();
	
	// First we try to kill the internal momenta mixing term in the propagator.
	// We have 8 diferent topologies where topolgy 5-6-7-8 are factorizable.

	
	int t1 = ex_to<numeric>(topo).to_int();
	ex c;
	switch(t1)
		{
		case 1:
			
			break;
		case 2:
			
			break;
		case 3:
			
			break;
		case 4:
			
			break;
			
		case 5:
			if ( !s.info(info_flags::odd)){
				for(int i =0; i<= p1_half; i++){
					c +=  binomial(i, numeric(p1_int,2))
						*OneLoop3Pt(p0, 2*i, p1-2*i, q, 0, 0, m.op(0), m.op(1), m.op(2),
									t.op(0), t.op(1), t.op(2),rho);
				}
				c *= OneLoop1Pt(r0+r1, m.op(3), t.op(3), rho)*CoeffDim(r0_int, r1_int, -1, eps)
					*numeric(2, 1+tz);
			} else{
				c = 0;
			}
			break;
		case 6:
			if ( !s.info(info_flags::odd)){
				c = OneLoop2Pt(p0, p1, q, m.op(0), m.op(1),t.op(0), t.op(1),rho)
					*OneLoop2Pt(r0, r1, -q, m.op(2), m.op(3),t.op(2), t.op(3),rho)*numeric(2,1+tz);
			}else{
				c = 0;
			}
			break;
		case 7:
			if ( !s.info(info_flags::odd)){
				c = OneLoop2Pt(p0, p1, q, m.op(0), m.op(1), t.op(0), t.op(1),rho)
					*OneLoop1Pt(r0+r1, m.op(2), t.op(2), rho)*CoeffDim(r0_int, r1_int, -1, eps)
					*numeric(2,1+tz);
			} else{
				c = 0;
			}
			break;
		case 8:
			if ( !s.info(info_flags::odd)){
				if (!q.is_zero()){
					for (int i =0; i <= p0_int; i++){
						c += pow(-q,i)*binomial(i, p0)*CoeffDim(i, p1_int, -1, eps)*
							OneLoop1Pt(i+p1,m.op(0),t.op(0),rho);

					}
					 c *= CoeffDim(r0_int, r1_int, -1, eps)*OneLoop1Pt(r0+r1,m.op(1),t.op(1),rho)
						 *numeric(2,1+tz);
				} else{
					c = 0;
				}
				
			} else{
				c = 0;
			}
			
			break;
			
		default:
			break;			  			  
		}
	
	
	
	return c.series(r, order);

	/*
	ex c = reduce_OneLoop2Pt(
	ex_to<numeric>(p0).to_int(), ex_to<numeric>(p1).to_int(),
	q,
	m1, m2,
	ex_to<numeric>(t1).to_int(), ex_to<numeric>(t2).to_int(),
	rho,
	eps
	);
	*/
	// Expand into series in eps

}

REGISTER_FUNCTION(TwoLoop2Pt, eval_func(TwoLoop2Pt_eval).
				  series_func(TwoLoop2Pt_series));



} // namespace xloops
