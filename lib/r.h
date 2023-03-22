/** @file r.h
 *
 *  Definition of necessary R functions. */

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

#ifndef __XLOOPS_R_H__
#define __XLOOPS_R_H__

#include <ginac/ginac.h>

namespace xloops {


/*
 *  Two-dimensional R function R2(t, b1, b2, z1, z2)
 *
 *  t:    index
 *  b1,2: parameters
 *  z1,2: arguments
 */

DECLARE_FUNCTION_6P(R2)
DECLARE_FUNCTION_2P(R2ex1)
DECLARE_FUNCTION_2P(R2ex2)


/*
 *  Three-dimensional R function R3(t, b1, b2, b3, z1, z2, z3)
 *
 *  t:     index
 *  b1..3: parameters
 *  z1..3: arguments
 */

DECLARE_FUNCTION_8P(R3)
DECLARE_FUNCTION_3P(R3ex2)
DECLARE_FUNCTION_3P(R3ex3)


/** Substitute any RXexY() functions by the real expressions. */
extern GiNaC::ex subs_Rex(const GiNaC::ex &e);


} // namespace xloops

#endif // __XLOOPS_R_H__
