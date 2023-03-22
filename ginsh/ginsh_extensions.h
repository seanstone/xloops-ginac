/** @file ginsh_extensions.h
 *
 *  The contents of this file are included in the ginsh parser. This makes
 *  it possible to create a customized version of ginsh that adds xloops-
 *  specific functions. */

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

#include "oneloop.h"
#include "r.h"

static ex f_tensor_reduction(const exprseq &e)
{
	return xloops::tensor_reduction(e[0], e[1]);
}

static ex f_subs_Rex(const exprseq &e)
{
	return xloops::subs_Rex(e[0]);
}

// Table of names and descriptors of functions to be added
static const fcn_init extended_fcns[] = {
	{"tensor_reduction", fcn_desc(f_tensor_reduction, 2)},
	{"subs_Rex", fcn_desc(f_subs_Rex, 1)},
	{NULL, fcn_desc(f_dummy, 0)} // End marker
};

// Table of help strings for functions
static const fcn_help_init extended_help[] = {
	{"tensor_reduction", "reduce OneLoopNPt() functions to ScalarNPt() functions"},
	{"subs_Rex", "substitute RXExY functions by their analytic equivalent"},
	{NULL, NULL} // End marker
};
