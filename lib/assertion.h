/** @file assertion.h
 *
 *  Assertion macro definition. */

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

#ifndef __XLOOPS_ASSERTION_H__
#define __XLOOPS_ASSERTION_H__

#include <assert.h>

#if !defined(ASSERT)
#if defined(DO_ASSERT)
#define ASSERT(X) assert(X)
#else
#define ASSERT(X) ((void)0)
#endif
#endif

#endif // ndef __XLOOPS_ASSERTION_H__
