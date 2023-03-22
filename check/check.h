/** @file check.h
 *
 *  Prototypes for all individual checks.
 *
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

#ifndef CHECK_H
#define CHECK_H

// fcntimer is defined in timer.cpp and used for timing check functions only:
unsigned fcntimer(unsigned fcn());

// prototypes for all individual checks must be unsigned fcn() in order to be
// able to use fcntimer() as a wrapper:
unsigned one_loop_two_point();
unsigned one_loop_three_point();
unsigned po_redux_test(void);

#endif // ndef CHECK_H
