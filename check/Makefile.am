## Process this file with automake to produce Makefile.in

TESTS = run_checks
check_PROGRAMS = check_loop
check_loop_SOURCES = one_loop_two_point.cpp one_loop_three_point.cpp \
  po_redux_test.cpp fcntimer.cpp main.cpp check.h
check_loop_LDADD = ../lib/libxloops.a
INCLUDES = -I$(srcdir)/../lib
CLEANFILES = result.out
EXTRA_DIST = result.ref run_checks
