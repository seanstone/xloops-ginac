## Process this file with automake to produce Makefile.in

bin_PROGRAMS = ginsh-xloops
ginsh_xloops_SOURCES = ginsh_parser.yy ginsh_lexer.ll ginsh.h ginsh_extensions.h
ginsh_xloops_LDADD = @GINSH_LIBS@ ../lib/libxloops.a

INCLUDES = -I$(srcdir)/../lib
YFLAGS = -d

#man_MANS = ginsh.1

CLEANFILES = ginsh_fcn_help.h ginsh_op_help.h
EXTRA_DIST = ginsh_parser.h ginsh_fcn_help.sed ginsh_op_help.sed ginsh.1

# files created by sed scripts
ginsh_fcn_help.h: ginsh.1 $(srcdir)/ginsh_fcn_help.sed
	sed -n -f $(srcdir)/ginsh_fcn_help.sed <$< >$@

ginsh_op_help.h: ginsh.1 $(srcdir)/ginsh_op_help.sed
	sed -n -f $(srcdir)/ginsh_op_help.sed <$< >$@

ginsh_parser.o: ginsh_fcn_help.h ginsh_op_help.h
