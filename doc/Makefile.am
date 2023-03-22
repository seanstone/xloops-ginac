## Process this file with automake to produce Makefile.in

xloopsdocdir = $(datadir)/doc/xloops
xloopsdoc_DATA = @DOC_TARGETS@

xloops.dvi: xloops.tex
	TEXINPUTS=$(srcdir):$$TEXINPUTS \
	latex $<
	TEXINPUTS=$(srcdir):$$TEXINPUTS \
	latex $<

xloops.ps: xloops.dvi
	TEXPICTS=$(srcdir):$$TEXPICTS \
	dvips -o $@ $<

CLEANFILES = *.ps *.dvi *.aux *.log
