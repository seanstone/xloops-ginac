## Process this file with automake to produce Makefile.in
SUBDIRS = lib check xloops ginsh doc

# Requires automake 1.4
AUTOMAKE_OPTIONS = 1.4

# All the rest of the distributed files
EXTRA_DIST = autogen.sh depcomp

# Rule to build tar-gzipped distribution package
$(PACKAGE)-$(VERSION).tar.gz: dist

# Rule to build RPM distribution package
rpm: $(PACKAGE)-$(VERSION).tar.gz
	rpm -ta --clean $(PACKAGE)-$(VERSION).tar.gz
