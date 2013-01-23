# Flux makefile
#
# Author : Lester. O. Hedges
# Email  : lester.hedges+flux@gmail.com
# Date   : January 21st 2013

# Tell make that these are phony targets
.PHONY: all clean install uninstall

include config.mk

# Build all of the executable files
all:
	$(MAKE) -C src

# Clean up the executable files
clean:
	$(MAKE) -C src clean

# Install the executable, man page, and shared library
install:
	$(MAKE) -C src
	$(INSTALL) -d $(IFLAGS_EXEC) $(PREFIX)/bin
	$(INSTALL) -d $(IFLAGS_EXEC) $(PREFIX)/man
	$(INSTALL) -d $(IFLAGS_EXEC) $(PREFIX)/man/man1
	$(INSTALL) $(IFLAGS_EXEC) src/flux $(PREFIX)/bin
	$(INSTALL) $(IFLAGS) man/flux.1 $(PREFIX)/man/man1

# Uninstall the executable, man page, and shared library
uninstall:
	rm -f $(PREFIX)/bin/flux
	rm -f $(PREFIX)/man/man1/flux.1
