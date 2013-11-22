# Flux makefile
#
# Author : Lester. O. Hedges
# Email  : lester.hedges+flux@gmail.com
# Date   : January 21st 2013

# Load the common configuration file
include ../config.mk

all: flux

flux: lodepng.cpp flux.cpp
	$(CXX) $(CXXFLAGS) lodepng.cpp flux.cpp -o flux

clean:
	rm -f flux

.PHONY: all
