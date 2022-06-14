##
## When compiling with the optional switches, OpenMP can be enabled
## with either initialization or without of FFTW threads. When using
## switches, call make as
##
## make openmp=1 fftw_threads=1
##
## Optional switches:
##
## openmp: if set, openmp support is enabled
## fftw_threads: if set, FFTW fftw_init_threads and fftw_plan_with_nthreads
##


######################################################

OBJS	= $(patsubst $(SRCSUBDIR)/%.cpp,$(OBJSUBDIR)/%.o,$(wildcard $(SRCSUBDIR)/*.cpp))

######################################################

all: cppapi help

clean:
	make -C cpp clean

cppapi:
	+make -C cpp

.PHONY: help
help:
	doxygen Doxygen
