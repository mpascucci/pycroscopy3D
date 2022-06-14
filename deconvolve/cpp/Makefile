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
# Compiler and libraries
CXX := g++

CXX_EXTRA_OPTIONS += -std=c++11

LD_EXTRA_OPTIONS += -pthread 
LD	 = g++

ifdef openmp
 OPENMP_CPP = -fopenmp
 OPENMP_LD = -lfftw3_omp -fopenmp
 ifdef fftw_threads
  OPENMP_CPP += -DUSE_FFTW_THREADS
 endif
else
 OPENMP_CPP = 
 OPENMP_LD = 
endif



######################################################

LIBNAME := libdeconvolve.so

SRCSUBDIR := src
OBJSUBDIR := obj

INCLUDE  = -I$(SRCSUBDIR) -Iinclude
LIBRARIES += $(OPENMP_LD) -lfftw3 -lfftw3f 

CXXFLAGS := $(OPENMP_CPP) -Ofast -march=native -fPIC -g $(EXTRA_OPTIONS) $(CXX_EXTRA_OPTIONS) $(INCLUDE)

######################################################

OBJS	= $(patsubst $(SRCSUBDIR)/%.cpp,$(OBJSUBDIR)/%.o,$(wildcard $(SRCSUBDIR)/*.cpp))

######################################################

all: $(OBJSUBDIR) $(LIBNAME)

clean:
	rm -rf core* $(LIBNAME) $(OBJSUBDIR)

$(LIBNAME): $(OBJS)
	@echo
	@echo "--------- LINKING --- $@ "
	rm -f $(LIBNAME)
	$(LD) -shared -o $@ $^ $(LIBRARIES) $(LD_EXTRA_OPTIONS)
	@echo
	@echo '--------- Make done '
	@echo

$(OBJSUBDIR):
	@echo
	@echo "--------- Making dir: $@ "
	mkdir -p $(OBJSUBDIR)
	@echo

$(OBJSUBDIR)/%.o: $(SRCSUBDIR)/%.cpp 
	@echo
	@echo "------------ $< "
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $<
	@echo

.PHONY: help
help:
	doxygen Doxygen
