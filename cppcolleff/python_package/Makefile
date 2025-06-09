
#### COMPILATION OPTIONS ####
C = gcc
CXX = g++
SWIG = swig
MACHINE = -m64
OPT_FLAG = -O4
SOURCES_I = cppcolleff.i
WRAPPEROBJS = cppcolleff_wrap.o
# INTERFACEOBJS = interface.o
CPPCOLLEFFLIB = ../build/libcppcolleff.a
PACKAGE = cppcolleff
PREFIX ?=
PIP ?= pip
PYTHON ?= python3
ifeq ($(CONDA_PREFIX),)
	PREFIX = sudo -H
	PIP = pip-sirius
	PYTHON = python-sirius
endif

PYTHON_INCLUDES := $(shell $(PYTHON)-config --includes)

# OBJECTS = $(WRAPPEROBJS) $(INTERFACEOBJS)
OBJECTS = $(WRAPPEROBJS)

SWIG_FLAGS= -c++ -python -threads
GEN_FLAG  = -pthread -std=c++17 -fPIC -fopenmp
DFLAGS = -DOPENMP
CPPFLAGS = $(MACHINE) $(OPT_FLAG) $(GEN_FLAG) $(DFLAGS)
LDFLAGS = $(MACHINE) $(GEN_FLAG) $(DFLAGS)

LIBS = $(shell gsl-config --libs)
LIBS += -L../build -lcppcolleff -lblas -lfftw3 -lfftw3_omp
INC = $(shell gsl-config --cflags)
INC += -I../include $(PYTHON_INCLUDES)

# INTERFACESRCS = $(INTERFACEOBJS:.o=.cpp)
# INTERFACEHDRS = $(INTERFACEOBJS:.o=.h)

.NOTPARALLEL:

all: cppcolleff/_cppcolleff.so

# $(shell $(CXX) -MM $(INTERFACESRCS) > .depend)
# -include .depend

cppcolleff/_cppcolleff.so: $(CPPCOLLEFFLIB) $(OBJECTS)
	$(CXX) $(LDFLAGS) -shared -Wl,--whole-archive -Wl,--no-as-needed $(CPPCOLLEFFLIB) -Wl,--no-whole-archive $(OBJECTS) $(LIBS) -o cppcolleff/_cppcolleff.so

# cppcolleff_wrap.cxx: $(SOURCES_I) $(INTERFACESRCS) $(INTERFACEHDRS)
cppcolleff_wrap.cxx: $(SOURCES_I)
	$(SWIG) $(SWIG_FLAGS) $(INC) cppcolleff.i && cp -f cppcolleff.py cppcolleff

install: uninstall all
	$(PREFIX) $(PIP) install --no-deps ./

uninstall:
	$(PREFIX) $(PIP) uninstall -y $(PACKAGE)

develop-install: develop-uninstall all
	$(PIP) install --no-deps -e ./

develop-uninstall:
	$(PIP) uninstall -y $(PACKAGE)

clean:
	-rm -rf *.o cppcolleff/*.so *_wrap.cxx cppcolleff/__pycache__ cppcolleff/cppcolleff.py cppcolleff.py build dist cppcolleff.egg-info

# $(INTERFACESRCS): $(INTERFACEHDRS)

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(INC) $< -o $@

%.o: %.cxx
	$(CXX) -c $(CPPFLAGS) $(INC) $< -o $@
