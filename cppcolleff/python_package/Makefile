
#### COMPILATION OPTIONS ####
C         = gcc
CXX       = g++
SWIG      = swig3.0
PYTHON    = python3.4

MACHINE   = -m64
OPT_FLAG  = -O4
SWIG_FLAGS= -c++ -python -threads
GEN_FLAG  = -fopenmp -pthread -std=c++11 -fPIC
DFLAGS    = -DOPENMP

SOURCES_I = cppcolleff.i

WRAPPEROBJS = cppcolleff_wrap.o
OBJECTS = $(WRAPPEROBJS)


COMP_FLAGS = $(MACHINE) $(OPT_FLAG) $(GEN_FLAG) $(DFLAGS)
LDFLAGS    = $(MACHINE) $(GEN_FLAG) $(DFLAGS)
LIBS = -lgsl -lblas -L../build -lcppcolleff
CPPCOLLEFFLIB = ../build/libcppcolleff.a
INC = -I../include -I/usr/include/$(PYTHON)

SETUPARGS =

all: cppcolleff/_cppcolleff.so

cppcolleff/_cppcolleff.so: $(CPPCOLLEFFLIB) $(OBJECTS)
	$(CXX) $(LDFLAGS) -shared -Wl,--whole-archive -Wl,--no-as-needed $(CPPCOLLEFFLIB) -Wl,--no-whole-archive $(OBJECTS) $(LIBS) -o cppcolleff/_cppcolleff.so

cppcolleff_wrap.cxx: $(SOURCES_I)
	$(SWIG) $(SWIG_FLAGS) $(INC) cppcolleff.i && cp -f cppcolleff.py cppcolleff

install: all
	$(PYTHON) setup.py install $(SETUPARGS)

develop: all
	$(PYTHON) setup.py develop $(SETUPARGS)

clean:
	-rm -rf *.o cppcolleff/*.so *_wrap.cxx cppcolleff/__pycache__ cppcolleff/cppcolleff.py build dist cppcolleff.egg-info

%.o: %.cpp
	$(CXX) -c $(COMP_FLAGS) $(INC) $< -o $@

%.o: %.cxx
	$(CXX) -c $(COMP_FLAGS) $(INC) $< -o $@
