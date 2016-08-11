
#### COMPILATION OPTIONS ####
C         = gcc
CXX       = g++
SWIG      = swig3.0
MACHINE   = -m64
OPT_FLAG  = -O3
SOURCES_I =     cppcolleff.i
WRAPPEROBJS =   cppcolleff_wrap.o
CPPCOLLEFFLIB = ../build/libcppcolleff.a

OBJECTS = $(WRAPPEROBJS)

CPPFLAGS = -std=c++11 -fPIC $(OPT_FLAG)
LIBS = -lgsl -lblas -L../build -lcppcolleff
INC = -I/usr/include/python3.4 -I../include
PYTHON = python3
SETUPARGS =

.NOTPARALLEL:

all: cppcolleff/_cppcolleff.so

cppcolleff/_cppcolleff.so: $(CPPCOLLEFFLIB) $(OBJECTS)
	$(CXX) -shared -Wl,--whole-archive $(CPPCOLLEFFLIB) -Wl,--no-whole-archive $(OBJECTS) $(LIBS) -o cppcolleff/_cppcolleff.so

cppcolleff_wrap.cxx: $(SOURCES_I)
	$(SWIG) -c++ -python $(INC) cppcolleff.i && cp -f cppcolleff.py cppcolleff

install: all
	$(PYTHON) setup.py install $(SETUPARGS)

develop: all
	$(PYTHON) setup.py develop $(SETUPARGS)

clean:
	-rm -rf *.o cppcolleff/*.so *_wrap.cxx cppcolleff/__pycache__ cppcolleff/cppcolleff.py build dist cppcolleff.egg-info

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(INC) $< -o $@

%.o: %.cxx
	$(CXX) -c $(CPPFLAGS) $(INC) $< -o $@
