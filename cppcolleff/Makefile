## CPPCOLLEFF
## ========
## Author:      Accelerator Physics Group - LNLS
## contact:     xresende@gmail.com
## affiliation: Laboratorio Nacional de Luz Sincrotron
##
## The MIT License (MIT)
##
## Copyright (c) <year> <copyright holders>
##
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in
## all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
## THE SOFTWARE.

#### READS LIB VERSION ####

FILE=VERSION
VERSION=$(shell cat ${FILE})

#### COMPILATION OPTIONS ####
CC          = gcc
CXX         = g++
AR          = ar
MACHINE     = -m64
OPT_FLAG    = -O3 -fopenmp -std=c++11 -fPIC #-DOPENMP
DBG_FLAG    = -O0 -g3 -std=c++11 -fPIC
ARFLAGS     = rcs
DFLAGS      = -DVERSION=$(VERSION)
LIBSOURCES_CPP = essentials.cpp \
				 Bunch.cpp \
				 Ring.cpp \
				 Wake.cpp \
				 Results.cpp \
				 Feedback.cpp \
				 cppcolleff.cpp
BINSOURCES_CPP = exec.cpp

AUXFILES  = VERSION

LIBS = -lgsl -lgslcblas -lm -fopenmp
INC  = -I./include
BINDEST_DIR = /usr/local/bin
LIBDEST_DIR = /usr/local/lib
INCDEST_DIR = /usr/local/include

OBJDIR = build
SRCDIR = src
INCDIR = include

PYTHON_PACKAGE_DIR = python_package

# this is so that compilation is forced
$(shell touch $(SRCDIR)/cppcolleff.cpp)

ifeq ($(MAKECMDGOALS),cppcolleff-debug)
  CFLAGS    = $(MACHINE) $(DBG_FLAG) $(DFLAGS)
else
  CFLAGS    = $(MACHINE) $(OPT_FLAG) $(DFLAGS)
endif

LIBOBJECTS  = $(addprefix $(OBJDIR)/, $(LIBSOURCES_CPP:.cpp=.o))
BINOBJECTS  = $(addprefix $(OBJDIR)/, $(BINSOURCES_CPP:.cpp=.o))
LDFLAGS    = $(MACHINE)

#### DERIVED CONDITIONALS AND VARIABLES ####
ifeq ($(shell hostname), uv100)
    CC          = /progs/users/gcc-4.7.3/bin/gcc
    CXX         = /progs/users/gcc-4.7.3/bin/g++
    LIBS        = -Wl,-rpath,/progs/users/gcc-4.7.3/lib64 ~/bin/GSL/lib/libgsl.a ~/bin/GSL/lib/libgslcblas.a
    INC         = -I../bin/GSL/include/
endif


.PHONY: all alllibs cppcolleff clean cleanall

#### TARGETS ####

all:  libcppcolleff cppcolleff python_package

#### GENERATES DEPENDENCY FILE ####
$(shell $(CXX) -MM $(CFLAGS) $(addprefix $(SRCDIR)/, $(LIBSOURCES_CPP)) $(addprefix $(SRCDIR)/, $(BINSOURCES_CPP)) | sed 's/.*\.o/$(OBJDIR)\/&/' > .depend)
-include .depend

libcppcolleff: $(OBJDIR)/libcppcolleff.a

cppcolleff: $(OBJDIR)/cppcolleff

python_package: $(PYTHON_PACKAGE_DIR)/cppcolleff/_cppcolleff.so

$(PYTHON_PACKAGE_DIR)/cppcolleff/_cppcolleff.so: libcppcolleff
	$(MAKE) -C $(PYTHON_PACKAGE_DIR)

$(OBJDIR)/libcppcolleff.a: $(LIBOBJECTS)
	$(AR) $(ARFLAGS) $@ $^

$(OBJDIR)/cppcolleff: libcppcolleff $(BINOBJECTS)
	$(CXX) $(LDFLAGS) $(BINOBJECTS) $(OBJDIR)/libcppcolleff.a $(LIBS) -o $@

$(LIBOBJECTS): | $(OBJDIR)

$(BINOBJECTS): | $(OBJDIR)

$(OBJDIR):
	mkdir $(OBJDIR)

install: uninstall all
	cp $(OBJDIR)/cppcolleff $(BINDEST_DIR)
	cp $(OBJDIR)/libcppcolleff.a $(LIBDEST_DIR)
	cp -r $(INCDIR)/cppcolleff $(INCDEST_DIR)
	$(MAKE) install -C $(PYTHON_PACKAGE_DIR)

develop: uninstall all
	ln -srf $(OBJDIR)/cppcolleff $(BINDEST_DIR)
	ln -srf $(OBJDIR)/libcppcolleff.a $(LIBDEST_DIR)
	ln -srf $(INCDIR)/cppcolleff $(INCDEST_DIR)
	$(MAKE) develop -C $(PYTHON_PACKAGE_DIR)

$(BINDEST_DIR):
	mkdir $(BINDEST_DIR)

$(LIBDEST_DIR):
	mkdir $(LIBDEST_DIR)

$(INCDEST_DIR):
	mkdir $(INCDEST_DIR)

clean:
	-rm -rf $(OBJDIR) cppcolleff cppcolleff-debug .depend *.out *.dat *~ *.o *.a
	$(MAKE) clean -C $(PYTHON_PACKAGE_DIR)

uninstall:
	-rm -rf $(BINDEST_DIR)/cppcolleff
	-rm -rf $(LIBDEST_DIR)/libcppcolleff.a
	-rm -rf $(INCDEST_DIR)/cppcolleff


#### RULES ####

*.cpp: VERSION
	touch $(SRCDIR)/*.cpp
*.cc: VERSION
	touch $(SRCDIR)/*.cc
*.c: VERSION
	touch $(SRCDIR)/*.c

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) -c $(CFLAGS) $(INC) -I./$(SRCDIR) $< -o $@;
