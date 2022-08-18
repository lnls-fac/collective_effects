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
GEN_FLAG    = -pthread -std=c++17 -fPIC -fopenmp
OPT_FLAG    = -O4 $(GEN_FLAG)
DBG_FLAG    = -O0 -g3 $(GEN_FLAG)
ARFLAGS     = rcs
DFLAGS      = -DVERSION=$(VERSION) -DOPENMP
LDFLAGS     = $(MACHINE) $(GEN_FLAG) $(DFLAGS)
WL_FLAGS    = -Wl,--no-as-needed

LIBSOURCES_CPP = essentials.cpp \
				 Bunch.cpp \
				 Ring.cpp \
				 Wake.cpp \
				 Results.cpp \
				 Feedback.cpp \
				 cppcolleff.cpp
BINSOURCES_CPP = exec.cpp

AUXFILES  = VERSION

LIBS = $(shell gsl-config --libs)
LIBS += -lfftw3 -lfftw3_omp
INC = $(shell gsl-config --cflags)
INC += -I./include

ifeq ($(CONDA_PREFIX),)
    PREFIX = /usr/local
else
    PREFIX = $(CONDA_PREFIX)
endif
BINDEST_DIR = $(PREFIX)/bin
LIBDEST_DIR = $(PREFIX)/lib
INCDEST_DIR = $(PREFIX)/include

OBJDIR = build
SRCDIR = src
INCDIR = include

PYTHON_PACKAGE_DIR = python_package

# this is so that compilation is forced
$(shell touch $(SRCDIR)/cppcolleff.cpp)

WARNINGS_CFLAGS += -Wall
# reasonable and standard:
WARNINGS_CFLAGS += -Wextra
# warn the user if a variable declaration shadows one from a parent context:
WARNINGS_CFLAGS += -Wshadow
# warn the user if a class with virtual functions has a non-virtual destructor:
WARNINGS_CFLAGS += -Wnon-virtual-dtor
# # This helps catch hard to track down memory errors
# # warn if float is implicit promoted to double:
# WARNINGS_CFLAGS += -Wdouble-promotion
# # warn on security issues around functions that format output (ie printf):
# WARNINGS_CFLAGS += -Wformat=2
# # warn for c-style casts:
# WARNINGS_CFLAGS += -Wold-style-cast
# # warn on sign conversions
# WARNINGS_CFLAGS += -Wsign-conversion
# warn for potential performance problem casts
WARNINGS_CFLAGS += -Wcast-align
# warn on type conversions that may lose data
WARNINGS_CFLAGS += -Wconversion
# # warn if if/else branches have duplicated code. not available in g++ 6.3.0
# WARNINGS_CFLAGS += -Wduplicated-branches
# warn if if/else chain has duplicated conditions
WARNINGS_CFLAGS += -Wduplicated-cond
# # warn on statements that fallthrough without an explicit annotation
# # not available in g++ 6.3.0
# WARNINGS_CFLAGS += -Wimplicit-fallthrough
# warn about logical operations being used where bitwise were probably wanted
WARNINGS_CFLAGS += -Wlogical-op
# warn if indentation implies blocks where blocks do not exist
WARNINGS_CFLAGS += -Wmisleading-indentation
# warn if a null dereference is detected:
WARNINGS_CFLAGS += -Wnull-dereference
# warn if you overload (not override) a virtual function
WARNINGS_CFLAGS += -Woverloaded-virtual
# warn if non-standard C++ is used:
WARNINGS_CFLAGS += -Wpedantic
# warn on anything being unused:
WARNINGS_CFLAGS += -Wunused
# warn if you perform a cast to the same type:
WARNINGS_CFLAGS += -Wuseless-cast

ifeq ($(MAKECMDGOALS),cppcolleff-debug)
  CFLAGS = $(MACHINE) $(DBG_FLAG) $(DFLAGS)
else
  CFLAGS = $(MACHINE) $(OPT_FLAG) $(DFLAGS)
endif
CFLAGS += $(WARNINGS_CFLAGS)

LIBOBJECTS  = $(addprefix $(OBJDIR)/, $(LIBSOURCES_CPP:.cpp=.o))
BINOBJECTS  = $(addprefix $(OBJDIR)/, $(BINSOURCES_CPP:.cpp=.o))

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
	$(CXX) $(LDFLAGS) $(WL_FLAGS) $(BINOBJECTS) $(OBJDIR)/libcppcolleff.a $(LIBS) -o $@

tests: libcppcolleff
	$(CXX) $(LDFLAGS) $(WL_FLAGS) src/tests.cpp $(OBJDIR)/libcppcolleff.a $(LIBS) -o tests.x

$(LIBOBJECTS): | $(OBJDIR)

$(BINOBJECTS): | $(OBJDIR)

$(OBJDIR):
	mkdir $(OBJDIR)

install-cpp: uninstall-cpp all
	cp $(OBJDIR)/cppcolleff $(BINDEST_DIR)
	cp $(OBJDIR)/libcppcolleff.a $(LIBDEST_DIR)
	cp -r $(INCDIR)/cppcolleff $(INCDEST_DIR)

uninstall-cpp:
	-rm -rf $(BINDEST_DIR)/cppcolleff
	-rm -rf $(LIBDEST_DIR)/libcppcolleff.a
	-rm -rf $(INCDEST_DIR)/cppcolleff

install-py: uninstall-py
	$(MAKE) install -C $(PYTHON_PACKAGE_DIR)

uninstall-py:
	$(MAKE) uninstall -C $(PYTHON_PACKAGE_DIR)

install: clean install-cpp install-py

uninstall: uninstall-cpp uninstall-py

develop-install-py: develop-uninstall-py all
	$(MAKE) develop-install -C $(PYTHON_PACKAGE_DIR)

develop-uninstall-py:
	$(MAKE) develop-uninstall -C $(PYTHON_PACKAGE_DIR)

$(BINDEST_DIR):
	mkdir $(BINDEST_DIR)

$(LIBDEST_DIR):
	mkdir $(LIBDEST_DIR)

$(INCDEST_DIR):
	mkdir $(INCDEST_DIR)

clean:
	-rm -rf $(OBJDIR) cppcolleff cppcolleff-debug .depend *.out *.dat *~ *.o *.a
	$(MAKE) clean -C $(PYTHON_PACKAGE_DIR)


#### RULES ####

*.cpp: VERSION
	touch $(SRCDIR)/*.cpp
*.cc: VERSION
	touch $(SRCDIR)/*.cc
*.c: VERSION
	touch $(SRCDIR)/*.c

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) -c $(CFLAGS) $(INC) -I./$(SRCDIR) $< -o $@;
