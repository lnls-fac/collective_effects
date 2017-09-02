// cppcolleff - Particle tracking code
// Copyright (C) 2015  LNLS Accelerator Physics Group
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

%module cppcolleff

%{
#include <cppcolleff/essentials.h>
#include <cppcolleff/Bunch.h>
#include <cppcolleff/Ring.h>
#include <cppcolleff/Results.h>
#include <cppcolleff/Feedback.h>
#include <cppcolleff/Wake.h>
#include <cppcolleff/cppcolleff.h>
%}

%include "../include/cppcolleff/essentials.h"
%include "../include/cppcolleff/Bunch.h"
%include "../include/cppcolleff/Ring.h"
%include "../include/cppcolleff/Results.h"
%include "../include/cppcolleff/Feedback.h"
%include "../include/cppcolleff/Wake.h"
%include "../include/cppcolleff/cppcolleff.h"


%include "std_vector.i"
%include "stl.i"

namespace std {
    %template(my_Ivector) vector<int>;
    %template(my_Dvector) vector<double>;
    %template(my_Dmatrix) vector<my_Dvector>;
    %template(my_Cvector) vector<complex<double>>;
    %template(my_PartVector) vector<Particle_t>;
    %template(my_PartMatrix) vector<my_PartVector>;
}
