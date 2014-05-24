/**
 * Boundary.h/Boundary.cc
 * Created at 2014.5.23
 *
 * This file defines boundary conditions.
 *
 * Copyright (C) 2014 Yi-Xin Liu <liuyxpp@gmail.com>
 *
 * This file is part of cheb++
 *
 * cheb++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version.
 *
 * cheb++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with cheb++. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef cheb_boundary_h
#define cheb_boundary_h

enum class BC {DBC, NBC, RBC};

/**
 * The boundary condition can be generally written as
 *      alpha * du/dx + beta * u = gamma
 */
class Boundary{
public:
    Boundary(): bc(BC::DBC), a1(0), a2(1.0), a3(0) {}
    Boundary(double, double, double);
    BC kind() {return bc;}
    double alpha() {return a1;}
    double beta() {return a2;}
    double gamma() {return a3;}
private:
    BC bc;
    double a1, a2, a3;
};

#endif

