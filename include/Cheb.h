/**
 * Cheb.h/Cheb.cc
 * Created at 2014.5.23
 *
 * This file provide a set of subroutines for producing
 * Chebyshev differential matrix
 * and evaluating derivatives of Chebyshev expanded functions.
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

#ifndef cheb_cheb_h
#define cheb_cheb_h

#include "armadillo"
#include "Boundary.h"

/**
 * The index of the Chebyshev-Gauss-Lobatto nodes is
 *          0, 1, 2, ..., N
 * Number of nodes is N+1
 */
class Cheb{
public:
    Cheb(arma::uword num_nodes): N(num_nodes-1){}
    int size() {return N+1;}
    arma::colvec x();  // Chebyshev-Gauss-Lobatto nodes
    arma::colvec w();  // Chebyshev-Gauss-Lobatto weights
    arma::colvec c();  // Chebyshev coefficients
    arma::mat D();  // 1st order Chebyshev differential matrix
    arma::mat D2();   // 2nd order Chebyshev differential matrix
    arma::mat D(Boundary, Boundary);  // subject to boundary condition
    arma::mat D2(Boundary, Boundary);  // subject to boundary condition
    arma::colvec clencurt_weights_fft();
    double quadrature_clencurt(const arma::colvec &f);
    double quadrature_clencurt(const arma::colvec &f, const arma::colvec &w);
    arma::colvec barycentric_weights();
    arma::mat barycentric_matrix(const arma::colvec &y);
    arma::colvec interpolate(const arma::colvec &y, const arma::colvec &f);

private:
    arma::uword N;
    // computed when first calling quadrature_clencurt().
    arma::colvec clencurt_weights;

    arma::mat D_dbc_dbc();
    arma::mat D_dbc_rbc(Boundary);
    arma::mat D_rbc_dbc(Boundary);
    arma::mat D_rbc_rbc(Boundary, Boundary);
    arma::mat D2_dbc_dbc();
    arma::mat D2_dbc_rbc(Boundary);
    arma::mat D2_rbc_dbc(Boundary);
    arma::mat D2_rbc_rbc(Boundary, Boundary);
};

#endif

