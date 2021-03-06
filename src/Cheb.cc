#include "Cheb.h"
#include "Boundary.h"
#include "armadillo"

#include <iostream>
using namespace std;

using namespace arma;

const double PI=3.14159265358979323846264338327950288;

/**
 * Generate a set of Chebyshev-Gauss-Lobatto nodes.
 */
colvec Cheb::x(){
    colvec i = linspace(0, N, N+1);
    return cos(PI * i / N);
}

/**
 * Generate a set of Chebyshev-Gauss-Lobatto weights.
 */
colvec Cheb::w(){
    colvec ww(N+1);
    ww.fill(PI / N);
    ww(0) = 0.5 * ww(0);
    ww(N) = 0.5 * ww(N);
    return ww;
}

/**
 * Generate a set of Chebyshev coefficients. Often denoted as $\gamma_n$.
 */
colvec Cheb::c(){
    colvec cc(N+1, fill::ones);
    cc(0) = 2.0;
    cc(N) = 2.0;
    return cc;
}

/**
 * Generate a set of Barycentric weights for Chebyshev-Gauss-Lobatto grid.
 * The weigts are
 *      w_j = (-1)^j * d_j, with
 *      d_j = 1/2,  j=0 or j=N
 *      d_j = 1,    j = 1, 2, ..., N-1
 */
colvec Cheb::barycentric_weights(){
    colvec ww = ones<colvec>(N+1);
    for(uword i=0; i<=N; i+=2)
        ww(i) *= -1;
    ww(0) *= 0.5;
    ww(N) *= 0.5;
    return ww;
}

/**
 * Compute matrix T_{kj} for interpolation between two sets of points.
 * y is a colvec which contains a set of locations to be interpolated.
 * y should be in the range of [-1, 1]
 * k is the index of colvec y.
 * j is the index of colvec x.
 * The size of T is M x (N+1), where M = size(y).
 */
 mat Cheb::barycentric_matrix(const colvec &y){
    uword M = y.n_elem;
    colvec xx = x();
    colvec ww = barycentric_weights();
    mat T = zeros<mat>(M, N+1);
    bool row_has_match;

    for(uword k=0; k<M; k++){
        row_has_match = false;
        for(uword j=0; j<=N; j++){
            double yt = y(k);
            double xt = xx(j);
            if(yt == xt || abs(yt-xt) < abs(min(yt, xt))*datum::eps){
                row_has_match = true;
                T(k, j) = 1.0;
                continue;
            }
        }
        if(!row_has_match){
            colvec t = ww / (y(k) - xx);
            T.row(k) = t.t() / sum(t);
        }
    }

    return T;
 }

 /**
  * Interpolate f on the set of locations y.
  */
colvec Cheb::interpolate(const colvec &y, const colvec &f){
    mat T = barycentric_matrix(y);
    return T * f;
}

/**
 * Generate 1st order Chebyshev differentiation matrix.
 */
mat Cheb::D(){
    colvec cc = c();
    for(uword i=1; i<=N; i+=2)
        cc(i) = -cc(i);

    colvec xx = x();
    mat X = repmat(xx, 1, N+1);
    mat dX = X - X.t();

    mat I = eye<mat>(N+1, N+1);
    mat DD = (cc * (1/cc).t()) / (dX + I);
    DD.diag() -= sum(DD, 1);
    return DD;
}

/**
 * Generate 2nd order Chebyshev differentiation matrix.
 */
mat Cheb::D2(){
    mat DD = D();
    mat DD2 = DD * DD;
    DD2.diag() = zeros<colvec>(N+1);
    DD2.diag() -= sum(DD2, 1);
    return DD2;
}

/**
 * Generate 1st order Chebyshev differentiation matrix
 * subject to boundary condition.
 * Note that this method requires N >= 2, i.e. at least 3 Chebyshev points.
 */
mat Cheb::D(Boundary lb, Boundary rb){
    BC lbc = lb.kind();
    BC rbc = rb.kind();
    if(lbc==BC::DBC && rbc==BC::DBC)
        return D_dbc_dbc();
    else if(lbc==BC::DBC && (rbc==BC::RBC || rbc==BC::NBC))
        return D_dbc_rbc(rb);
    else if((lbc==BC::RBC || lbc==BC::NBC) && rbc==BC::DBC)
        return D_rbc_dbc(lb);
    else
        // can be either RBC or NBC for each boundary
        return D_rbc_rbc(lb, rb);
}

/**
 * Generate 2nd order Chebyshev differentiation matrix
 * subject to boundary condition.
 * Note that this method requires N >= 2, i.e. at least 3 Chebyshev points.
 */
mat Cheb::D2(Boundary lb, Boundary rb){
    BC lbc = lb.kind();
    BC rbc = rb.kind();
    if(lbc==BC::DBC && rbc==BC::DBC)
        return D2_dbc_dbc();
    else if(lbc==BC::DBC && (rbc==BC::RBC || rbc==BC::NBC))
        return D2_dbc_rbc(rb);
    else if((lbc==BC::RBC || lbc==BC::NBC) && rbc==BC::DBC)
        return D2_rbc_dbc(lb);
    else
        // can be either RBC or NBC for each boundary
        return D2_rbc_rbc(lb, rb);
}

/**
 * Compute the (N+1) weights for Clenshaw-Curtis quadrature
 * on the interval [-1, 1].
 * Using FFT algorithm, the weights are computed in linear time.
 */
colvec Cheb::clencurt_weights_fft(){
    colvec c = zeros<colvec>(N+1);
    // N even: c(0), c(2), ..., c(N)
    // N odd:  c(0), c(2), ..., c(N-1)
    for(uword i=0; i<N+1; i+=2){
        c(i) = 2.0 / (1.0 - i*i);
    }
    colvec cc = join_vert(c, flipud(c.subvec(1, N-1)));

    colvec v = real(ifft(cc*cx_double(1, 0)));

    colvec w = 2.0 * v.subvec(0, N);
    w(0) *= 0.5;
    w(N) *= 0.5;

    return w;
}

/**
 * Compute Clenshaw-Curtis quadrature.
 * Avoid recomputing Clenshaw-Curtis weights by storing them.
 */
double Cheb::quadrature_clencurt(const colvec &f){
    if(clencurt_weights.is_empty()){
        clencurt_weights.set_size(N+1);
        clencurt_weights = clencurt_weights_fft();
    }

    return quadrature_clencurt(f, clencurt_weights);
}

double Cheb::quadrature_clencurt(const colvec &f, const colvec &w){
    return dot(w, f);
}

/***************************************************************************
 *
 *  Private Member Functions
 *
 ***************************************************************************/

/**
 * Generate 1st order Chebyshev differentiation matrix
 * subject to homogeneous DBCs
 *      u(x=-1) = 0
 *      u(x=1) = 0
 * i.e. alpha = 0, gamma = 0.
 */
mat Cheb::D_dbc_dbc(){
    mat DD = D();

    return DD.submat(span(1, N-1), span(1, N-1));
}

/**
 * Generate 2nd order Chebyshev differentiation matrix
 * subject to homogeneous DBCs.
 */
mat Cheb::D2_dbc_dbc(){
    mat DD2 = D2();

    return DD2.submat(span(1, N-1), span(1, N-1));
}

/**
 * Generate 1st order Chebyshev differentiation matrix
 * subject to homogeneous RBCs.
 *      du/dx + beta * u = 0
 * i.e. alpha = 1, gamma = 0.
 */
mat Cheb::D_rbc_rbc(Boundary lb, Boundary rb){
    double lbeta = lb.beta();
    double rbeta = rb.beta();
    mat DD = D();
    colvec xx = x();
    mat I = eye<mat>(N+1, N+1);

    rowvec xjrow = 1 / (1 - xx.subvec(1, N-1).t() % xx.subvec(1, N-1).t());
    colvec xkcol0 = 1 - xx % xx;
    colvec xkcol1 = -2 * xx;

    mat fac0 = xkcol0 * xjrow;
    mat fac1 = xkcol1 * xjrow;

    mat D1t = DD;
    D1t.cols(1, N-1) = fac0 % DD.cols(1, N-1) + fac1 % I.cols(1, N-1);

    colvec omx = 0.5 * (1 - xx);
    colvec opx = 0.5 * (1 + xx);

    colvec r0 = opx + (0.5 + DD(0, 0) + rbeta) * xkcol0 / 2;
    colvec r1 = 0.5 - (0.5 + DD(0, 0) + rbeta) * xx;
    D1t.col(0) = r0 % DD.col(0) + r1 % I.col(0);

    colvec l0 = omx + (0.5 - DD(N, N) - lbeta) * xkcol0 / 2;
    colvec l1 = -0.5 + (DD(N, N) + lbeta - 0.5) * xx;
    D1t.col(N) = l0 % DD.col(N) + l1 % I.col(N);

    return D1t;
}

/**
 * Generate 2nd order Chebyshev differentiation matrix
 * subject to homogeneous RBCs.
 */
mat Cheb::D2_rbc_rbc(Boundary lb, Boundary rb){
    double lbeta = lb.beta();
    double rbeta = rb.beta();
    mat DD = D();
    mat DD2 = D2();
    colvec xx = x();
    mat I = eye<mat>(N+1, N+1);

    rowvec xjrow = 1 / (1 - xx.subvec(1, N-1).t() % xx.subvec(1, N-1).t());
    colvec xkcol0 = 1 - xx % xx;
    colvec xkcol1 = -2 * xx;
    colvec xkcol2 = -2 * ones<colvec>(N+1);

    mat fac0 = xkcol0 * xjrow;
    mat fac1 = xkcol1 * xjrow;
    mat fac2 = xkcol2 * xjrow;

    mat D2t = DD2;
    D2t.cols(1, N-1) = fac0 % DD2.cols(1, N-1)
                       + 2 * fac1 % DD.cols(1, N-1)
                       + fac2 % I.cols(1, N-1);

    colvec omx = 0.5 * (1 - xx);
    colvec opx = 0.5 * (1 + xx);

    colvec r0 = opx + (0.5 + DD(0, 0) + rbeta) * xkcol0 / 2;
    colvec r1 = 0.5 - (0.5 + DD(0, 0) + rbeta) * xx;
    double r2 = -0.5 - DD(0, 0) - rbeta;
    D2t.col(0) = r0 % DD2.col(0) + 2 * r1 % DD.col(0) + r2 * I.col(0);

    colvec l0 = omx + (0.5 - DD(N, N) - lbeta) * xkcol0 / 2;
    colvec l1 = -0.5 + (DD(N, N) + lbeta - 0.5) * xx;
    double l2 = -0.5 + DD(N, N) + lbeta;
    D2t.col(N) = l0 % DD2.col(N) + 2 * l1 % DD.col(N) + l2 * I.col(N);

    return D2t;
}

/**
 * Generate 1st order Chebyshev differentiation matrix
 * subject to DBC at x=-1 and RBC (or NBC) at x=1.
 */
mat Cheb::D_dbc_rbc(Boundary rb){
    double rbeta = rb.beta();
    mat DD = D();
    colvec xx = x();
    mat I = eye<mat>(N+1, N+1);

    rowvec xjrow = 1 - xx.subvec(1, N-1).t();
    colvec xkcol = 1 - xx.subvec(0, N-1);
    colvec oner = ones<colvec>(N);

    mat fac0 = oner * (1 / xjrow);
    mat fac1 = xkcol * (1 / xjrow);

    mat D1t(N, N);
    D1t.cols(1, N-1) = fac1 % DD(span(0, N-1), span(1, N-1))
                       - fac0 % I(span(0, N-1), span(1, N-1));

    double cfac = DD(0, 0) + rbeta;
    D1t.col(0) = -cfac * I(span(0, N-1), 0)
                 + (1 + cfac * xkcol) % DD(span(0, N-1), 0);

    return D1t;
}

/**
 * Generate 2nd order Chebyshev differentiation matrix
 * subject to DBC at x=-1 and RBC (or NBC) at x=1.
 */
mat Cheb::D2_dbc_rbc(Boundary rb){
    double rbeta = rb.beta();
    mat DD = D();
    mat DD2 = D2();
    colvec xx = x();
    //mat I = eye<mat>(N+1, N+1);

    rowvec xjrow = 1 - xx.subvec(1, N-1).t();
    colvec xkcol = 1 - xx.subvec(0, N-1);
    colvec oner = ones<colvec>(N);

    mat fac0 = oner * (1 / xjrow);
    mat fac1 = xkcol * (1 / xjrow);

    mat D2t(N, N);
    D2t.cols(1, N-1) = fac1 % DD2(span(0, N-1), span(1, N-1))
                       - 2 * fac0 % DD(span(0, N-1), span(1, N-1));

    double cfac = DD(0, 0) + rbeta;
    D2t.col(0) = -2 * cfac * DD(span(0, N-1), 0)
                 + (1 + cfac * xkcol) % DD2(span(0, N-1), 0);

    return D2t;
}

/**
 * Generate 1st order Chebyshev differentiation matrix
 * subject to RBC at x=-1 and DBC (or NBC) at x=1.
 */
mat Cheb::D_rbc_dbc(Boundary lb){
    double lbeta = lb.beta();
    mat DD = D();
    colvec xx = x();
    mat I = eye<mat>(N+1, N+1);

    rowvec xjrow = 1 + xx.subvec(1, N-1).t();
    colvec xkcol = 1 + xx.subvec(1, N);
    colvec oner = ones<colvec>(N);

    mat fac0 = oner * (1 / xjrow);
    mat fac1 = xkcol * (1 / xjrow);

    mat D1t(N, N);
    D1t.cols(0, N-2) = fac1 % DD(span(1, N), span(1, N-1))
                       + fac0 % I(span(1, N), span(1, N-1));

    double cfac = DD(N, N) + lbeta;
    D1t.col(N-1) = -cfac * I(span(1, N), N)
                 + (1 - cfac * xkcol) % DD(span(1, N), N);

    return D1t;
}

/**
 * Generate 2nd order Chebyshev differentiation matrix
 * subject to RBC at x=-1 and DBC (or NBC) at x=1.
 */
mat Cheb::D2_rbc_dbc(Boundary lb){
    double lbeta = lb.beta();
    mat DD = D();
    mat DD2 = D2();
    colvec xx = x();
    //mat I = eye<mat>(N+1, N+1);

    rowvec xjrow = 1 + xx.subvec(1, N-1).t();
    colvec xkcol = 1 + xx.subvec(1, N);
    colvec oner = ones<colvec>(N);

    mat fac0 = oner * (1 / xjrow);
    mat fac1 = xkcol * (1 / xjrow);

    mat D2t(N, N);
    D2t.cols(0, N-2) = fac1 % DD2(span(1, N), span(1, N-1))
                       + 2 * fac0 % DD(span(1, N), span(1, N-1));

    double cfac = DD(N, N) + lbeta;
    D2t.col(N-1) = -2 * cfac * DD(span(1, N), N)
                 + (1 - cfac * xkcol) % DD2(span(1, N), N);

    return D2t;
}
