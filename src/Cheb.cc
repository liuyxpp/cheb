#include "Cheb.h"
#include "armadillo"

using namespace arma;

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
 * Generate a set of Chebyshev coefficients.
 */
colvec Cheb::c(){
    colvec cc(N+1, fill::ones);
    cc(0) = 2.0;
    cc(N) = 2.0;
    return cc;
}

/**
 * Generate 1st order Chebyshev differentiation matrix.
 */
mat Cheb::D(){
    colvec cc = c();
    for(int i=1; i<=N; i+=2){
        cc(i) = -cc(i);
    }

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

