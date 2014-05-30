#include <iostream>
#include "armadillo"
#include "Cheb.h"

using namespace std;
using namespace arma;

void test_clencurt_weights_fft(){
    uword N = 8;
    Cheb cheb(N+1);
    colvec w = cheb.clencurt_weights_fft();
    w.print("w =");

    N = 7;
    Cheb cheb1(N+1);
    colvec w1 = cheb1.clencurt_weights_fft();
    w1.print("w1 =");
}

/**
 * Four funcitons in [-1, 1] are tested:
        f(x) = |x|^3,           I = .5
        f(x) = exp(-x^(-2)),    I = 2*(exp(-1) + sqrt(pi)*(erf(1) - 1))
        f(x) = 1/(1+x^2),       I = pi/2
        f(x) = x^10,            I = 2/11
 *
 */
void test_quadrature_clencurt(){
    uword Nmax = 25;
    colvec e1 = zeros<colvec>(Nmax-1);
    colvec e2 = zeros<colvec>(Nmax-1);
    colvec e3 = zeros<colvec>(Nmax-1);
    colvec e4 = zeros<colvec>(Nmax-1);
    double I1 = 0.5;
    double I2 = 2 * (exp(-1) + sqrt(datum::pi) * (erf(1.0) - 1.0));
    double I3 = 0.5 * datum::pi;
    double I4 = 2.0 / 11;
    for(uword N=2; N<Nmax+1; N++){
        Cheb cheb(N+1);
        colvec x = cheb.x();
        colvec f1 = abs(x) % abs(x) % abs(x);
        e1(N-2) = abs(cheb.quadrature_clencurt(f1) - I1);
        colvec f2 = exp(-1.0 / (x%x));
        e2(N-2) = abs(cheb.quadrature_clencurt(f2) - I2);
        colvec f3 = 1.0 / (1 + x % x);
        e3(N-2) = abs(cheb.quadrature_clencurt(f3) - I3);
        colvec f4 = pow(x, 10);
        e4(N-2) = abs(cheb.quadrature_clencurt(f4) - I4);
    }
    e1.print("e1 =");
    e2.print("e2 =");
    e3.print("e3 =");
    e4.print("e4 =");
}

int main(){
    //test_clencurt_weights_fft();
    test_quadrature_clencurt();

    return 0;
}
