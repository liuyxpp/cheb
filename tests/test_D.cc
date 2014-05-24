#include <iostream>
#include "armadillo"
#include "Cheb.h"

using namespace std;
using namespace arma;

int main(){
    int N = 3;

    Cheb cheb(N+1);
    colvec x = cheb.x();
    x.print("x:");
    colvec w = cheb.w();
    w.print("w:");
    colvec c = cheb.c();
    c.print("c:");
    mat D = cheb.D();
    D.print("D:");
    mat D2 = cheb.D2();
    D2.print("D2:");

    return 0;
}
