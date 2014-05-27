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

    // DBC-DBC
    Boundary dbc = Boundary(); // default: DBC
    mat Ddbc = cheb.D(dbc, dbc);
    Ddbc.print("D with homo DBC:");
    mat D2dbc = cheb.D2(dbc, dbc);
    D2dbc.print("D2 with homo DBC:");

    // RBC-RBC
    Boundary rbc = Boundary(1, 1, 0); // RBC
    mat Drbc = cheb.D(rbc, rbc);
    Drbc.print("D with homo RBC:");
    mat D2rbc = cheb.D2(rbc, rbc);
    D2rbc.print("D2 with homo RBC:");

    // NBC-NBC
    Boundary nbc = Boundary(1, 0, 0); // NBC
    mat Dnbc = cheb.D(nbc, nbc);
    Dnbc.print("D with homo NBC:");
    mat D2nbc = cheb.D2(nbc, nbc);
    D2nbc.print("D2 with homo NBC:");

    // DBC-RBC
    mat Ddr = cheb.D(dbc, rbc);
    Ddr.print("D with DBC-RBC:");
    mat D2dr = cheb.D2(dbc, rbc);
    D2dr.print("D2 with DBC-RBC:");

    // DBC-NBC
    mat Ddn = cheb.D(dbc, nbc);
    Ddn.print("D with DBC-NBC:");
    mat D2dn = cheb.D2(dbc, nbc);
    D2dn.print("D2 with DBC-NBC:");

    // RBC-DBC
    mat Drd = cheb.D(rbc, dbc);
    Drd.print("D with RBC-DBC:");
    mat D2rd = cheb.D2(rbc, dbc);
    D2rd.print("D2 with RBC-DBC:");

    // NBC-DBC
    mat Dnd = cheb.D(nbc, dbc);
    Dnd.print("D with NBC-DBC:");
    mat D2nd = cheb.D2(nbc, dbc);
    D2nd.print("D2 with NBC-DBC:");

    return 0;
}
