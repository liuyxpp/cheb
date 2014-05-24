#include "Boundary.h"

/**
 * Initialization of Boundary by three coefficients.
 */
Boundary::Boundary(double k1, double k2, double k3){
    if(k1 == 0)
        bc = BC::DBC;
    else if(k2 == 0)
        bc = BC::NBC;
    else
        bc = BC::RBC;
    a1 = k1;
    a2 = k2;
    a3 = k3;
}

