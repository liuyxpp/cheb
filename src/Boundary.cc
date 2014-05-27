#include "Boundary.h"

/**
 * Initialization of Boundary by three coefficients.
 */
Boundary::Boundary(const double k1, const double k2, const double k3){
    if(k1 == 0 && k2 != 0)
        bc = BC::DBC;
    else if(k1 != 0 && k2 == 0)
        bc = BC::NBC;
    else if(k1 == 0 && k2 == 0)
        bc = BC::PBC;
    else
        bc = BC::RBC;
    a1 = k1;
    a2 = k2;
    a3 = k3;
}

