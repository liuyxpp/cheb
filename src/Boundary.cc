#include "Boundary.h"

/**
 * Initialization of Boundary by three coefficients.
 * NOTE: in this implementation,
 *          DBC: u = 0
 *          NBC: du/dx = 0
 * Therefore, it does not support
 *          DBC: u = c (c is a constant)
 *          NBC: du/dx = c (c is a constant)
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

Boundary& Boundary::operator=(const Boundary &rhs){
    bc = rhs.bc;
    a1 = rhs.a1;
    a2 = rhs.a2;
    a3 = rhs.a3;
}

