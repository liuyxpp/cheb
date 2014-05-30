#include <iostream>
#include "armadillo"
#include "Cheb.h"
#include "CMatFile.h"

using namespace std;
using namespace arma;

void save_mat(colvec &x, colvec &y, colvec &f, colvec &f0, colvec &fi){
    mwSize N1 = (mwSize) x.n_elem;
    mwSize N2 = (mwSize) y.n_elem;
    mwSize n_bytes_1 = N1 * sizeof(double);
    mwSize n_bytes_2 = N2* sizeof(double);
    mwSize dim1[1] = {N1};
    mwSize dim2[1] = {N2};

    CMatFile mat;
    mat.matInit("test_interpolation.mat", "u");
    if(!mat.queryStatus()){
        mat.matPut("x", x.memptr(), n_bytes_1, 1, dim1,
                mxDOUBLE_CLASS, mxREAL);
        mat.matPut("y", y.memptr(), n_bytes_2, 1, dim2,
                mxDOUBLE_CLASS, mxREAL);
        mat.matPut("f", f.memptr(), n_bytes_1, 1, dim1,
                mxDOUBLE_CLASS, mxREAL);
        mat.matPut("f0", f0.memptr(), n_bytes_2, 1, dim2,
                mxDOUBLE_CLASS, mxREAL);
        mat.matPut("fi", fi.memptr(), n_bytes_2, 1, dim2,
                mxDOUBLE_CLASS, mxREAL);
    }
}

colvec fx(colvec &x){
    //return abs(x) + 0.5 * x - x%x;
    return 1 / (1 + 16 * x % x);
}

/**
 * Test function:
    f(x) = |x| + x/2 - x^2
 *      N       M       error
 *      20      1000    0.297868
 *      64      1000    0.009313
 *      128     1000    0.004631
 *      64      100     0.007500
    f(x) = 1/(1+16x^2)
 *      N       M       error
 *      20      1000    6.67107e-3
 *      64      1000    1.24193e-7
 *      128     1000    1.68754e-14
 *      64      100     1.18828e-7
 *
 */
void test_interpolation_1d(){
    uword N = 128;  // Number of CGL grid nodes.
    uword M = 1000;  // Number of regular grid nodes.

    Cheb cheb(N+1);
    colvec x = cheb.x();
    colvec f = fx(x);  // function to be interpolated.

    colvec y = linspace(-1, 1, M);
    colvec f0 = fx(y);

    colvec fi = cheb.interpolate(y, f);

    cout<<"Max error ="<<max(abs(fi - f0))<<endl;

    save_mat(x, y, f, f0, fi);
}

int main(){
    test_interpolation_1d();

    return 0;
}
