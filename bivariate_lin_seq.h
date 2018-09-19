#ifndef BIVARIATE_LIN_SEQ
#define BIVARIATE_LIN_SEQ

#include <NTL/lzz_pX.h>

NTL_CLIENT

// give min poly P , returns D-th element
// note that P is the reverse of the denominator
zz_p get_elem(const long& D, const zz_pX &P, const Vec<zz_p>& init);

// returns the n initial conditions of num/den
Vec<zz_p> get_init(const long& n, const zz_pX &num, const zz_pX &den);

class bivariate_lin_seq{
    Vec<zz_pX> polX_num; // coefficients in X for numerator
    Vec<zz_pX> polX_den; // coefficients in X for denominator
    
    const int d1,d2; // deg in x,y
    
    void eval_x (zz_pX &res, const zz_p& x, const Vec<zz_pX> &poly);
    
public:
    // num and den represent bivariate polynomials, where each coefficient
    // is a polynomial in x. For example, (1+x) + (-x)y + y^2 => [[1,1], [0,-1], [1]]
    // d1 is degree in x, d2 is degree in y
    bivariate_lin_seq(const Vec<zz_pX>& num, Vec<zz_pX>& den,
                      int d1, int d2);
                      
    // find the generating series for the D-th row
    void find_row(zz_pX &num, zz_pX &den, const long& D);
    
};
#endif 
