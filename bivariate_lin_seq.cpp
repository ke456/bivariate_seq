#include "bivariate_lin_seq.h"
#include <iostream>
#include <NTL/BasicThreadPool.h>

Vec<int> get_binary(const long &t){
  Vec<int> result;
	auto n = t;
	while (n != 0){
		result.append(n % 2);
		n = n/2;
	}
	return result;
} 

void repeated_sq_mod(zz_pX &result, const zz_pX &p,
                              const zz_pX &h, const long &k){
	auto bin = get_binary(k);
	zz_pX b = p % h;
	zz_pX p_mod = b;
	for (long i = bin.length()-2; i>=0; i--){
		MulMod(b,b,b,h);
		if (bin[i] == 1)
			MulMod(b,b,p_mod,h);
	}
	result = b;
}

zz_p get_elem (const long& D, const zz_pX &P, const Vec<zz_p>& init){
    zz_pX mod;
	SetCoeff(mod,1,1);
	repeated_sq_mod(mod, mod, P, D);
	zz_p result{0};
	for (long i = 0; i < init.length(); i++)
			result = result + coeff(mod,i) * init[i];
    return result;
}

Vec<zz_p> get_init(const long& n, const zz_pX &num, const zz_pX &den){
    zz_pX partial_series;
    MulTrunc(partial_series, num,
            InvTrunc(den, n),n);
   Vec<zz_p> init;
   for (long t = 0; t <= deg(partial_series); t++)
        init.append(coeff(partial_series,t));
   return init;
}

bivariate_lin_seq::bivariate_lin_seq(const Vec<Vec<long>>& num,Vec<Vec<long>>& den, int d1, int d2): 
num_coeffs{num}, den_coeffs{den},d1{d1},d2{d2}{}

void create_poly(Vec<zz_pX> &res, const Vec<Vec<long>> &coeffs){
    res.SetLength(coeffs.length());

    for (long i = 0; i < coeffs.length(); i++){
        zz_pX tmp;
        for (long j = 0; j < coeffs[i].length(); j++){
            SetCoeff(tmp,j,zz_p(coeffs[i][j]));
        }
        res[i] = tmp;
    }
}

void bivariate_lin_seq::eval_x(zz_pX &res, const zz_p& x, const Vec<zz_pX> &poly){
    for (int i = 0; i < poly.length(); i++){
        zz_p val;
        eval(val,poly[i],x);
        SetCoeff(res, i, val);
    }  
}

void bivariate_lin_seq::find_row(zz_pX &num, zz_pX &den, const long& D){
    long degree = (D+1) * d1; // total degree on the bottom
    zz_p x_i = zz_p(0);
    Vec<zz_p> pointsX;
    Vec<zz_p> pointsY;
    
    pointsX.SetLength(degree);
    pointsY.SetLength(degree);
    
    Vec<zz_pX> polX_num, polX_den;
    
    create_poly(polX_num, num_coeffs);
    create_poly(polX_den, den_coeffs);
    
    for (long i = 0; i < degree; i++){
        zz_pX eval_num;
        zz_pX eval_den;
        do{
            eval_x(eval_num, x_i, polX_num);
            eval_x(eval_den, x_i, polX_den);
            x_i += zz_p(1);
        }while(ConstTerm(eval_den) == zz_p(0));
        
        // find initial conditions
        Vec<zz_p> init = get_init(d2,eval_num,eval_den);
        
        // find the point
        auto rp = get_elem(D,reverse(eval_den),init);
        auto p_pow = power(ConstTerm(eval_den), D+1);
        
        pointsX[i] = (x_i - zz_p(1));
        pointsY[i] = (rp*p_pow);
    }
    
    // interpolate
    interpolate(num, pointsX, pointsY);
    power(den,polX_den[0],D+1);
}

void bivariate_lin_seq::get_entry_sq_ZZ 
                        (Vec<ZZ> &entries_num, 
                         Vec<ZZ> &entries_den,
                         const long a, 
                         const long b, 
                         const long L,
                         const long nbits,
                         const long steps){
     Vec<long> primes;
     primes.SetLength(steps);
     
     long init_prime = GenPrime_long(nbits);
     primes[0] = init_prime;
     
     Vec<Vec<Vec<zz_p>>> coeffs;
     coeffs.SetLength(steps);
     
     for (long i = 1; i < steps; i++){
        primes[i] = NextPrime(primes[i-1]+1);
     }
        
    for (long i = 0; i < steps; i++){
        cout << "prime: " << primes[i] << endl;
        zz_p::init(primes[i]);
    
        zz_pX n,d;
        
        coeffs[i].SetLength(2*L+1);
        long at_D = 0;
        for (long D = b-L; D < b+L+1; D++, at_D++){
            find_row(n,d,D);
            cout << "D: " << D << endl;
            cout << "n: " << n << endl;
            cout << "d: " << d << endl;
            Vec<zz_p> init = get_init(deg(d), n, d);
            coeffs[i][at_D].SetLength(2*L+1);
            
             long at_N = 0;
            for (long N = a-L; N < a+L+1; N++, at_N++)
                coeffs[i][at_D][at_N] = get_elem(N, reverse(d), init);
            
            cout << "coeffs[" << i << "][" << at_D << "]=" 
                 << coeffs[i][at_D] << endl;
            
        }   
        cout << endl;     
    }
    
     
}



// g++ -std=c++14 *.cpp -g -lntl -lgmp -lm -lpthread -march=native














