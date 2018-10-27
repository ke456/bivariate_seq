#include "ZZ_bivariate_seq.h"
#include <NTL/BasicThreadPool.h>

#include <NTL/lzz_pX.h>
NTL_CLIENT

int main(){
    // set this to be the prime you want to work over
    zz_p::init(13);
    SetNumThreads(4);
 
    // this are "bivariate" polynomials
    Vec<Vec<long>> num;
    Vec<Vec<long>> den;
    
    Vec<long> tmp;
    /**** NUMERATOR ********************************/
    // encode num: 1-x
    tmp.append(1); tmp.append(-1);
    num.append(tmp);
    tmp = Vec<long>();
    
    /**** DENOMINATOR *******************************/
    // encode den: (1+x)^2 + (-x)y + y^2
    
    // (1+x)^2
    tmp.append(1); tmp.append(2); tmp.append(1);
    den.append(tmp);
    tmp = Vec<long>();
    
    // (-x)
    tmp.append(0); tmp.append(-1);
    den.append(tmp); 
    tmp = Vec<long>();
    
    // 1
    tmp.append(1);
    den.append(tmp);
    
    cout << "num: " << num << endl;
    cout << "den: " << den << endl;
    
    // this creates the object
    bivariate_lin_seq bls{num,den,2,2};
    
    long a = 3;
    long b = 10;
    long L = 1;
    
    cout << "a: " << a << endl;
    cout << "b: " << b << endl;
    cout << "L: " << L << endl;
    
    //Vec<ZZ> enZZ, edZZ;
		//double tw = GetWallTime();
    //bls.get_entry_sq_ZZ(enZZ,edZZ,a,b,L,60,3);
		//cout << "regular: " << GetWallTime() - tw << endl;
		//tw = GetWallTime();

		double tw = GetWallTime();
		Vec<Vec<ZZ>> entries;
		bls.get_entry_sq_ZZ_geometric(entries,a,b,L,15,3);
		cout << "geo: " << GetWallTime() - tw << endl;
    cout << "entries: " << entries << endl; 
}





















