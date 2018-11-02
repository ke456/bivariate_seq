#include "ZZ_bivariate_seq.h"
#include <NTL/BasicThreadPool.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <NTL/lzz_pX.h>
NTL_CLIENT

void read (Vec<Vec<long>> &num, 
           Vec<Vec<long>> &den,
					 long &a, long &b, long &L,
					 long &nbits, long &iter,
					 int &d1, int &d2,
					 istream &in){

    string line;

    // read numerator
    long coeff;
    while(getline(in, line)){
        if (line ==  "") break;
        num.append(Vec<long> ());
        istringstream iss{line};
        while (iss >> coeff) num[num.length()-1].append(coeff);
    }

    // read denominator
    while(getline(in, line)){
        if (line ==  "") break;
        den.append(Vec<long> ());
        istringstream iss{line};
        while (iss >> coeff) den[den.length()-1].append(coeff);
    }

    //cout << "num: " << num << endl;
    //cout << "den: " << den << endl;

		in >> a >> b >> L >> nbits >> iter >> d1 >> d2;

		//cout << "a: " << a << endl;
		//cout << "b: " << b << endl;
		//cout << "L: " << L << endl;
		//cout << "d1: " << d1 << endl;
		//cout << "d2: " << d2 << endl;
}

int main(int argc, char *argv[]){
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

    long a = 3;
    long b = 2000;
    long L = 1;
    long nbits = 60;
    long iter = 5;
		int d1 = 2;
		int d2 = 2;
		if (argc == 2){
			num = Vec<Vec<long>>();
			den = Vec<Vec<long>>();
			ifstream fin{argv[1]};
			read(num,den,a,b,L,nbits,iter,d1,d2,fin);
		}
    if (argc >= 4){
        a = atoi(argv[1]);
        b = atoi(argv[2]);
        L = atoi(argv[3]);
        if (argc >= 5) nbits = atoi(argv[4]);
        if (argc >= 6) iter = atoi(argv[5]);
    }
    cout << "num: " << num << endl;
		cout << "den: " << den << endl;
		cout << "a: " << a << endl;
    cout << "b: " << b << endl;
    cout << "L: " << L << endl;
		cout << "d1: " << d1 << endl;
		cout << "d2: " << d2 << endl;

    bivariate_lin_seq bls{num,den,d1,d2};
    
    double tw = GetWallTime();
    Vec<Vec<ZZ>> entries;
    bls.get_entry_sq_ZZ_geometric(entries,a,b,L,nbits,iter);
    cout << "geo: " << GetWallTime() - tw << endl;
    cout << "entries: " << entries << endl; 
}





















