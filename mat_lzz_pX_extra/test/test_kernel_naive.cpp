#include <NTL/BasicThreadPool.h>
#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <numeric>
#include <random>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

std::ostream &operator<<(std::ostream &out, const std::vector<long> &s)
{
    out << "[ ";
    for (auto &i: s)
        out << i << " ";
    return out << "]";
}

void check(long m, long n, long d){
    Mat<zz_pX> pmat;
    random(pmat, m, n, d+1);
    std::cout << "Computing kernel basis, dimensions " << m << " x " << n << ", degree " << d << std::endl;

    double t1w,t2w;

    // declare shifts
    Shift shift1(m,0); // uniform [0,...,0]
    Shift shift2(m); // increasing [0,1,2,..,m-1]
    std::iota(shift2.begin(), shift2.end(),0);
    Shift shift3(m); // decreasing [m,..,3,2,1]
    for (long i = 0; i < m; ++i)
        shift3[i] = m - i;
    Shift shift4(m); // random shuffle of [0,1,...,m-1]
    std::iota(shift4.begin(), shift4.end(),0);
    std::shuffle(shift4.begin(), shift4.end(), std::mt19937{std::random_device{}()});
    Shift shift5(m); // Hermite shift
    for (long i = 0; i < m; ++i)
        shift5[i] = n*d*i;
    Shift shift6(m); // reverse Hermite shift
    for (long i = 0; i < m; ++i)
        shift6[i] = (n*d+1)*(m-1-i);
    Shift shift7(m);
    for (long i = 0; i < m; ++i)
        if (i>=m/2)
            shift7[i] = n*d;
    Shift shift8(m,0);
    shift8[m-1] = n*d;

    std::vector<Shift> shifts = {shift1, shift2, shift3, shift4, shift5, shift6, shift7, shift8};

    for (auto shift: shifts)
    {
        if (m<40)
            std::cout << "--shift =\t" << shift << std::endl;
        else
            std::cout << "--shift =\t" << "<length " << m << ">" << std::endl;
        Mat<zz_pX> kerbas;
        t1w = GetWallTime();
        kernel_basis_via_approximation(kerbas, pmat, shift);
        t2w = GetWallTime();
        cout << "time (kernel_naive): " << t2w-t1w << endl;

        t1w = GetWallTime();
        bool correct = is_kernel_basis(kerbas, pmat, shift, ORD_WEAK_POPOV, true, false);
        t2w = GetWallTime();
        std::cout << "Verification: " << (correct ? "correct" : "wrong") << ", time " << (t2w-t1w) << std::endl << std::endl;
    }
}

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    SetNumThreads(4);

    zz_p::init(NTL::GenPrime_long(60));

    if (argc==4)
    {
        check(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
    }
    else
        throw std::invalid_argument("Usage: ./test_kernel_naive rdim cdim degree");


    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s