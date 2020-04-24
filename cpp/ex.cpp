#include <iostream>
#include "affine_point.cpp"
#include <chrono>
#include "projective_point.cpp"
#include "json/json.h"
#include <fstream>
#include <iomanip>
#include <cassert>

using namespace std;

mpz_t modulo;
mpz_t a;
mpz_t b;
mpz_t curveOrder;
gmp_randclass r(gmp_randinit_default);

typedef AffinePoint<a, b, modulo> Point;
typedef ProjectivePoint<a, b, modulo> PPoint;
typedef Field<curveOrder> CurveField;

unsigned long long int getSeed() {
    unsigned long long int seed = 0;
    ifstream urandom("/dev/urandom", ios::in|ios::binary);
    urandom.read(reinterpret_cast<char*>(&seed), sizeof(seed));
    return seed;
}

void init(mpz_class a_v, mpz_class b_v, mpz_class mod, mpz_class curve) {
   mpz_init_set(modulo, mod.get_mpz_t());
   mpz_init_set(a, a_v.get_mpz_t());
   mpz_init_set(b, b_v.get_mpz_t());
   mpz_init_set(curveOrder, curve.get_mpz_t());
}

void step(Point gen, Point y, Point &A, CurveField  &alpha, CurveField &beta){
    mpz_t tmp;
    mpz_init_set(tmp,A.getX().getValue().get_mpz_t());
    mpz_mod(tmp,tmp,mpz_class(3).get_mpz_t());
    CurveField one = CurveField(1);
    switch (mpz_cmp_si(tmp, 1))
    {
    case -1:
        if(A.getY().getValue().get_mpz_t() == 0){
            A = A + y;
            alpha = alpha;
            beta = beta + one;
        }
        else{
            A = A + A;
            alpha = alpha + alpha;
            beta = beta + beta;
        }
        break;
    case 0:
        A = A + y;
        alpha = alpha;
        beta = beta + one;
        break;
    case 1:
        A = A + gen;
        alpha = alpha + one;
        beta = beta;
        break;
    }
    mpz_clear(tmp);
}


void step1(PPoint gen, PPoint y, PPoint &A, CurveField  &alpha, CurveField &beta){
    mpz_t tmp;
    mpz_init_set(tmp,A.getX().getValue().get_mpz_t());
    mpz_mod(tmp,tmp,mpz_class(3).get_mpz_t());
    CurveField one = CurveField(1);
    switch (mpz_cmp_si(tmp, 1))
    {
    case -1:
        if(A.getY().getValue().get_mpz_t() == 0){
            A = A + y;
            alpha = alpha;
            beta = beta + one;
        }
        else{
            A = A + A;
            alpha = alpha + alpha;
            beta = beta + beta;
        }
        break;
    case 0:
        A = A + y;
        alpha = alpha;
        beta = beta + one;
        break;
    case 1:
        A = A + gen;
        alpha = alpha + one;
        beta = beta;
        break;
    }
    mpz_clear(tmp); 
}

CurveField pollard_rho(Point gen, Point y){
    Point a = gen;
    Point b = gen;
    CurveField alpha_a = CurveField(1), alpha_b = CurveField(1);
    CurveField beta_a = CurveField(0), beta_b = CurveField(0);
    while(true){
        step(gen,y,a,alpha_a, beta_a);
        step(gen,y,b,alpha_b, beta_b);
        step(gen,y,b,alpha_b, beta_b);
        if (a==b){
            CurveField x = (alpha_b - alpha_a) / (beta_a - beta_b);
            return x;
        }
    }
}

CurveField pollard_rho1(PPoint gen, PPoint y){
    PPoint a = gen;
    PPoint b = gen;
    CurveField alpha_a = CurveField(1), alpha_b = CurveField(1);
    CurveField beta_a = CurveField(0), beta_b = CurveField(0);
    while(true){
        step1(gen,y,a,alpha_a, beta_a);
        step1(gen,y,b,alpha_b, beta_b);
        step1(gen,y,b,alpha_b, beta_b);
        if (a==b){
            CurveField x = (alpha_b - alpha_a) / (beta_a - beta_b);
            return x;
        }
    }
}

int main(int argc, char const *argv[]) {    
    if (argc < 2){
        cerr << "Usage: ./ex <curve_params.json>" << endl;
        return 0;
    }
    Json::Value params;
    ifstream paramsFile(argv[1]);
    paramsFile >> params;
    cout << params << endl;
    mpz_class a = mpz_class(params["invariants"][0].asString());
    mpz_class b = mpz_class(params["invariants"][1].asString());
    mpz_class fieldOrder = mpz_class(params["fieldOrder"].asString());
    mpz_class curveOrder = mpz_class(params["curveOrder"].asString());
    mpz_class x = mpz_class(params["basePoint"][0].asString());
    mpz_class y = mpz_class(params["basePoint"][1].asString());
    mpz_class z = mpz_class("1");

    init(a, b, fieldOrder, curveOrder);
    Point P = Point(x, y);
    PPoint PP = PPoint(x, y, z);
    cout << P << endl;
    gmp_randclass r(gmp_randinit_default);
    r.seed(getSeed());
    mpz_class scalar = r.get_z_range(curveOrder);
    cout << scalar << endl;
    // Point Q = P * scalar;
    PPoint QQ = PP * scalar;
    // cout << Q << endl;
    cout << QQ << endl;
    auto start = chrono::system_clock::now();
    // CurveField ans = pollard_rho(P,Q);
    CurveField ans = pollard_rho1(PP,QQ);
    auto end = std::chrono::system_clock::now();
    auto t = chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout.fill('0');
    cout << "Answer: " <<  ans.getValue() << endl;
    assert (CurveField(scalar) == ans);
    cout << "Time: "
        << setw(2) << t / (60 * 60 * 1000) << ":"
        << setw(2) << (t / ( 60 * 1000)) % 60 << ":"
        << setw(2) << (t / 1000) % 60 << "."
        << t % 1000
        << endl;
}
