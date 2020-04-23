#include <iostream>
#include "affine_point.cpp"
#include <ctime> 
#include "json/json.h"
#include <fstream>

using namespace std;

mpz_t modulo;
mpz_t a;
mpz_t b;
mpz_t curveOrder;


typedef AffinePoint<a, b, modulo> Point;
typedef Field<curveOrder> CurveField;


void init(mpz_class a_v, mpz_class b_v, mpz_class mod, mpz_class curve) {
   mpz_init_set(modulo, mod.get_mpz_t());
   mpz_init_set(a, a_v.get_mpz_t());
   mpz_init_set(b, b_v.get_mpz_t());
   mpz_init_set(curveOrder, curve.get_mpz_t());
}

void step(Point gen, Point y, Point &A, CurveField   &alpha, CurveField &beta){
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

    init(a, b, fieldOrder, curveOrder);
    Point P = Point(x, y);
    P.print();
    Point Q = P * mpz_class(4422);
    Q.print();
    time_t now = time(0);
    CurveField ans = pollard_rho(P,Q);
    time_t end = time(0);
    cout <<difftime(end, now)<<endl;
    cout<<ans.getValue()<<endl;
}
