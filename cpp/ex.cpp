#include <iostream>
#include "affine_point.cpp"
#include <ctime> 

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
    init(755166058714695501, 491852097657725974, 757211378386428139, 757211378543584379);
    Point test = Point(mpz_class(90884636720006406), mpz_class(182740157931376758));
    // init(125054026421, 435854508262, 673996353473, 673995624277);
    // Point test = Point(mpz_class(388974051914), mpz_class(660883941535));
    test.print();

    Point jd = test * mpz_class(4422);
    jd.print();
    time_t now = time(0);
    CurveField ans = pollard_rho(test,jd);
    time_t end = time(0);
    cout << difftime(end, now)<<endl;
    cout<<ans.getValue()<<endl;
}
