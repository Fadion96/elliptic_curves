#include "field.cpp"
#include <iostream>


template<mpz_t a, mpz_t b, mpz_t mod>
class AffinePoint
{
private:
    typedef Field<mod> Number;
    Number x;
    Number y;
public:
    AffinePoint(mpz_t x, mpz_t y){
       this->x = Number(x);
       this->y = Number(y);
    }
    AffinePoint(mpz_class x, mpz_class y){
       this->x = Number(x);
       this->y = Number(y);
    }
    AffinePoint(long long x, long long y){
       this->x = Number(x);
       this->y = Number(y);
    }
    AffinePoint(int x, int y){
       this->x = Number(x);
       this->y = Number(y);
    }
    AffinePoint(Number x, Number y){
       this->x = x;
       this->y = y;
    }

    bool is_zero(){
        return x == Number(0);
    }

    Number getX(){
        return x;
    }

    Number getY(){
        return y;
    }


    AffinePoint dbl() {
        if (this->is_zero() || this->y == Number(0)){
            return AffinePoint(0,0);
        }
        else{ 
            Number temp = (x * x  * Number(3) + Number(a)) / (y * Number(2));
            Number rx = temp * temp - x * Number(2);
            Number ry = temp * (x - rx) - y;
            return AffinePoint(rx, ry);
        }
    }

    AffinePoint operator+ (AffinePoint const &other) {
        if(this->is_zero()) {
            return other;
        }
        else if (Number(other.x) == Number(0)) {
            return *this;
        }
        else if (this->x == other.x) {
            if (this->y == other.y){
                return this->dbl();
            }
            else{
                return AffinePoint(0,0);
            }
        }
        else {
            Number temp = (Number(other.y) - y) / (Number(other.x) - x);
            Number rx = temp * temp - x - other.x;
            Number ry = temp * (x - rx) - y;
            return AffinePoint(rx, ry);
        }
    }

    AffinePoint operator* (mpz_class n){
        AffinePoint result = AffinePoint(0,0);
        AffinePoint temp = *this;
        while(n != 0){
            if (n % 2 != 0){
                result = result + temp;
            }
            temp =  temp.dbl();
            n = n >> 1;
        }
        return result;
    }

    bool operator== (AffinePoint const &other){
         if(this->is_zero() || Number(other.x) == Number(0)){
             return this->is_zero() && Number(other.x) == Number(0);
         }
         else {
             return (x == other.x) && (y == other.y);
         }
    }
    
    void print() {
        std::cout << "(" << x.getValue()<< ", " << y.getValue() <<")"<<std::endl; // actual output done here
 
    }
     
};

