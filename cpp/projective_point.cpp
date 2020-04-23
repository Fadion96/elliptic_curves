//#include "field.cpp"
#include <iostream>




template<mpz_t a, mpz_t b, mpz_t mod>
class ProjectivePoint
{
private:
    typedef Field<mod> Number;
    Number x;
    Number y;
    Number z;
public:
    ProjectivePoint(mpz_t x, mpz_t y, mpz_t z){
       this->x = Number(x);
       this->y = Number(y);
       this->z = Number(z);
    }
    ProjectivePoint(mpz_class x, mpz_class y, mpz_class z){
       this->x = Number(x);
       this->y = Number(y);
       this->z = Number(z);
    }
    ProjectivePoint(long long x, long long y, long long z){
       this->x = Number(x);
       this->y = Number(y);
       this->z = Number(z);
    }
    ProjectivePoint(int x, int y, int z){
       this->x = Number(x);
       this->y = Number(y);
       this->z = Number(z);
    }
    ProjectivePoint(Number x, Number y ,Number z){
       this->x = x;
       this->y = y;
       this->z = z;
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

    Number getZ(){
        return z;
    }

    ProjectivePoint dbl() {
        if (this->is_zero() || this->y == Number(0)){
            return ProjectivePoint(0, 0, 0);
        }
        else{
            Number two = Number(2);
            Number t = this->x * this->x * Number(3) + Number(a) * this->z * this->z;
            Number u = this->y * this->z * two;
            Number v = u * this->x * this->y * two;
            Number w = t * t - v * two;
            Number rx = u * w;
            Number ry = t * (v - w) - u * u * this->y * this->y * two;
            Number rz = u * u * u;
            return ProjectivePoint(rx, ry, rz);
        }
    }

    ProjectivePoint operator + (ProjectivePoint const &other) {
        if(this->is_zero()) {
            return other;
        }
        else {
            if (Number(other.x) == Number(0)) {
                return *this;
            }
        }
        Number t0 = this->y * other.z;
        Number t1 = this->z * other.y;
        Number u0 = this->x * other.z;
        Number u1 = this->z * other.x;
        if (u0 == u1) {
            if (t0 == t1) {
                return this->dbl();
            } else {
                return ProjectivePoint(0, 0, 0);
            }
        } else {
            Number t = t0 - t1;
            Number u = u0 - u1;
            Number u2 = u * u;
            Number v = this->z * other.z;
            Number w = t * t * v - u2 * (u0 + u1);
            Number u3 = u * u2;
            Number rx = u * w;
            Number ry = t * (u0 * u2 - w) - t0 * u3;
            Number rz = u3 * v;
            return ProjectivePoint(rx, ry, rz);
        }
    }

    ProjectivePoint operator - (ProjectivePoint const &other) {
        return *this + -other;
    }

    ProjectivePoint operator - () {
        if (this->is_zero()) {
            return *this;
        } else {
            return ProjectivePoint(this->x, -this->y, this->z);
        }
    }

    ProjectivePoint operator * (mpz_class n){
        if (n < 0) {
            return -(*this) * -n;
        }
        ProjectivePoint result = ProjectivePoint(0, 0, 0);
        ProjectivePoint temp = *this;
        while(n != 0){
            if (n % 2 != 0){
                result = result + temp;
            }
            temp =  temp.dbl();
            n = n >> 1;
        }
        return result;
    }

    bool operator != (ProjectivePoint const &other) {
        return !(this == other);
    }

    bool operator == (ProjectivePoint const &other){
         if(this->is_zero() || Number(other.x) == Number(0)){
             return this->is_zero() && Number(other.x) == Number(0);
         }
         else {
             return (this->x * other.z == this->z * other.x ) && (this->y * other.z == this->z * other.y );
         }
    }

    friend std::ostream & operator << (std::ostream &out, ProjectivePoint &point){
        out << "(" << point.x.getValue()<< ", " << point.y.getValue() << ", " << point.z.getValue() <<")";
        return out;
    }

};

