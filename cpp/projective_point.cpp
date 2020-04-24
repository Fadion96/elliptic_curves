// #include "field.cpp"
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
        return z == Number(0);
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
        if (this->is_zero() || this->y == Number(1)){
            return ProjectivePoint(0, 0, 0);
        }
        else{
            Number two = Number(2);
            Number w = x * x * Number(3) + Number(a) * z * z;
            Number s = y * z * two;
            Number B = s * x * y * two;
            Number h = w * w - B * two;
            Number rx = h * s ;
            Number ry = w * (B - h) - s * s * y * y * two;
            Number rz = s * s * s;
            return ProjectivePoint(rx, ry, rz);
        }
    }

    ProjectivePoint operator + (ProjectivePoint const &other) {
        if(this->is_zero()) {
            return other;
        }
        else {
            if (Number(other.z) == Number(0)) {
                return *this;
            }
        }
        Number s_1 = this->y * other.z;
        Number s_2 = this->z * other.y;
        Number u_1 = this->x * other.z;
        Number u_2 = this->z * other.x;
        if (u_1 == u_2) {
            if (s_1 == s_2) {
                return this->dbl();
            } else {
                return ProjectivePoint(0, 1, 0);
            }
        } else {
            Number w = z * other.z;
            Number r = s_2 - s_1; 
            Number p = u_2 - u_1;   
            Number p_sqr = p * p;
            Number rrw = r * r * w;
            Number p_sqr_u1_u2 = p_sqr * (u_1 + u_2);
            Number p_cube = p * p_sqr;
            Number rx = p * (rrw - p_sqr_u1_u2);
            Number ry = r * (Number(-2)*rrw + Number(3) * p_sqr_u1_u2) - p_cube * (s_1 + s_2);
            Number rz = p_cube * w;
            return ProjectivePoint(Number(2) * rx, ry, Number(2) * rz);
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

