#include "gmpxx.h"

template<mpz_t mod>
class Field {
private:
    mpz_class value;
public:
    Field(){}

    Field(mpz_t value){
        mpz_mod(this->value.get_mpz_t(), value, mod);
    }

    Field(mpz_class value){
        mpz_mod(this->value.get_mpz_t(), value.get_mpz_t(), mod);
    }

    Field(int value){
        mpz_t tmp;
        mpz_init_set_si(tmp, value);
        mpz_mod(this->value.get_mpz_t(), tmp, mod);
        mpz_clear(tmp);
    }

    Field(long long value){
        mpz_t tmp;
        mpz_init_set_si(tmp, value);
        mpz_mod(this->value.get_mpz_t(), tmp, mod);
        mpz_clear(tmp);
    }

    void getValue(mpz_t ret){
        mpz_init_set (ret,this->value.get_mpz_t());
    }

    mpz_class getValue(){
        return value;
    }

    Field operator + (Field const &obj) {
        return Field(value + obj.value);
    }

    Field operator - (Field const &obj) {
        return Field(value - obj.value);
    }

    Field operator - () {
        return Field(-value);
    }

    Field operator * (Field const &obj) {
        return Field(value * obj.value);
    }

    Field operator / (Field const &obj) {
        mpz_t tmp;
        mpz_init(tmp);
        mpz_invert(tmp ,obj.value.get_mpz_t(),mod);
        mpz_mul(tmp,value.get_mpz_t(), tmp);
        Field ret = Field(tmp);
        mpz_clear(tmp);
        return ret;
    }

    Field operator % (Field const &obj) {
        return Field(value % obj.value);
    }

    bool operator == (Field const &rhs) {
        return mpz_cmp(value.get_mpz_t(),rhs.value.get_mpz_t()) == 0;
    }

    bool operator != (Field const &rhs) {
        return mpz_cmp(value.get_mpz_t(),rhs.value.get_mpz_t()) != 0;
    }

};
