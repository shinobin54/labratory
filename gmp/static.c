#include<stdio.h>
#include<gmp.h>
typedef struct{
  mpz_t x0;
}Fp;
typedef struct{
  Fp x;
  Fp y;
  int infinity;
}EFp;
typedef struct{
  Fp x0;
  Fp x1;
}Fp2;
typedef struct{
  Fp2 x;
  Fp2 y;
  int infinity;
}EFp2;
void Fp_init(Fp *rop);
void Fp_set(Fp *rop, Fp *op);
void Fp_set_ui(Fp *rop, int op);
void Fp_print(Fp *rop);
void Fp_clear(Fp *rop);
void Fp_random(Fp *rop);
void Fp_add(Fp *rop, const Fp *op1, const Fp *op2);
void Fp_sub(Fp *rop, const Fp *op1, const Fp *op2);
void Fp_mul(Fp *rop, const Fp *op1, const Fp *op2);
void Fp_inv(Fp *rop, const Fp *op1);
void Fp_add_ui(Fp *rop, const Fp *op1, int op2);
void Fp_mul_ui (Fp *rop, const Fp *op1, int ui);
void Fp_sub_ui (Fp *rop, const Fp *op1, int ui);
int Fp_cmp (const Fp *op1, const Fp *op2);
int Fp_cmp_ui (const Fp *op1, int op2);
void Fp_neg(Fp *rop,Fp *op);
void Fp_pow_gmp(Fp *rop, const Fp *base, const Fp *exp);
void Fp_pow_ui(Fp *rop, const Fp *base, int exp);
void Fp_div(Fp *rop, const Fp *op1, const Fp *op2);
void Fp_div_ui(Fp *rop, const Fp *op1, int op);
int Fp_legendre(Fp *op);
void Fp_sqrt(Fp *rop,Fp *a);
void Fp_pow(Fp *rop, Fp *op1,mpz_t scalar);
//static mpz_t tmp;
void DBL(mpz_t a);
void init();
int main(){
  init();
  mpz_t a;
  mpz_init(a);
  mpz_set_ui(a,2);
  return 0;
}
void init(){
  static mpz_t tmp;
  mpz_init(tmp);
}
void DBL(mpz_t a){
  static mpz_t tmp;
  mpz_mul_ui(tmp,a,2);
  mpz_set(a,tmp);
}
