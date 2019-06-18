#include<stdio.h>
#include<gmp.h>
#include<time.h>
#include<unistd.h>
#include<math.h>
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
gmp_randstate_t state;
/*Fp,EFp*/
mpz_t Fp_p,SS_r,EFp_s,Fp_x;
Fp Fp_unity,Fp_zero;
Fp Fp_b,Fp_a;
mpz_t mpz_tmp1,mpz_tmp2,mpz_tmp3;
Fp Fp_tmp1,Fp_tmp2,Fp_tmp3;
EFp EFp_tmp1,EFp_tmp2,EFp_tmp3,EFp_tmp4;
/*Fp2,EFp2*/
Fp Fp2_i;
mpz_t Fp2_p,EFp2_s;
Fp2 Fp2_unity,Fp2_zero,Fp2_b,Fp2_a;
Fp2 Fp2_tmp1,Fp2_tmp2,Fp2_tmp3,Fp2_tmp4,Fp2_tmp5,Fp2_tmp6,Fp2_tmp7;
EFp2 EFp2_tmp1,EFp2_tmp2,EFp2_tmp3,EFp2_tmp4;

EFp2 P_ex,Q_ex;
void billinear(EFp2 *op1,EFp2 *op2);
/*Fp*/
void Fp_testinit();
void Fp_make();
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
/*EFp*/
void EFp_random(EFp *rop);
void EFp_DBL(EFp *rop, EFp *op);
void EFp_ADD(EFp *rop, EFp *op1, EFp *op2);
void EFp_SCM(EFp *rop, EFp *op, mpz_t s);
void EFp_print(EFp *rop);
void EFp_init(EFp *rop);
void EFp_set(EFp *rop, EFp *op);
void EFp_clear(EFp *rop);
int EFp_cmp(EFp *op1, EFp *op2);
/*Fp2*/
void Fp2_testinit();
void Fp2_make();
void Fp2_init(Fp2 *op);
int Fp2_cmp(Fp2 *op1,Fp2 *op2);
void Fp2_random(Fp2 *rop);
void Fp2_clear(Fp2 *op);
void Fp2_sub(Fp2 *rop, Fp2 *op1, Fp2 *op2);
void Fp2_add(Fp2 *rop, Fp2 *op1, Fp2 *op2);
void Fp2_mul(Fp2 *rop, Fp2 *op1, Fp2 *op2);
void Fp2_conju(Fp2 *rop,Fp2 *op);
void Fp2_inv(Fp2 *rop, Fp2 *op);
void Fp2_print(Fp2 *rop);
void Fp2_set(Fp2 *rop,Fp2 *op);
void Fp2_pow(Fp2 *rop, Fp2 *op1, mpz_t scalar);
void Fp2_pow_ui(Fp2 *rop,Fp2 *op1,unsigned int op2);
void Fp2_div(Fp2 *rop,Fp2 *op1,Fp2 *op2);
void Fp2_sqrt(Fp2 *rop,Fp2 *op);
void Fp2_mul_ui(Fp2 *rop,Fp2 *op1,unsigned int op2);
int Fp2_legendre(Fp2 *op);
int Fp2_legendre3(Fp2 *op);
/*EFp2*/
void EFp2_random(EFp2 *rop);
void EFp2_DBL(EFp2 *rop, EFp2 *op);
void EFp2_ADD(EFp2 *rop, EFp2 *op1, EFp2 *op2);
void EFp2_SCM(EFp2 *rop, EFp2 *op, mpz_t s);
void EFp2_print(EFp2 *rop);
void EFp2_init(EFp2 *rop);
void EFp2_set(EFp2 *rop, EFp2 *op);
void EFp2_clear(EFp2 *rop);
int EFp2_cmp(EFp2 *op1, EFp2 *op2);
void EFp2_wail();

void EFp2_lineTP(Fp2 *rop,EFp2 *q,EFp2 *p,EFp2 *t);
void EFp2_lineTT(Fp2 *rop,EFp2 *q,EFp2 *t);
void Fp2_Miller(Fp2 *rop,mpz_t s,EFp2 *p,EFp2 *q);
void Fp2paring(Fp2 *rop,EFp2 *p,EFp2 *q);
void EFp2_set_P(EFp2 *P);
void EFp2_set_Q(EFp2 *Q);
void inversiontest();
void paringtest();
void EFp2_set_gQ(EFp2 *Q);
void EFp2_set_gP(EFp2 *P);

int main(){
  Fp_testinit();
  Fp_make();
  Fp2_testinit();
  Fp2_make();
  //paringtest();

  //paringtest();
  inversiontest();
  printf("ANS=(109,25)\n");

  Fp2 Fp2_1,Fp2_2,Fp2_3;
  EFp EFp_1,EFp_2,EFp_3;
  EFp2 EFp2_1,EFp2_2,EFp2_3;
  int i;
  mpz_t mpz_s;
  Fp2_init(&Fp2_1);
  Fp2_init(&Fp2_2);
  Fp2_init(&Fp2_3);
  EFp2_init(&EFp2_1);
  EFp2_init(&EFp2_2);
  EFp2_init(&EFp2_3);
  EFp_init(&EFp_1);
  EFp_init(&EFp_2);
  EFp_init(&EFp_3);
  mpz_init(mpz_s);
  mpz_set_ui(mpz_s,35);
/*
  EFp2_random(&EFp2_1);
  EFp2_random(&EFp2_3);
  EFp2_ADD(&EFp2_2,&EFp2_1,&EFp2_3);
  EFp2_print(&EFp2_3);

  EFp2_ADD(&EFp2_1,&EFp2_1,&EFp2_3);
  EFp2_print(&EFp2_3);
*/
  /*
  EFp_random(&EFp_1);
  EFp_SCM(&EFp_1,&EFp_1,EFp_s);
  printf("aaa=\n");
  EFp_print(&EFp_1);
  EFp2_random(&EFp2_1);
  EFp2_SCM(&EFp2_2,&EFp2_1,SS_r);
  printf("bbb=\n");
  EFp2_print(&EFp2_2);
*/
  /*
  mpz_set_ui(mpz_s,4);
  EFp2_random(&EFp2_1);
  EFp2_2.infinity=1;
  EFp2_set(&EFp2_3,&EFp2_1);

  for(i=0;i<4;i++){
    EFp2_ADD(&EFp2_2,&EFp2_2,&EFp2_1);
  }
  EFp2_print(&EFp2_2);
  for(i=0;i<2;i++){
    EFp2_DBL(&EFp2_3,&EFp2_3);
  }
  EFp2_print(&EFp2_3);
  EFp2_SCM(&EFp2_1,&EFp2_1,mpz_s);
  EFp2_print(&EFp2_1);
  //EFp2_random(&EFp2_2);
  if(EFp2_cmp(&EFp2_1,&EFp2_2)==0 && EFp2_cmp(&EFp2_1,&EFp2_3)==0){
    printf("p*4=p+p+p+p=p*2*2:success!!\n");
  }else{
    printf("false\n");
  }
  */
  /*
  EFp2 a,b,c;
  mpz_t s;
  mpz_init(s);
  EFp2_init(&a);
  EFp2_init(&b);

  EFp2_init(&c);
  EFp2_random(&a);
  EFp2_SCM(&a,&a,EFp2_s);
  EFp2_print(&a);
  */
  /*
  EFp2 p,q;
  EFp2_init(&p);
  EFp2_init(&q);
  EFp2_set_P(&p);
  EFp2_set_Q(&q);
  EFp2_print(&p);
  EFp2_print(&q);
  inversiontest();
  */
  /*
  EFp a,b,c;
  EFp_init(&a);
  EFp_init(&b);
  EFp_init(&c);

  EFp_random(&a);
  EFp_DBL(&b,&a);
  EFp_DBL(&b,&b);
  EFp_print(&b);

  EFp_DBL(&b,&a);
  EFp_set(&c,&b);
  EFp_ADD(&b,&b,&c);
  EFp_print(&b);

  EFp_ADD(&b,&a,&a);
  EFp_ADD(&c,&b,&a);
  EFp_ADD(&b,&c,&a);
  EFp_print(&b);
  */
  /*
  EFp_DBL(&b,&a);
  EFp_set(&c,&b);
  EFp_ADD(&b,&b,&c);
  EFp_print(&b);
*/

  return 0;
}

/*Fp*/
void Fp_testinit(){
  mpz_init(Fp_p);
  mpz_init(SS_r);
  mpz_init(EFp_s);
  mpz_init(Fp_x);
  mpz_init(mpz_tmp1);
  mpz_init(mpz_tmp2);
  mpz_init(mpz_tmp3);
  Fp_init(&Fp_unity);
  Fp_init(&Fp_zero);
  Fp_init(&Fp_b);
  Fp_init(&Fp_a);
  Fp_init(&Fp_tmp1);
  Fp_init(&Fp_tmp2);
  Fp_init(&Fp_tmp3);
  EFp_init(&EFp_tmp1);
  EFp_init(&EFp_tmp2);
  EFp_init(&EFp_tmp3);
  EFp_init(&EFp_tmp4);

  gmp_randinit_mt(state);
  gmp_randseed_ui(state, (unsigned int) time(NULL));
}
void Fp_make(){
  mpz_set_ui(Fp_p,139);
  mpz_set_ui(SS_r,140);
  mpz_add_ui(EFp_s,Fp_p,1);
  Fp_set_ui(&Fp_unity,1);
  Fp_sub_ui(&Fp_a,&Fp_zero,13);
  Fp_sub_ui(&Fp_b,&Fp_zero,7);
}
void Fp_random(Fp *rop){
  mpz_urandomm(rop->x0, state, Fp_p);
}
void Fp_init(Fp *rop){
  mpz_init(rop->x0);
}
void Fp_set_ui(Fp *rop, int op){
  mpz_set_ui(rop->x0,op);
}
void Fp_set(Fp *rop, Fp *op){
  mpz_set(rop->x0, op->x0);
}
void Fp_clear(Fp *rop){
  mpz_clear(rop->x0);
}
void Fp_add(Fp *rop, const Fp *op1, const Fp *op2){ //Fp_pを法とした足し算の定義
  mpz_add(rop->x0, op1->x0, op2->x0);
  mpz_mod(rop->x0, rop->x0, Fp_p);
}
void Fp_add_ui(Fp *rop, const Fp *op1, int op2){ //moduleを法とした足し算の定義
  mpz_add_ui(rop->x0, op1->x0, op2);
  mpz_mod(rop->x0, rop->x0, Fp_p);
}
void Fp_sub(Fp *rop, const Fp *op1, const Fp *op2){
  mpz_sub(rop->x0, op1->x0, op2->x0);
  mpz_mod(rop->x0, rop->x0, Fp_p);
}
void Fp_sub_ui (Fp *rop, const Fp *op1, int ui){
  mpz_sub_ui (rop->x0, op1->x0, ui);
  mpz_mod(rop->x0, rop->x0, Fp_p);
}
void Fp_mul(Fp *rop, const Fp *op1, const Fp *op2){ //moduleを法とした掛け算の定義
  mpz_mul(rop->x0, op1->x0, op2->x0);
  mpz_mod(rop->x0, rop->x0, Fp_p);
}
void Fp_mul_ui (Fp *rop, const Fp *op1, int ui){
  mpz_mul_ui (rop->x0, op1->x0, ui);
  mpz_mod(rop->x0, rop->x0, Fp_p);
}
void Fp_inv(Fp *rop, const Fp *op1){
  mpz_invert(rop->x0, op1->x0, Fp_p);
}
void Fp_neg(Fp *rop,Fp *op){
  mpz_neg(rop->x0,op->x0);
  mpz_mod(rop->x0, rop->x0, Fp_p);
}
void Fp_pow(Fp *rop, Fp *op1, mpz_t scalar){//Fp_tmp1
  int i;
  Fp_set_ui(&Fp_tmp1,1);
  Fp_set(&Fp_tmp2,op1);
  for(i=(int)(mpz_sizeinbase(scalar,2)-1);i>=0;i--){
    Fp_mul(&Fp_tmp1,&Fp_tmp1,&Fp_tmp1);
    if(mpz_tstbit(scalar, i)==1){
      Fp_mul(&Fp_tmp1,&Fp_tmp1,&Fp_tmp2);
    }
  }
  Fp_set(rop,&Fp_tmp1);
}
int Fp_legendre(Fp *op){//Fp_tmp1,Fp_tmp2,mpz_tmp1
  mpz_sub_ui(mpz_tmp1,Fp_p,1);
  mpz_divexact_ui(mpz_tmp1,mpz_tmp1,2);
  Fp_pow(&Fp_tmp2,op,mpz_tmp1);
  if(Fp_cmp(&Fp_tmp1,&Fp_unity)==0||Fp_cmp(&Fp_tmp1,&Fp_zero)==0){
    return 1;
  }else{
    return -1;
  }
}
void Fp_sqrt(Fp *rop, Fp *a){
  mpz_t q,mpz_temp;
  Fp n,x,y,b,t,Fp_temp,Fp_temp2;
  unsigned int e,r,i,m;
  Fp_init(&n);
  Fp_init(&x);
  Fp_init(&b);
  Fp_init(&y);
  Fp_init(&t);
  Fp_init(&Fp_temp);
  Fp_init(&Fp_temp2);

  mpz_init(q);
  mpz_init(mpz_temp);

  //1
  while(1){
    Fp_random(&n);
    if(Fp_legendre(&n)==-1){
      break;
    }
  }
  //2
  mpz_sub_ui(mpz_temp,Fp_p,1);
  e = mpz_scan1(mpz_temp,0);
  mpz_fdiv_q_2exp(q,mpz_temp,e);
  //3
  Fp_pow(&y,&n,q);
  r=e;
  mpz_sub_ui(mpz_temp,q,1);
  mpz_divexact_ui(mpz_temp,mpz_temp,2);
  Fp_pow(&x,a,mpz_temp);
  //4
  Fp_pow_ui(&b,&x,2);
  Fp_mul(&b,&b,a);
  Fp_mul(&x,a,&x);
  //5
  while(Fp_cmp_ui(&b,1)!=0){
    m=0;
    Fp_set(&Fp_temp,&b);
    while(Fp_cmp_ui(&Fp_temp,1)!=0){
      Fp_pow_ui(&Fp_temp,&Fp_temp,2);
      m++;
    }
    r=r-m-1;
    mpz_set_ui(mpz_temp,2);
    mpz_pow_ui(mpz_temp,mpz_temp,r);
    Fp_pow(&t,&y,mpz_temp);
    //ok
    Fp_pow_ui(&y,&t,2);
    r=m;
    Fp_mul(&x,&x,&t);
    Fp_mul(&b,&b,&y);
  }
  Fp_set(rop,&x);
  mpz_clear(mpz_temp);
  Fp_clear(&Fp_temp2);
  Fp_clear(&Fp_temp);
  Fp_clear(&n);
  mpz_clear(q);
  Fp_clear(&x);
  Fp_clear(&b);
  Fp_clear(&y);
  Fp_clear(&t);
}
int Fp_cmp (const Fp *op1, const Fp *op2){
  return mpz_cmp (op1->x0, op2->x0);
}
int Fp_cmp_ui (const Fp *op1, int op2){
  return mpz_cmp_ui (op1->x0, op2);
}
void Fp_pow_ui(Fp *rop, const Fp *base, int exp){
  mpz_powm_ui(rop->x0, base->x0, exp, Fp_p);
}
void Fp_div(Fp *rop, const Fp *op1, const Fp *op2){//Fp_tmp1
  Fp_inv(&Fp_tmp1, op2);
  Fp_mul(rop, op1, &Fp_tmp1);
}
void Fp_div_ui(Fp *rop, const Fp *op1, int op){//Fp_tmp1
  Fp_set_ui(&Fp_tmp2,op);
  Fp_inv(&Fp_tmp1, &Fp_tmp2);
  Fp_mul(rop, op1, &Fp_tmp1);
}
void Fp_print(Fp *rop){
  mpz_out_str(stdout,10,rop->x0);
  printf("\n");
}

/*EFp*/
void EFp_random(EFp *rop){/*関数内定義*/
  Fp tmp1,tmp2;
  Fp_init(&tmp1);
  Fp_init(&tmp2);
  /*tmp1=x^3+ax+b*/
  while(1){
    Fp_random(&rop->x);
    Fp_mul(&tmp1,&rop->x,&rop->x);
    Fp_mul(&tmp1,&tmp1,&rop->x);
    Fp_mul(&tmp2,&rop->x,&Fp_a);
    Fp_add(&tmp1,&tmp1,&tmp2);
    Fp_add(&tmp1,&tmp1,&Fp_b);
    if(Fp_legendre(&tmp1)==1) break;
  }
  Fp_sqrt(&rop->y,&tmp1);
  rop->infinity=0;
  Fp_clear(&tmp1);
  Fp_clear(&tmp2);
}
void EFp_DBL(EFp *rop, EFp *op){
  if(op->infinity==1){
    rop->infinity = 1;
  }else if(Fp_cmp_ui(&(op->y),0)==0){
    rop->infinity = 1;
  }else{
    EFp_set(&EFp_tmp1,op);

    Fp_add(&Fp_tmp1,&EFp_tmp1.y,&EFp_tmp1.y);
    Fp_inv(&Fp_tmp1,&Fp_tmp1);

    Fp_mul(&Fp_tmp2,&EFp_tmp1.x,&EFp_tmp1.x);
    Fp_add(&Fp_tmp3,&Fp_tmp2,&Fp_tmp2);
    Fp_add(&Fp_tmp2,&Fp_tmp2,&Fp_tmp3);
    Fp_add(&Fp_tmp2,&Fp_tmp2,&Fp_a);

    Fp_mul(&Fp_tmp3,&Fp_tmp1,&Fp_tmp2);
    Fp_mul(&Fp_tmp1,&Fp_tmp3,&Fp_tmp3);

    Fp_add(&Fp_tmp2,&EFp_tmp1.x,&EFp_tmp1.x);
    Fp_sub(&rop->x,&Fp_tmp1,&Fp_tmp2);

    Fp_sub(&Fp_tmp1,&EFp_tmp1.x,&rop->x);
    Fp_mul(&Fp_tmp2,&Fp_tmp3,&Fp_tmp1);
    Fp_sub(&rop->y,&Fp_tmp2,&EFp_tmp1.y);
  }
}
void EFp_ADD(EFp *rop, EFp *op1, EFp *op2){
  if(op1->infinity==1){
      EFp_set(rop,op2);
      return;
  }else if(op2->infinity==1){
      EFp_set(rop,op1);
      return;
  }else if(Fp_cmp(&op1->x,&op2->x)==0){
      if(Fp_cmp(&op1->y,&op2->y)!=0){
          rop->infinity=1;
          return;
      }else{
          EFp_DBL(rop,op1);
          return;
      }
  }
  EFp_set(&EFp_tmp1,op1);
  EFp_set(&EFp_tmp2,op2);

  Fp_sub(&Fp_tmp1,&EFp_tmp2.x,&EFp_tmp1.x);
  Fp_inv(&Fp_tmp1,&Fp_tmp1);
  Fp_sub(&Fp_tmp2,&EFp_tmp2.y,&EFp_tmp1.y);
  Fp_mul(&Fp_tmp3,&Fp_tmp1,&Fp_tmp2);
  Fp_mul(&Fp_tmp1,&Fp_tmp3,&Fp_tmp3);

  Fp_sub(&Fp_tmp2,&Fp_tmp1,&EFp_tmp1.x);
  Fp_sub(&rop->x,&Fp_tmp2,&EFp_tmp2.x);

  Fp_sub(&Fp_tmp1,&EFp_tmp1.x,&rop->x);
  Fp_mul(&Fp_tmp2,&Fp_tmp3,&Fp_tmp1);
  Fp_sub(&rop->y,&Fp_tmp2,&EFp_tmp1.y);
}
void EFp_SCM(EFp *rop, EFp *op, mpz_t scalar){
  if(EFp_tmp3.infinity==1){
    rop->infinity=1;
    return;
  }
  int i;
  EFp_set(&EFp_tmp3,op);
  EFp_tmp4.infinity=1;
  for(i=(int)(mpz_sizeinbase(scalar,2)-1);i>=0;i--){
    EFp_DBL(&EFp_tmp4,&EFp_tmp4);
    if(mpz_tstbit(scalar, i)==1){
      EFp_ADD(&EFp_tmp4,&EFp_tmp4,&EFp_tmp3);
    }
  }
  EFp_set(rop,&EFp_tmp4);
}
int EFp_cmp(EFp *op1, EFp *op2){
  if(Fp_cmp(&op1->x,&op2->x)==0 && Fp_cmp(&op1->y,&op2->y)==0){
    return 0;
  }else if(op1->infinity==1&&op2->infinity==1){
    return 0;
  }else{
  return 1;
  }
}
void EFp_init(EFp *rop){
  Fp_init(&(rop->x));
  Fp_init(&(rop->y));
  rop->infinity = 0;
}
void EFp_clear(EFp *rop){
  Fp_clear(&(rop->x));
  Fp_clear(&(rop->y));
}
void EFp_print(EFp *rop){
  if(rop->infinity==0){
    printf("(");
    mpz_out_str(stdout,10,rop->x.x0);
    printf(",");
    mpz_out_str(stdout,10,rop->y.x0);
    printf(")\n");
  }else{
    printf("0\n");
  }
}
void EFp_set(EFp *rop, EFp *op){
  Fp_set(&(rop->x), &(op->x));
  Fp_set(&(rop->y), &(op->y));
  rop->infinity = op->infinity;
}

/*Fp2*/
void Fp2_testinit(){
  Fp_init(&Fp2_i);
  mpz_init(Fp2_p);
  mpz_init(EFp2_s);

  Fp2_init(&Fp2_unity);
  Fp2_init(&Fp2_zero);

  Fp2_init(&Fp2_b);
  Fp2_init(&Fp2_a);

  Fp2_init(&Fp2_tmp1);
  Fp2_init(&Fp2_tmp2);
  Fp2_init(&Fp2_tmp3);
  Fp2_init(&Fp2_tmp4);
  Fp2_init(&Fp2_tmp5);
  Fp2_init(&Fp2_tmp6);
  Fp2_init(&Fp2_tmp7);

  EFp2_init(&EFp2_tmp1);
  EFp2_init(&EFp2_tmp2);
  EFp2_init(&EFp2_tmp3);
  EFp2_init(&EFp2_tmp4);
}
void Fp2_make(){
  Fp_sub_ui(&Fp2_i,&Fp_zero,4);
  mpz_mul(Fp2_p,Fp_p,Fp_p);
  mpz_mul(EFp2_s,EFp_s,EFp_s);
  Fp_set_ui(&Fp2_unity.x0,1);
  Fp_set(&Fp2_a.x0,&Fp_a);
  Fp_set(&Fp2_b.x0,&Fp_b);

  EFp2_init(&P_ex);

  EFp2_init(&Q_ex);

  Fp_set_ui(&P_ex.x.x0,67);
  Fp_set_ui(&P_ex.x.x1,38);
  /*Q=(59,-54)*/
  Fp_set_ui(&Q_ex.x.x0,59);
  Fp_sub_ui(&Q_ex.y.x0,&Fp_zero,54);
}
int Fp2_legendre(Fp2 *op){
  mpz_sub_ui(mpz_tmp1,Fp2_p,1);
  mpz_divexact_ui(mpz_tmp1,mpz_tmp1,2);
  Fp2_pow(&Fp2_tmp1,op,mpz_tmp1);
  if(Fp2_cmp(&Fp2_tmp1,&Fp2_unity)==0||Fp2_cmp(&Fp2_tmp1,&Fp2_zero)==0){
    return 1;
  }else{
    return -1;
  }
}
void Fp2_pow_ui(Fp2 *rop,Fp2 *op1,unsigned int op2){
  mpz_set_ui(mpz_tmp1,op2);
  Fp2_pow(rop,op1,mpz_tmp1);
}
void Fp2_mul_ui(Fp2 *rop,Fp2 *op1,unsigned int op2){/*スカラー倍算*/
  Fp_mul_ui(&(rop->x0),&(op1->x0),op2);
  Fp_mul_ui(&(rop->x1),&(op1->x1),op2);
}
void Fp2_sqrt(Fp2 *rop, Fp2 *a){
  mpz_t q,mpz_temp,mpz_temp2;
  Fp2 n,x,y,b,t,Fp2_temp,Fp2_temp2;
  unsigned int e,r,i,m;
  Fp2_init(&n);
  Fp2_init(&x);
  Fp2_init(&b);
  Fp2_init(&y);
  Fp2_init(&t);
  Fp2_init(&Fp2_temp);
  Fp2_init(&Fp2_temp2);

  mpz_init(q);
  mpz_init(mpz_temp);
  mpz_init(mpz_temp2);

  //1
  while(1){
    Fp2_random(&n);
    if(Fp2_legendre(&n)==-1){
      break;
    }
  }
  //2
  mpz_sub_ui(mpz_temp,Fp2_p,1);
  e = mpz_scan1(mpz_temp,0);
  mpz_fdiv_q_2exp(q,mpz_temp,e);
  //3
  Fp2_pow(&y,&n,q);
  r=e;
  mpz_sub_ui(mpz_temp,q,1);
  mpz_divexact_ui(mpz_temp,mpz_temp,2);
  Fp2_pow(&x,a,mpz_temp);
  //4
  Fp2_pow_ui(&b,&x,2);
  Fp2_mul(&b,&b,a);
  Fp2_mul(&x,a,&x);
  //5
  while(Fp2_cmp(&b,&Fp2_unity)!=0){
    m=0;
    Fp2_set(&Fp2_temp,&b);
    while(Fp2_cmp(&Fp2_temp,&Fp2_unity)!=0){
      Fp2_pow_ui(&Fp2_temp,&Fp2_temp,2);
      m++;
    }
    r=r-m-1;
    mpz_set_ui(mpz_temp,2);
    mpz_pow_ui(mpz_temp,mpz_temp,r);
    Fp2_pow(&t,&y,mpz_temp);
    //ok
    Fp2_pow_ui(&y,&t,2);
    r=m;
    Fp2_mul(&x,&x,&t);
    Fp2_mul(&b,&b,&y);
  }
  Fp2_set(rop,&x);
  mpz_clear(mpz_temp);
  Fp2_clear(&Fp2_temp2);
  Fp2_clear(&Fp2_temp);
  Fp2_clear(&n);
  mpz_clear(q);
  Fp2_clear(&x);
  Fp2_clear(&b);
  Fp2_clear(&y);
  Fp2_clear(&t);
}
void Fp2_pow(Fp2 *rop, Fp2 *op1, mpz_t scalar){
  int i;
  Fp2_set(&Fp2_tmp4,op1);
  Fp2_set(&Fp2_tmp3,&Fp2_unity);
  for(i=(int)(mpz_sizeinbase(scalar,2)-1);i>=0;i--){
    Fp2_mul(&Fp2_tmp3,&Fp2_tmp3,&Fp2_tmp3);
    if(mpz_tstbit(scalar, i)==1){
      Fp2_mul(&Fp2_tmp3,&Fp2_tmp3,&Fp2_tmp4);
    }
  }
  Fp2_set(rop,&Fp2_tmp3);
}
void Fp2_init(Fp2 *op){
  Fp_init(&(op->x0));
  Fp_init(&(op->x1));
}
void Fp2_clear(Fp2 *op){
  Fp_clear(&(op->x0));
  Fp_clear(&(op->x1));
}
void Fp2_random(Fp2 *rop){
  Fp_random(&(rop->x0));
  Fp_random(&(rop->x1));
}
void Fp2_add(Fp2 *rop, Fp2 *op1, Fp2 *op2){
  Fp_add(&(rop->x0),&(op1->x0),&(op2->x0));
  Fp_add(&(rop->x1),&(op1->x1),&(op2->x1));
}
void Fp2_sub(Fp2 *rop, Fp2 *op1, Fp2 *op2){
  Fp_sub(&(rop->x0),&(op1->x0),&(op2->x0));
  Fp_sub(&(rop->x1),&(op1->x1),&(op2->x1));
}
void Fp2_mul(Fp2 *rop, Fp2 *op1, Fp2 *op2){
  /*use Fp2_tmp1,Fp2_tmp2*/
  Fp_mul(&(Fp2_tmp1.x0),&(op1->x0),&(op2->x0));
  Fp_mul(&(Fp2_tmp2.x0),&(op1->x1),&(op2->x1));
  Fp_mul(&(Fp2_tmp2.x0),&(Fp2_tmp2.x0),&Fp2_i);

  Fp_mul(&(Fp2_tmp1.x1),&(op1->x0),&(op2->x1));
  Fp_mul(&(Fp2_tmp2.x1),&(op1->x1),&(op2->x0));

  Fp2_add(rop,&Fp2_tmp1,&Fp2_tmp2);
}
void Fp2_conju(Fp2 *rop,Fp2 *op){
  Fp_set(&(rop->x0),&(op->x0));
  Fp_neg(&(rop->x1),&(op->x1));
}
void Fp2_neg(Fp2 *rop,Fp2 *op){
  Fp2_sub(rop,&Fp2_zero,op);
}
void Fp2_inv(Fp2 *rop, Fp2 *op){
  /*use Fp2_tmp3,Fp2_tmp4,*/
  Fp2_conju(&Fp2_tmp3,op);
  //Fp2_tmp4=(a+bi)*(a-bi)
  Fp2_mul(&Fp2_tmp4,op,&Fp2_tmp3);
  Fp_inv(&Fp2_tmp4.x0,&Fp2_tmp4.x0);
  Fp_mul(&rop->x0,&Fp2_tmp3.x0,&Fp2_tmp4.x0);
  Fp_mul(&rop->x1,&Fp2_tmp3.x1,&Fp2_tmp4.x0);
}
void Fp2_div(Fp2 *rop,Fp2 *op1,Fp2 *op2){
  Fp2 inv;
  Fp2_init(&inv);
  Fp2_inv(&inv,op2);
  Fp2_mul(rop,op1,&inv);
  Fp2_clear(&inv);
}
int Fp2_cmp(Fp2 *op1,Fp2 *op2){
  if(Fp_cmp(&(op1->x0),&(op2->x0))==0 && Fp_cmp(&(op1->x1),&(op2->x1))==0){
    return 0;
  }else{
    return 1;
  }
}
void Fp2_print(Fp2 *rop){
  printf("(");
  mpz_out_str(stdout,10,rop->x0.x0);
  printf(",");
  mpz_out_str(stdout,10,rop->x1.x0);
  printf(")\n");
}
void Fp2_set(Fp2 *rop,Fp2 *op){
  Fp_set(&(rop->x0),&(op->x0));
  Fp_set(&(rop->x1),&(op->x1));
}
/*EFp2 function*/
void EFp2_random(EFp2 *rop){
    /*Fp2_tmp6=x^3+ax+b*/
    while(1){
      Fp2_random(&rop->x);
      Fp2_mul(&Fp2_tmp6,&rop->x,&rop->x);
      Fp2_mul(&Fp2_tmp6,&Fp2_tmp6,&rop->x);
      Fp2_mul(&Fp2_tmp5,&rop->x,&Fp2_a);
      Fp2_add(&Fp2_tmp6,&Fp2_tmp6,&Fp2_tmp5);
      Fp2_add(&Fp2_tmp6,&Fp2_tmp6,&Fp2_b);
      if(Fp2_legendre(&Fp2_tmp6)==1) break;
    }
    Fp2_sqrt(&rop->y,&Fp2_tmp6);
    rop->infinity=0;
}
void EFp2_DBL(EFp2 *rop, EFp2 *op){
  if(op->infinity==1){
    rop->infinity = 1;
  }else if(Fp2_cmp(&(op->y),&Fp2_zero)==0){
    rop->infinity = 1;
  }else{
    EFp2_set(&EFp2_tmp1,op);

    Fp2_add(&Fp2_tmp7,&EFp2_tmp1.y,&EFp2_tmp1.y);
    Fp2_inv(&Fp2_tmp7,&Fp2_tmp7);

    Fp2_mul(&Fp2_tmp5,&EFp2_tmp1.x,&EFp2_tmp1.x);
    Fp2_add(&Fp2_tmp6,&Fp2_tmp5,&Fp2_tmp5);
    Fp2_add(&Fp2_tmp5,&Fp2_tmp5,&Fp2_tmp6);
    Fp2_add(&Fp2_tmp5,&Fp2_tmp5,&Fp2_a);

    Fp2_mul(&Fp2_tmp6,&Fp2_tmp7,&Fp2_tmp5);
    Fp2_mul(&Fp2_tmp7,&Fp2_tmp6,&Fp2_tmp6);

    Fp2_add(&Fp2_tmp5,&EFp2_tmp1.x,&EFp2_tmp1.x);
    Fp2_sub(&rop->x,&Fp2_tmp7,&Fp2_tmp5);
    Fp2_sub(&Fp2_tmp7,&EFp2_tmp1.x,&rop->x);
    Fp2_mul(&Fp2_tmp5,&Fp2_tmp6,&Fp2_tmp7);
    Fp2_sub(&rop->y,&Fp2_tmp5,&EFp2_tmp1.y);
  }
}

void EFp2_ADD(EFp2 *rop, EFp2 *op1, EFp2 *op2){
  if(op1->infinity==1){
      EFp2_set(rop,op2);
      return;
  }else if(op2->infinity==1){
      EFp2_set(rop,op1);
      return;
  }else if(Fp2_cmp(&op1->x,&op2->x)==0){
      if(Fp2_cmp(&op1->y,&op2->y)!=0){
          rop->infinity=1;
          return;
      }else{
          EFp2_DBL(rop,op1);
          return;
      }
  }
  EFp2_set(&EFp2_tmp1,op1);
  EFp2_set(&EFp2_tmp2,op2);

  Fp2_sub(&Fp2_tmp7,&EFp2_tmp2.x,&EFp2_tmp1.x);
  Fp2_inv(&Fp2_tmp7,&Fp2_tmp7);
  Fp2_sub(&Fp2_tmp5,&EFp2_tmp2.y,&EFp2_tmp1.y);
  Fp2_mul(&Fp2_tmp6,&Fp2_tmp7,&Fp2_tmp5);
  Fp2_mul(&Fp2_tmp7,&Fp2_tmp6,&Fp2_tmp6);


  Fp2_sub(&Fp2_tmp5,&Fp2_tmp7,&EFp2_tmp1.x);
  Fp2_sub(&rop->x,&Fp2_tmp5,&EFp2_tmp2.x);

  Fp2_sub(&Fp2_tmp7,&EFp2_tmp1.x,&rop->x);
  Fp2_mul(&Fp2_tmp5,&Fp2_tmp6,&Fp2_tmp7);
  Fp2_sub(&rop->y,&Fp2_tmp5,&EFp2_tmp1.y);
}
void EFp2_SCM(EFp2 *rop, EFp2 *op, mpz_t scalar){
  if(EFp_tmp1.infinity==1){
    rop->infinity=1;
    return;
  }
  int i;
  EFp2_set(&EFp2_tmp3,op);
  EFp2_tmp4.infinity=1;
  for(i=(int)(mpz_sizeinbase(scalar,2)-1);i>=0;i--){
    EFp2_DBL(&EFp2_tmp4,&EFp2_tmp4);
    if(mpz_tstbit(scalar, i)==1){
      EFp2_ADD(&EFp2_tmp4,&EFp2_tmp4,&EFp2_tmp3);
    }
  }
  EFp2_set(rop,&EFp2_tmp4);
}
int EFp2_cmp(EFp2 *op1, EFp2 *op2){
  if(Fp2_cmp(&op1->x,&op2->x)==0 && Fp2_cmp(&op1->y,&op2->y)==0){
    return 0;
  }else if(op1->infinity==1&&op2->infinity==1){
    return 0;
  }else{
  return 1;
  }
}
void EFp2_init(EFp2 *rop){
  Fp2_init(&(rop->x));
  Fp2_init(&(rop->y));
  rop->infinity = 0;
}
void EFp2_clear(EFp2 *rop){
  Fp2_clear(&(rop->x));
  Fp2_clear(&(rop->y));
}
void EFp2_print(EFp2 *rop){
  if(rop->infinity==0){
    printf("(");
    mpz_out_str(stdout,10,rop->x.x0.x0);
    printf(",");
    mpz_out_str(stdout,10,rop->x.x1.x0);
    printf(",");
    mpz_out_str(stdout,10,rop->y.x0.x0);
    printf(",");
    mpz_out_str(stdout,10,rop->y.x1.x0);
    printf(")\n");
  }else{
    printf("0\n");
  }
}
void EFp2_set(EFp2 *rop, EFp2 *op){
  Fp2_set(&(rop->x), &(op->x));
  Fp2_set(&(rop->y), &(op->y));
  rop->infinity = op->infinity;
}
void EFp2_inv(EFp2 *rop,EFp2 *op){
  if(op->infinity==1){
    rop->infinity=1;
    return;
  }else{
    Fp2_set(&rop->x,&op->x);
    Fp2_neg(&rop->y,&op->y);
    rop->infinity=0;
  }
}
void EFp2_sub(EFp2 *rop,EFp2 *op1,EFp2 *op2){
  EFp2_inv(&EFp2_tmp3,op2);
  EFp2_ADD(rop,op1,&EFp2_tmp3);
}

/*=============================================*/

void EFp2_lineTP(Fp2 *rop,EFp2 *q,EFp2 *p,EFp2 *t){/*L_t,p(q)*/
  if(Fp2_cmp(&(p->x),&(t->x))==0){/*Xp=Xtの時の例外処理*/
    Fp2_sub(rop,&(q->x),&(t->x));
  }else if(q->infinity==1){
    printf("\nTP::q.infinity==1\n");
    getchar();
  }else if(p->infinity==1){
    printf("\nTP::p.infinity==1\n");
    getchar();
  }else if(t->infinity==1){
    printf("\nTP::t.infinity==1\n");
    getchar();
  }else{
    EFp2 temp;
    EFp2_init(&temp);
    Fp2 slope;
    Fp2_init(&slope);

    //傾き計算
    Fp2_sub(&(temp.y), &(p->y), &(t->y));
    Fp2_sub(&(temp.x), &(p->x), &(t->x));
    Fp2_div(&slope, &(temp.y), &(temp.x));
    //temp.x=slope*(xq-xp)
    Fp2_sub(&(temp.x), &(q->x), &(p->x));
    Fp2_mul(&(temp.x), &slope, &(temp.x));
    //temp.y=(yq-yp)
    Fp2_sub(&(temp.y), &(q->y), &(p->y));
    Fp2_sub(rop, &(temp.y), &(temp.x));
    Fp2_clear(&slope);
    EFp2_clear(&temp);
  }
}
void EFp2_lineTT(Fp2 *rop,EFp2 *q,EFp2 *t){/*l_T,T(Q)*/
  if(Fp2_cmp(&(t->y),&Fp2_zero)==0){
    /*
    printf("\n=====check line ttf cmp====\n");
    getchar();
    */
    Fp2_sub(rop,&(q->x),&(t->x));
  }else if(q->infinity==1){
    printf("\nTT::q.infinity==1\n");
    getchar();
  }else if(t->infinity==1){
    printf("\nTT::t.infinity==1\n");
    getchar();
  }else{
    EFp2 temp;
    EFp2_init(&temp);
    Fp2 slope;
    Fp2_init(&slope);

    //傾き計算
    //temp.x=3*xt^2
    Fp2_pow_ui(&(temp.x), &(t->x), 2);
    Fp2_mul_ui(&(temp.x), &(temp.x), 3);

    //Fp2_add(&temp.x,&temp.x,&Fp2_a);
    //tenp.y=2*yt
    Fp2_mul_ui(&(temp.y), &(t->y), 2);
    //slope=temp.x/temp.y
    Fp2_div(&slope, &(temp.x), &(temp.y));

    //temp.x=slope*(xq-xt)
    Fp2_sub(&(temp.x), &(q->x), &(t->x));
    Fp2_mul(&(temp.x), &slope, &(temp.x));

    //temp.y=yq-yt
    Fp2_sub(&(temp.y), &(q->y), &(t->y));
    Fp2_sub(rop,&(temp.y),&(temp.x));

    EFp2_clear(&temp);
    Fp2_clear(&slope);
  }
}
void Fp2_Miller(Fp2 *rop,mpz_t s,EFp2 *p,EFp2 *q){
  Fp2 f,Fp2_temp2;
  EFp2 t;
  Fp2_init(&f);
  Fp2_init(&Fp2_temp2);
  EFp2_init(&t);
  int i;
  //f<-1,T<-P
  Fp2_set(&f,&Fp2_unity);
  EFp2_set(&t,p);

  for(i=(int)(mpz_sizeinbase(s,2))-2;i>=0;i--){
    //Fp2_pow_ui(&f,&f,2);
    Fp2_mul(&f,&f,&f);
    //void EFp2_lineTT(Fp2 *rop,EFp2 *q,EFp2 *t)
    EFp2_lineTT(&Fp2_temp2,q,&t);
    Fp2_mul(&f,&f,&Fp2_temp2);
    EFp2_DBL(&t,&t);
    if(mpz_tstbit(s,i)==1){
      //void EFp2_lineTP(Fp2 *rop,EFp2 *q,EFp2 *p,EFp2 *t)
      EFp2_lineTP(&Fp2_temp2,q,p,&t);
      Fp2_mul(&f,&f,&Fp2_temp2);
      EFp2_ADD(&t,&t,p);
    }
    if(i==0) break;
  }
  Fp2_set(rop,&f);

  /*-------------print---------------*/
  printf("Miller:: ANS=");
  Fp2_print(rop);

  Fp2_clear(&f);
  Fp2_clear(&Fp2_temp2);
  EFp2_clear(&t);
}
/*===================================================*/
void Fp2_Miller_v(Fp2 *rop,mpz_t s,EFp2 *p,EFp2 *q){
  Fp2 f,Fp2_temp2,Vtt,Vtp;
  EFp2 t;
  Fp2_init(&f);
  Fp2_init(&Fp2_temp2);
  EFp2_init(&t);
  Fp2_init(&Vtt);
  Fp2_init(&Vtp);
  mp_bitcnt_t i;
  //f<-1,T<-P
  Fp2_set(&f,&Fp2_unity);
  EFp2_set(&t,p);

  for(i=(int)(mpz_sizeinbase(s,2))-2;i>=0;i--){
    Fp2_pow_ui(&f,&f,2);
    //void EFp2_lineTT(Fp2 *rop,EFp2 *q,EFp2 *t)
    EFp2_lineTT(&Fp2_temp2,q,&t);
    Fp2_mul(&f,&f,&Fp2_temp2);
    EFp2_DBL(&t,&t);

    /*=======add========*/

    if(t.infinity==1){
      Fp_set_ui(&Vtt.x0,1);
      getchar();
      printf("infinity\n");
    }else{
      Fp2_sub(&Vtt,&q->x,&t.x);
    }
    Fp2_div(&f,&f,&Vtt);

    /*==================*/

    if(mpz_tstbit(s,i)==1){
      //void EFp2_lineTP(Fp2 *rop,EFp2 *q,EFp2 *p,EFp2 *t)
      EFp2_lineTP(&Fp2_temp2,q,p,&t);
      Fp2_mul(&f,&f,&Fp2_temp2);
      EFp2_ADD(&t,&t,p);

      /*=======add========*/
      if(t.infinity==1){
        Fp_set_ui(&Vtp.x0,1);
        getchar();
        printf("infinity2\n");
      }else{
        Fp2_sub(&Vtp,&q->x,&t.x);
      }
      Fp2_div(&f,&f,&Vtp);
      /*==================*/
    }
    if(i==0) break;
  }
  Fp2_set(rop,&f);

  /*-------------print---------------*/
  printf("Miller:: ANS=");
  Fp2_print(rop);

  Fp2_clear(&f);
  Fp2_clear(&Fp2_temp2);
  EFp2_clear(&t);
}
/*===================================================*/
void Fp2paring(Fp2 *rop,EFp2 *p,EFp2 *q){
  mpz_t temp;
  mpz_init(temp);
  //最終べき
  mpz_sub_ui(temp,Fp_p,1);
  Fp2_Miller_v(rop,SS_r,p,q);
  Fp2_pow(rop,rop,temp);
}
int check_generator(EFp2 *op){
  int flag;
  mpz_t s1,s2,s3;
  EFp2 tmp,check;
  mpz_init(s1);
  mpz_init(s2);
  mpz_init(s3);
  mpz_set_ui(s1,4);
  mpz_set_ui(s2,5);
  mpz_set_ui(s3,7);
  EFp2_init(&tmp);
  EFp2_init(&check);
  EFp2_set(&tmp,op);

  flag=0;
  /*  生成元かどうかのチェック*/
  EFp2_SCM(&check,&tmp,s1);
  if(check.infinity!=1){
    EFp2_SCM(&check,&tmp,s2);
    if(check.infinity!=1){
      EFp2_SCM(&check,&tmp,s3);
      if(check.infinity!=1){
        flag=1;
      }
    }
  }
  return flag;
}
void EFp2_set_gP(EFp2 *P){
  EFp tmp;
  EFp_init(&tmp);
  while(1){
    EFp_random(&tmp);
    Fp_set(&P->x.x0,&tmp.x);
    Fp_set(&P->y.x0,&tmp.y);
    if(check_generator(P)==1) break;
  }
}
void EFp2_set_gQ(EFp2 *Q){
  EFp2 tmp,tmp2;
  EFp2_init(&tmp);
  EFp2_init(&tmp2);
  while(1){
    EFp2_random(&tmp);
    Fp2_pow(&tmp2.x,&tmp.x,Fp_p);/*frobenius*/
    Fp2_pow(&tmp2.y,&tmp.y,Fp_p);/*frobenius*/
    EFp2_sub(Q,&tmp2,&tmp);
    if(check_generator(Q)==1) break;
  }
}
void EFp2_set_P(EFp2 *P){
  EFp tmp;
  EFp_init(&tmp);
  EFp_random(&tmp);
  Fp_set(&P->x.x0,&tmp.x);
  Fp_set(&P->y.x0,&tmp.y);
}

void EFp2_set_Q(EFp2 *Q){
  EFp2 tmp,tmp2;
  EFp2_init(&tmp);
  EFp2_init(&tmp2);
  EFp2_random(&tmp);
  Fp2_pow(&tmp2.x,&tmp.x,Fp_p);/*frobenius*/
  Fp2_pow(&tmp2.y,&tmp.y,Fp_p);/*frobenius*/
  EFp2_sub(Q,&tmp2,&tmp);
}
void inversiontest(){
  EFp2 P,Q,s;
  Fp2 rop;
  mpz_t Ps;
  mpz_init(Ps);
  EFp2_init(&P);
  EFp2_init(&Q);
  Fp2_init(&rop);
  EFp2_init(&s);
  /*
  EFp2_set_P(&P);
  EFp2_set_Q(&Q);
  */
  /*P=(67,38i)*/
  Fp_set_ui(&P.x.x0,67);
  Fp_set_ui(&P.x.x1,38);
  EFp2_print(&P);
  /*Q=(59,-54)*/
  Fp_set_ui(&Q.x.x0,59);
  Fp_sub_ui(&Q.y.x0,&Fp_zero,54);
  EFp2_print(&Q);

/*
  EFp2_SCM(&s,&P,SS_r);
  printf("P=");
  EFp2_print(&s);

  EFp2_SCM(&s,&Q,SS_r);
  printf("Q=");
  EFp2_print(&s);
*/
  mpz_set_ui(Ps,4);
  EFp2_SCM(&s,&P,Ps);
  EFp2_print(&s);
  Fp2_Miller_v(&rop,Ps,&P,&Q);
  //Fp2paring(&rop,&P,&Q);
  //Fp2_print(&rop);


  /*総線形テスト*/
  /*
  EFp2 tmp1,tmp2;
  mpz_t s1,s2,s3;
  Fp2 tmp3;
  EFp2_init(&tmp1);
  EFp2_init(&tmp2);
  Fp2_init(&tmp3);
  mpz_init(s1);
  mpz_init(s2);
  mpz_init(s3);

  mpz_set_ui(s1,2);
  mpz_set_ui(s2,3);
  mpz_mul(s3,s1,s2);

  EFp2_SCM(&tmp1,&P,s1);
  EFp2_SCM(&tmp2,&Q,s2);
  Fp2paring(&tmp3,&tmp1,&tmp2);
  Fp2_pow(&rop,&rop,s3);

  if(Fp2_cmp(&rop,&tmp3)==0){
    printf("Inversiontest:: ok\n");
  }else{
    printf("Inversiontest::error!!\n");
  }
  */
}
void paringtest(){
  EFp2 p,q,ap,bq,r,EFp2_temp;
  Fp2 rop1,rop2;
  mpz_t temp,temp2;
  EFp p_temp;
  EFp_init(&p_temp);
  EFp2_init(&p);
  EFp2_init(&q);
  EFp2_init(&r);
  EFp2_init(&EFp2_temp);
  Fp2_init(&rop1);
  Fp2_init(&rop2);
  EFp2_init(&ap);
  EFp2_init(&bq);
  mpz_init(temp);
  mpz_init(temp2);

  EFp2_set_gP(&p);
  EFp2_set_gQ(&q);

  printf("\n------\n");
  printf("p=\n");
  EFp2_print(&p);
  printf("\n------\n");
  printf("q=\n");
  EFp2_print(&q);

/*
  printf("\ntest p\n");
  EFp2_SCM(&r,&p,SS_r);
  EFp2_print(&r);
  printf("\ntest q\n");
  EFp2_SCM(&r,&q,SS_r);
  EFp2_print(&r);
*/
  /*双線形性を確かめる*/
  //temp=2,temp2=3
  mpz_set_ui(temp,3);
  mpz_set_ui(temp2,13);
  //ap=2*p,bp=3*q
  EFp2_SCM(&ap,&p,temp);
  EFp2_SCM(&bq,&q,temp2);
  //rop1=e(p,q),rop2=e(2*p,3*q)
  Fp2paring(&rop1,&p,&q);
  Fp2paring(&rop2,&ap,&bq);
  //temp=2*3
  mpz_mul(temp,temp2,temp);
  //rop1=e(p,q)^6
  Fp2_print(&rop1);
  Fp2_pow(&rop1,&rop1,temp);
  Fp2_print(&rop1);

  if(Fp2_cmp(&rop1,&rop2)==0){
    printf("CONGRATURATION!!!!!!!!\ne([a]p,[b]q)=e(p,q)^([a*b])\nLET'S A PARTY!!\n");
  }else{
    printf("error!!\n");
  }
}
void billinear(EFp2 *op1,EFp2 *op2){
  Fp2 tmp1,tmp2;
  EFp2 ap,bq;
  mpz_t s1,s2,s3;
  Fp2_init(&tmp1);
  Fp2_init(&tmp2);
  EFp2_init(&ap);
  EFp2_init(&bq);

  mpz_set_ui(s1,2);
  mpz_set_ui(s2,3);
  mpz_mul(s3,s1,s2);

  EFp2_SCM(&ap,op1,s1);
  EFp2_SCM(&bq,op2,s2);

  Fp2paring(&tmp1,op1,op2);
  Fp2paring(&tmp2,&ap,&bq);

  Fp2_pow(&tmp1,&tmp1,s3);
  if(Fp2_cmp(&tmp1,&tmp2)==0){
    printf("CONGRATURATION!!!!!!!!\ne([a]p,[b]q)=e(p,q)^([a*b])\nLET'S A PARTY!!\n");
  }else{
    printf("error!!\n");
  }
}
