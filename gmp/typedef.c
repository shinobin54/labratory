//https://na-inet.jp/na/gmp_ja/
//x = 1180591620717411305536
//p = 902561749388286096788557522098812120395206075319055706574520464280352613660532282580684031180264266084299760644204235424181611
//r = 1942668892225729084820684407073067588070130796063029175009226959140704863352059064321
//t = 1180591620717411305537
#include<stdio.h>
#include<gmp.h>
#include<time.h>
#include<unistd.h>
#include<stdlib.h>

#define p_value "902561749388286096788557522098812120395206075319055706574520464280352613660532282580684031180264266084299760644204235424181611"
#define r_value "1942668892225729084820684407073067588070130796063029175009226959140704863352059064321"
#define t_value "1180591620717411305537"
#define b_value 2
typedef struct{
  mpz_t x0;
}Fp;
typedef struct{
  Fp x;
  Fp y;
  int infinity;
}EFp;
mpz_t BLS_p;
mpz_t BLS_r;
mpz_t BLS_t;
gmp_randstate_t state;
Fp BLS_b;
void Fp_testinit();
void Fp_init(Fp *rop);
void Fp_set(Fp *rop, Fp *op);
void Fp_set_ui(Fp *rop, int op);
void Fp_print(Fp *rop);
void Fp_clear(Fp *rop);
void BLS_pr(mpz_t p, mpz_t r, mpz_t x);
void BLS_prt();
void Fp_random(Fp *rop);
void Fp_add(Fp *rop, const Fp *op1, const Fp *op2);
void Fp_add_ui(Fp *rop, const Fp *op1, int op2);
void Fp_mul(Fp *rop, const Fp *op1, const Fp *op2);
void Fp_inv(Fp *rop, const Fp *op1);
void Fp_mul_ui (Fp *rop, const Fp *op1, int ui);
void Fp_sub(Fp *rop, const Fp *op1, const Fp *op2);
void Fp_sub_ui (Fp *rop, const Fp *op1, int ui);
int Fp_cmp (const Fp *op1, const Fp *op2);
int Fp_cmp_ui (const Fp *op1, int op2);
void Fp_pow_gmp(Fp *rop, const Fp *base, const Fp *exp);
void Fp_pow_ui(Fp *rop, const Fp *base, int exp);
void Fp_pow(Fp *rop, Fp *op1, Fp *op2);
void Fp_division(Fp *rop, const Fp *op1, const Fp *op2);
void Fp_division_ui(Fp *rop, const Fp *op1, int op);
void mpz_print(mpz_t rop);
void BLS_search();
void fermat();
/*kadai3(i)*/
int Fp_legendre(Fp *a);
void Fp_sqrt(Fp *rop,Fp *a);
/*kadai4(i)*/
void EFp_random(EFp *rop);
void EFp_DBL(EFp *rop, EFp *op);
void EFp_ADD(EFp *rop, EFp *op1, EFp *op2);
void EFp_print(EFp *rop);
void EFp_init(EFp *rop);
void EFp_set(EFp *rop, EFp *op);
void EFp_clear(EFp *rop);
void EFp_SCM(EFp *rop, EFp *op, mpz_t s);
int main (){
  Fp_testinit();
  Fp a0,a1,a2;
  EFp b0,b1,b2;
  Fp_init(&a0);
  Fp_init(&a1);
  Fp_init(&a2);
  EFp_init(&b0);
  EFp_init(&b1);
  EFp_init(&b2);

  /*kadai 1 check(clear)*/
  /*
  BLS_prt();
  */

  /*kadai2 check(clear)*/

  /*
  Fp_random(&a0);
  sleep(1);
  Fp_random(&a1);
  Fp_add(&a2,&a1,&a0);
  Fp_print(&a1);
  printf("+\n");
  Fp_print(&a0);
  printf("=");
  Fp_print(&a2);
  */

  /*kadai3 check(clear)*/
  /*
  while(1){
    Fp_random(&a0);
    if(Fp_legendre(&a0)==1) break;
  }
  Fp_print(&a0);
  Fp_sqrt(&a0,&a0);
  Fp_print(&a0);
  */

  /*kadai4 check*/

  EFp_random(&b0);
  EFp_random(&b1);
  EFp_ADD(&b2,&b1,&b0);
  EFp_print(&b1);
  printf("+\n");
  EFp_print(&b0);
  printf("=");
  EFp_print(&b2);

  /*kadai5*/
  return 0;
}

void BLS_search(){
  mpz_t BLS_s,fin;
  EFp v,temp;
  mpz_init(fin);
  mpz_init(BLS_s);
  EFp_init(&v);
  EFp_init(&temp);
  mpz_set(BLS_s,BLS_p);
  mpz_add_ui(BLS_s,BLS_s,1);
  mpz_sub(BLS_s,BLS_s,BLS_t);
  mpz_sub_ui(fin,BLS_p,1);

  int flag;
  while(1){
    Fp_random(&BLS_b);
    flag = 0;
    Fp_set_ui(&(v.x),0);
    do{
      v.infinity=0;
      Fp_pow_ui(&(temp.y),&(v.x),3);
      Fp_add(&(temp.y),&(temp.y),&BLS_b);
      if(Fp_legendre(&(temp.y))==1){
        Fp_sqrt(&(v.y),&(temp.y));
        EFp_SCM(&v,&v,BLS_s);
        if(v.infinity==0){
          flag=0;
          break;
        }
      }
      if(mpz_cmp(temp.x.x0,fin)==0){
        flag=1;
        break;
      }
      Fp_add_ui(&(v.x),&(v.x),1);
    }while(mpz_cmp(v.x.x0,BLS_p)>=0);

    if(flag==1){  /*frag=1つまり全有理点をs倍したものが無限遠点になったということ*/
      break;
    }
  }
  printf("end");
  Fp_print(&BLS_b);
}

void Fp_testinit(){
  mpz_init(BLS_p);
  mpz_init(BLS_r);
  mpz_init(BLS_t);

  mpz_set_str(BLS_p, p_value, 10);
  mpz_set_str(BLS_r, r_value, 10);
  mpz_set_str(BLS_t, t_value, 10);

  Fp_init(&BLS_b);
  Fp_set_ui(&BLS_b, b_value);


}


/*kadai1 search x*/

void BLS_prt(){
  mpz_t x,p,r,t;
  mpz_init(x);
  mpz_init(p);
  mpz_init(r);
  mpz_init(t);

  mpz_set_ui(x, 2);
  mpz_pow_ui(x,x,70);

  while(1){
    BLS_pr(p,r,x);
    if( mpz_probab_prime_p(p,50)!=0 && mpz_probab_prime_p(r,50)!=0){
      break;
    }
    mpz_add_ui(x,x,1);
  }
  mpz_add_ui(t,x,1);

  printf("x = ");
  mpz_print(x);
  printf("p = ");
  mpz_print(p);
  printf("r = ");
  mpz_print(r);
  printf("t = ");
  mpz_print(t);
  printf("\n");

  mpz_clear(x);
  mpz_clear(p);
  mpz_clear(r);
  mpz_clear(t);

}

void BLS_pr(mpz_t p, mpz_t r, mpz_t x){
  mpz_t temp,temp2;

  mpz_init(temp);
  mpz_init(temp2);
  int i;

  mpz_pow_ui(temp, x, 4);
  mpz_pow_ui(temp2, x, 2);
  mpz_sub(r, temp, temp2);
  mpz_add_ui(r, r, 1);

  mpz_sub_ui(temp, x, 1);
  mpz_pow_ui(temp, temp, 2);
  mpz_mul(temp, r, temp);
  mpz_set_ui(temp2, 3);
  if(mpz_divisible_p(temp, temp2)!=0){
    mpz_cdiv_q(p, temp, temp2);
    mpz_add(p, p, x);
  }else{
    mpz_set_ui(p,1);
  }
  mpz_clear(temp);
  mpz_clear(temp2);
}

/*kadai2(i)*/
void Fp_random(Fp *rop){
  gmp_randinit_mt(state);
  gmp_randseed_ui(state, (unsigned int) time(NULL));
  mpz_urandomm(rop->x0, state, BLS_p);
}

/*kadai(ii)*/
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

void Fp_add(Fp *rop, const Fp *op1, const Fp *op2){ //BLS_pを法とした足し算の定義
  mpz_add(rop->x0, op1->x0, op2->x0);
  mpz_mod(rop->x0, rop->x0, BLS_p);
}

void Fp_add_ui(Fp *rop, const Fp *op1, int op2){ //moduleを法とした足し算の定義
  mpz_add_ui(rop->x0, op1->x0, op2);
  mpz_mod(rop->x0, rop->x0, BLS_p);
}

void Fp_sub(Fp *rop, const Fp *op1, const Fp *op2){
  mpz_sub(rop->x0, op1->x0, op2->x0);
  mpz_mod(rop->x0, rop->x0, BLS_p);
}

void Fp_sub_ui (Fp *rop, const Fp *op1, int ui){
  mpz_sub_ui (rop->x0, op1->x0, ui);
  mpz_mod(rop->x0, rop->x0, BLS_p);
}

void Fp_mul(Fp *rop, const Fp *op1, const Fp *op2){ //moduleを法とした掛け算の定義
  mpz_mul(rop->x0, op1->x0, op2->x0);
  mpz_mod(rop->x0, rop->x0, BLS_p);
}

void Fp_mul_ui (Fp *rop, const Fp *op1, int ui){
  mpz_mul_ui (rop->x0, op1->x0, ui);
  mpz_mod(rop->x0, rop->x0, BLS_p);
}

void Fp_inv(Fp *rop, const Fp *op1){
  mpz_invert(rop->x0, op1->x0, BLS_p);
}

/*kadai2(iii)*/
void Fp_pow(Fp *rop, Fp *op1, Fp *op2){
  mp_bitcnt_t i;
  Fp temp;
  Fp_init(&temp);
  Fp_set_ui(&temp,1);
  i=0;
  for(i=(mpz_sizeinbase(op2->x0,2)-1);;i--){
    Fp_mul(&temp,&temp,&temp);
    if(mpz_tstbit(op2->x0, i)==1){
      Fp_mul(&temp,&temp,op1);
    }
    if(i==0){
      break;
    }
  }
  Fp_set(rop,&temp);
  Fp_clear(&temp);
}


void EFp_SCM(EFp *rop, EFp *op, mpz_t s){
  mp_bitcnt_t i;
  EFp temp,sum;
  EFp_init(&temp);
  EFp_init(&sum);
  EFp_set(&temp,op);

  sum.infinity=1;
  for(i=(mpz_sizeinbase(s,2)-1);;i--){
    EFp_DBL(&sum,&sum);
    if(mpz_tstbit(s,i)==1){
      EFp_ADD(&sum,&sum,&temp);
    }
    if(i==0) break;
  }
  EFp_set(rop,&sum);
  EFp_clear(&sum);
  EFp_clear(&temp);
}

void fermat(){
  Fp a,p;
  Fp_init(&a);
  Fp_init(&p);

  Fp_random(&a);
  mpz_set(p.x0,BLS_p);
  Fp_sub_ui(&p,&p,1);
  Fp_pow(&a,&a,&p);
  printf("a^(p-1)=");
  Fp_print(&a);

  Fp_clear(&a);
  Fp_clear(&p);
}

/*kadai3(i)*/
int Fp_legendre(Fp *a){
  return mpz_legendre (a->x0, BLS_p);
}

void Fp_sqrt(Fp *rop, Fp *a){
  Fp n,q,x,b,y,t,temp,temp2;
  unsigned int e,r,i,m;
  Fp_init(&n);
  Fp_init(&q);
  Fp_init(&x);
  Fp_init(&b);
  Fp_init(&y);
  Fp_init(&t);

  Fp_init(&temp);
  Fp_init(&temp2);
  //1
  while(1){
    Fp_random(&n);
    if(Fp_legendre(&n)==-1){
      break;
    }
  }
  //2
  mpz_sub_ui(temp.x0,BLS_p,1);
  e = mpz_scan1(temp.x0,0);
  mpz_fdiv_q_2exp(q.x0,temp.x0,e);
  //3
  Fp_pow(&y,&n,&q);
  r=e;
  Fp_sub_ui(&temp,&q,1);
  Fp_division_ui(&temp,&temp,2);
  Fp_pow(&x,a,&temp);
  //4
  Fp_pow_ui(&b,&x,2);
  Fp_mul(&b,&b,a);
  Fp_mul(&x,a,&x);
  //5
  while(Fp_cmp_ui(&b,1)!=0){
    Fp_set_ui(&temp2,1);
    m=0;
    Fp_pow_ui(&temp,&b,2);
    do{
      Fp_mul(&temp2,&temp2,&temp);
      m++;
    }while(Fp_cmp_ui(&temp2,1)!=0);
    r=r-m-1;
    Fp_set_ui(&temp2,2);
    if(r>=0){
      Fp_pow_ui(&temp2,&temp2,r);
    }else{
      r=-r;
      Fp_inv(&temp2,&temp2);
      Fp_pow_ui(&temp2,&temp2,r);
    }
    Fp_pow(&t,&y,&temp2);
    //ok
    Fp_pow_ui(&y,&t,2);
    r=m;
    Fp_mul(&x,&x,&t);
    Fp_mul(&b,&b,&y);
  }
  Fp_set(rop,&x);
  Fp_clear(&temp);
  Fp_clear(&temp2);
  Fp_clear(&n);
  Fp_clear(&q);
  Fp_clear(&x);
  Fp_clear(&b);
  Fp_clear(&y);
  Fp_clear(&t);
}

//kadai4(i)
void EFp_random(EFp *rop){
  EFp temp;
  EFp_init(&temp);
  while(1){
    Fp_random(&(rop->x));
    Fp_pow_ui(&(temp.y),&(rop->x),3);
    Fp_add(&(temp.y),&(temp.y),&BLS_b);
    if(Fp_legendre(&(temp.y))==1){
      break;
    }
  }
  Fp_sqrt(&(rop->y),&(temp.y));

  rop->infinity = 0;

  EFp_print(rop);
}

//kadai4(ii)
void EFp_DBL(EFp *rop, EFp *op){
  if(op->infinity==1){
    rop->infinity = 1;
  }else if(Fp_cmp_ui(&(op->y),0)==0){
    rop->infinity = 1;
  }else{
    EFp temp;
    EFp_init(&temp);

    rop->infinity = 0;
    Fp slope;
    Fp_init(&slope);

    Fp slope2;
    Fp_init(&slope2);

    //傾き計算
    Fp_pow_ui(&(temp.x), &(op->x), 2);
    Fp_mul_ui(&(temp.x), &(temp.x), 3);

    Fp_mul_ui(&(temp.y), &(op->y), 2);
    Fp_division(&slope, &(temp.x), &(temp.y));

    //傾きの2乗計算
    Fp_pow_ui(&slope2, &slope, 2);

    //xrの計算
    Fp_mul_ui(&(temp.x), &(op->x), 2);
    Fp_sub(&(temp.x), &slope2, &(temp.x));

    //yrの計算
    Fp_sub(&(temp.y), &(op->x), &(temp.x));
    Fp_mul(&(temp.y), &slope, &(temp.y));
    Fp_sub(&(temp.y), &(temp.y), &(op->y));

    EFp_set(rop,&temp);
    EFp_clear(&temp);
  }
}

void EFp_ADD(EFp *rop, EFp *op1, EFp *op2){

  if((op1->infinity==1)&&(op2->infinity==1)){
    rop->infinity = 1;
  }else if(op1->infinity==1){
    EFp_set(rop, op2);
    rop->infinity = 0;
  }else if(op2->infinity==1){
    EFp_set(rop, op1);
    rop->infinity = 0;
  }else if(Fp_cmp(&(op1->x), &(op2->x))==0){
    if(Fp_cmp(&(op1->y), &(op2->y))==0){
      EFp_DBL(rop, op1);
    }else{
      rop->infinity = 1;
    }
  }else{
    rop->infinity = 0;

    EFp temp;
    EFp_init(&temp);

    Fp slope;
    Fp_init(&slope);

    Fp slope2;
    Fp_init(&slope2);

    //傾き計算
    Fp_sub(&(temp.y), &(op1->y), &(op2->y));
    Fp_sub(&(temp.x), &(op1->x), &(op2->x));
    Fp_division(&slope, &(temp.y), &(temp.x));
    //傾きの2乗計
    Fp_pow_ui(&slope2, &slope, 2);
    //xrの計算
    Fp_add(&(temp.x), &(op1->x), &(op2->x));
    Fp_sub(&(temp.x), &slope2, &(temp.x));
    //yrの計算
    Fp_sub(&(temp.y), &(op1->x), &(temp.x));
    Fp_mul(&(temp.y), &(temp.y), &slope);
    Fp_sub(&(temp.y), &(temp.y), &(op1->y));

    EFp_set(rop,&temp);
    Fp_clear(&slope);
    Fp_clear(&slope2);
    EFp_clear(&temp);
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
    Fp_print(&(rop->x));
    printf(",");
    Fp_print(&(rop->y));
    printf(")\n");
  }else{
    printf("infinity = %d\n",rop->infinity);
  }
}


void EFp_set(EFp *rop, EFp *op){
  Fp_set(&(rop->x), &(op->x));
  Fp_set(&(rop->y), &(op->y));
  rop->infinity = op->infinity;
}

void Fp_pow_gmp(Fp *rop, const Fp *base, const Fp *exp){
  mpz_powm(rop->x0, base->x0, exp->x0, BLS_p);
}

int Fp_cmp (const Fp *op1, const Fp *op2){
  return mpz_cmp (op1->x0, op2->x0);
}

int Fp_cmp_ui (const Fp *op1, int op2){
  return mpz_cmp_ui (op1->x0, op2);
}

void Fp_pow_ui(Fp *rop, const Fp *base, int exp){
  mpz_powm_ui(rop->x0, base->x0, exp, BLS_p);
}

void Fp_division(Fp *rop, const Fp *op1, const Fp *op2){
  Fp inv;
  Fp_init(&inv);
  Fp_inv(&inv, op2);
  Fp_mul(rop, op1, &inv);
  Fp_clear(&inv);
}

void Fp_division_ui(Fp *rop, const Fp *op1, int op){
  Fp inv,op2;
  Fp_init(&inv);
  Fp_init(&op2);
  Fp_set_ui(&op2,op);
  Fp_inv(&inv, &op2);
  Fp_mul(rop, op1, &inv);
  Fp_clear(&op2);
  Fp_clear(&inv);
}

void Fp_print(Fp *rop){
  mpz_out_str(stdout,10,rop->x0);
  printf("\n");
}

void mpz_print(mpz_t rop){
  mpz_out_str(stdout,10,rop);
  printf("\n");
}
