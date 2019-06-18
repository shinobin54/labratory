#include<stdio.h>
#include<gmp.h>
#include<time.h>
#include<unistd.h>
#include<math.h>


#define x_value "151115727451828646839744"
#define p_value "543225985036183036420758698714240607972453215193144672893818981352702661536580511866253158031231564392535113253181410731511"
#define r_value "13848450846756179314508504915203280026129928319360397923245259185956976824255752121"
#define t_value "151115727451828646839745"
#define s_value "543225985036183036420758698714240607972453215193144672893818981352702661536580511866253158031231564241419385801352763891767"
#define b_value 16
#define a_value 0
/*
#define p_value "902561749388286096788557522098812120395206075319055706574520464280352613660532282580684031180264266084299760644204235424181611"
#define r_value "1942668892225729084820684407073067588070130796063029175009226959140704863352059064321"
#define t_value "1180591620717411305537"
#define s_value "902561749388286096788557522098812120395206075319055706574520464280352613660532282580684031180264266084298580052583518012876075"
#define b_value 5
#define a_value 0
*/
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
typedef struct{
  Fp2 x0;
  Fp2 x1;
  Fp2 x2;
}Fp6;
typedef struct{
  Fp6 x;
  Fp6 y;
  int infinity;
}EFp6;
typedef struct{
  Fp6 x0;
  Fp6 x1;
}Fp12;
typedef struct{
  Fp12 x;
  Fp12 y;
  int infinity;
}EFp12;
gmp_randstate_t state;
/*Fp,EFp*/
mpz_t BLS_p,BLS_r,BLS_t,BLS_s,BLS_x;
Fp Fp_unity,Fp_zero;
Fp BLS_b,BLS_a;
/*Fp2,EFp2*/
Fp Fp2_i,Fp2_i2;
mpz_t Fp2_p,EFp2_r,EFp2_t;
Fp2 Fp2_unity,Fp2_zero,Fp2_b,Fp2_a;
/*Fp6,EFp6*/
Fp2 Fp6_i,Fp6_i2,Fp6_i3;
mpz_t Fp6_p,EFp6_r,EFp6_t;
Fp6 Fp6_unity,Fp6_zero,Fp6_b,Fp6_a;
/*Fp12,EFp12*/
Fp6 Fp12_i,Fp12_i2;
mpz_t Fp12_p,EFp12_r,EFp12_t;
Fp12 Fp12_unity,Fp12_zero,Fp12_b,Fp12_a;
/*frobenius*/
Fp c6,c3,c12,c4,c5;

/*time*/
clock_t start,end;
void clockstart(){
  start=clock();
}
void paring_opt(Fp12 *rop,EFp12 *p,EFp12 *q);
void Miller_opt(Fp12 *rop,mpz_t x,EFp12 *p,EFp12 *q);
void clockend(){
  end=clock();
  printf("%.2f秒かかりました\n",(double)(end-start)/CLOCKS_PER_SEC);
}
void paring(Fp12 *rop,EFp12 *p,EFp12 *q);
void mpz_print(mpz_t rop);
/*Fp*/
void Fp_testinit();
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
void BLS_make();
/*EFp*/
void EFp_random(EFp *rop);
void EFp_DBL(EFp *rop, EFp *op);
void EFp_ADD(EFp *rop, EFp *op1, EFp *op2);
void EFp_SCM(EFp *rop, EFp *op, mpz_t s);
void EFp_print(EFp *rop);
void EFp_init(EFp *rop);
void EFp_set(EFp *rop, EFp *op);
void EFp_clear(EFp *rop);
/*Fp2*/
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

/*Fp2*/
void Fp6_make();
void Fp6_init(Fp6 *op);
int Fp6_cmp(Fp6 *op1,Fp6 *op2);
void Fp6_random(Fp6 *rop);
void Fp6_clear(Fp6 *op);
void Fp6_sub(Fp6 *rop, Fp6 *op1, Fp6 *op2);
void Fp6_add(Fp6 *rop, Fp6 *op1, Fp6 *op2);
void Fp6_mul(Fp6 *rop, Fp6 *op1, Fp6 *op2);
void Fp6_conju(Fp6 *rop1,Fp6 *rop2,Fp6 *op);
void Fp6_inv(Fp6 *rop, Fp6 *op);
void Fp6_print(Fp6 *rop);
void Fp6_set(Fp6 *rop,Fp6 *op);
void Fp6_pow(Fp6 *rop, Fp6 *op1, mpz_t scalar);
int Fp6_legendre(Fp6 *op);
void Fp6_mul_ui(Fp6 *rop,Fp6 *op1,unsigned int op2);
void Fp6_div(Fp6 *rop,Fp6 *op1,Fp6 *op2);
void Fp6_sqrt(Fp6 *rop,Fp6 *op);
void Fp6_pow_ui(Fp6 *rop,Fp6 *op1,unsigned int op2);

/*EFp6*/
void EFp6_random(EFp6 *rop);
void EFp6_DBL(EFp6 *rop, EFp6 *op);
void EFp6_ADD(EFp6 *rop, EFp6 *op1, EFp6 *op2);
void EFp6_SCM(EFp6 *rop, EFp6 *op, mpz_t s);
void EFp6_print(EFp6 *rop);
void EFp6_init(EFp6 *rop);
void EFp6_set(EFp6 *rop, EFp6 *op);
void EFp6_clear(EFp6 *rop);
int EFp6_cmp(EFp6 *op1, EFp6 *op2);
/*Fp12*/
void Fp12_Frobenius(Fp12 *rop,Fp12 *op);
void Fp12_Frobenius_init();
void Fp12_make();
void Fp12_init(Fp12 *op);
int Fp12_cmp(Fp12 *op1,Fp12 *op2);
void Fp12_random(Fp12 *rop);
void Fp12_clear(Fp12 *op);
void Fp12_sub(Fp12 *rop, Fp12 *op1, Fp12 *op2);
void Fp12_add(Fp12 *rop, Fp12 *op1, Fp12 *op2);
void Fp12_mul(Fp12 *rop, Fp12 *op1, Fp12 *op2);
void Fp12_conju(Fp12 *rop,Fp12 *op);
void Fp12_inv(Fp12 *rop, Fp12 *op);
void Fp12_print(Fp12 *rop);
void Fp12_set(Fp12 *rop,Fp12 *op);

void Fp12_pow(Fp12 *rop, Fp12 *op1, mpz_t scalar);
void Fp12_pow_ui(Fp12 *rop,Fp12 *op1,unsigned int op2);
void Fp12_div(Fp12 *rop,Fp12 *op1,Fp12 *op2);
void Fp12_sqrt(Fp12 *rop,Fp12 *op);
void Fp12_mul_ui(Fp12 *rop,Fp12 *op1,unsigned int op2);
int Fp12_legendre(Fp12 *op);
int Fp12_legendre3(Fp12 *op);
/*EFp12*/
void EFp12_sub(EFp12 *rop,EFp12 *op1,EFp12 *op2);
void EFp12_random(EFp12 *rop);
void EFp12_DBL(EFp12 *rop, EFp12 *op);
void EFp12_ADD(EFp12 *rop, EFp12 *op1, EFp12 *op2);
void EFp12_SCM(EFp12 *rop, EFp12 *op, mpz_t s);
void EFp12_print(EFp12 *rop);
void EFp12_init(EFp12 *rop);
void EFp12_set(EFp12 *rop, EFp12 *op);
void EFp12_clear(EFp12 *rop);
int EFp12_cmp(EFp12 *op1, EFp12 *op2);
void EFp12_inv(EFp12 *rop,EFp12 *op);
/*BLS search*/
void BLS_pr(mpz_t p, mpz_t r);
void BLS_make_random();
void BLS_search();
void p_make(EFp *rop);
void fermat();

void Miller(Fp12 *rop,mpz_t s,EFp12 *p,EFp12 *q);
void EFp12_lineTP(Fp12 *rop,EFp12 *q,EFp12 *p,EFp12 *t);
void EFp12_lineTT(Fp12 *rop,EFp12 *q,EFp12 *t);

/*以下編集中*/
void test2(){
  EFp12 p,q,ap,bq,r,EFp12_temp;
  Fp12 rop1,rop2;
  mpz_t temp,temp2;
  EFp p_temp;
  EFp_init(&p_temp);
  EFp12_init(&p);
  EFp12_init(&q);
  EFp12_init(&r);
  EFp12_init(&EFp12_temp);
  Fp12_init(&rop1);
  Fp12_init(&rop2);
  EFp12_init(&ap);
  EFp12_init(&bq);
  mpz_init(temp);
  mpz_init(temp2);


  p_make(&p_temp);
  Fp_set(&(p.x.x0.x0.x0),&(p_temp.x));
  Fp_set(&(p.y.x0.x0.x0),&(p_temp.y));

  EFp12_random(&r);
  mpz_pow_ui(temp,BLS_r,2);
  mpz_divexact(temp,EFp12_r,temp);
  EFp12_SCM(&r,&r,temp);
  Fp12_pow(&(EFp12_temp.x),&(r.x),BLS_p);
  Fp12_pow(&(EFp12_temp.y),&(r.y),BLS_p);
  EFp12_sub(&q,&EFp12_temp,&r);
  printf("\n------\n");
  printf("p=\n");
  EFp12_print(&p);
  printf("\n------\n");
  printf("q=\n");
  EFp12_print(&q);

  /*双線形性を確かめる*/
  //temp=2,temp2=3
  mpz_set_ui(temp,2);
  mpz_set_ui(temp2,3);
  //ap=2*p,bp=3*q
  EFp12_SCM(&ap,&p,temp);
  EFp12_SCM(&bq,&q,temp2);
  //rop1=e(p,q),rop2=e(2*p,3*q)
  paring_opt(&rop1,&p,&q);
  paring_opt(&rop2,&ap,&bq);
  //temp=2*3
  mpz_mul(temp,temp2,temp);
  //rop1=e(p,q)^6
  Fp12_pow(&rop1,&rop1,temp);

  if(Fp12_cmp(&rop1,&rop2)==0){
    printf("CONGRATURATION!!!!!!!!\ne([a]p,[b]q)=e(p,q)^([a*b])\nLET'S A PARTY!!\n");
  }else{
    printf("error!!\n");
  }
}
void EFp12_lineTP(Fp12 *rop,EFp12 *q,EFp12 *p,EFp12 *t){/*L_t,p(q)*/
  if(Fp12_cmp(&(p->x),&(t->x))==0){/*Xp=Xtの時の例外処理*/
    Fp12_sub(rop,&(q->x),&(t->x));
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
    EFp12 temp;
    EFp12_init(&temp);
    Fp12 slope;
    Fp12_init(&slope);

    //傾き計算
    Fp12_sub(&(temp.y), &(p->y), &(t->y));
    Fp12_sub(&(temp.x), &(p->x), &(t->x));
    Fp12_div(&slope, &(temp.y), &(temp.x));
    //temp.x=slope*(xq-xp)
    Fp12_sub(&(temp.x), &(q->x), &(p->x));
    Fp12_mul(&(temp.x), &slope, &(temp.x));
    //temp.y=(yq-yp)
    Fp12_sub(&(temp.y), &(q->y), &(p->y));
    Fp12_sub(rop, &(temp.y), &(temp.x));
    Fp12_clear(&slope);
    EFp12_clear(&temp);
  }
}
void EFp12_lineTT(Fp12 *rop,EFp12 *q,EFp12 *t){/*l_T,T(Q)*/
  if(Fp12_cmp(&(t->y),&Fp12_zero)==0){
    /*
    printf("\n=====check line ttf cmp====\n");
    getchar();
    */
    Fp12_sub(rop,&(q->x),&(t->x));
  }else if(q->infinity==1){
    printf("\nTT::q.infinity==1\n");
    getchar();
  }else if(t->infinity==1){
    printf("\nTT::t.infinity==1\n");
    getchar();
  }else{
    EFp12 temp;
    EFp12_init(&temp);
    Fp12 slope;
    Fp12_init(&slope);

    //傾き計算
    //temp.x=3*xt^2
    Fp12_pow_ui(&(temp.x), &(t->x), 2);
    Fp12_mul_ui(&(temp.x), &(temp.x), 3);
    //tenp.y=2*yt
    Fp12_mul_ui(&(temp.y), &(t->y), 2);
    //slope=temp.x/temp.y
    Fp12_div(&slope, &(temp.x), &(temp.y));

    //temp.x=slope*(xq-xt)
    Fp12_sub(&(temp.x), &(q->x), &(t->x));
    Fp12_mul(&(temp.x), &slope, &(temp.x));

    //temp.y=yq-yt
    Fp12_sub(&(temp.y), &(q->y), &(t->y));
    Fp12_sub(rop,&(temp.y),&(temp.x));

    EFp12_clear(&temp);
    Fp12_clear(&slope);
  }

}
void paring(Fp12 *rop,EFp12 *p,EFp12 *q){
  mpz_t temp;
  mpz_init(temp);
  //最終べき
  mpz_sub_ui(temp,Fp12_p,1);
  mpz_divexact(temp,temp,BLS_r);
  Miller(rop,BLS_r,p,q);
  Fp12_pow(rop,rop,temp);
}
void paring_opt(Fp12 *rop,EFp12 *p,EFp12 *q){
  mpz_t temp;
  mpz_init(temp);
  mpz_sub_ui(temp,Fp12_p,1);
  mpz_divexact(temp,temp,BLS_r);
  Miller_opt(rop,BLS_x,p,q);
  Fp12_pow(rop,rop,temp);
}
void Miller(Fp12 *rop,mpz_t s,EFp12 *p,EFp12 *q){
  Fp12 f,Fp12_temp2;
  EFp12 t;
  Fp12_init(&f);
  Fp12_init(&Fp12_temp2);
  EFp12_init(&t);
  mp_bitcnt_t i;
  //f<-1,T<-P
  Fp12_set(&f,&Fp12_unity);
  EFp12_set(&t,p);

  for(i=(int)(mpz_sizeinbase(s,2))-2;i>=0;i--){
    Fp12_pow_ui(&f,&f,2);
    //void EFp12_lineTT(Fp12 *rop,EFp12 *q,EFp12 *t)
    EFp12_lineTT(&Fp12_temp2,q,&t);
    Fp12_mul(&f,&f,&Fp12_temp2);
    EFp12_DBL(&t,&t);
    if(mpz_tstbit(s,i)==1){
      //void EFp12_lineTP(Fp12 *rop,EFp12 *q,EFp12 *p,EFp12 *t)
      EFp12_lineTP(&Fp12_temp2,q,p,&t);
      Fp12_mul(&f,&f,&Fp12_temp2);
      EFp12_ADD(&t,&t,p);
    }
    if(i==0) break;
  }
  Fp12_set(rop,&f);
  Fp12_clear(&f);
  Fp12_clear(&Fp12_temp2);
  EFp12_clear(&t);
}
void Miller_opt(Fp12 *rop,mpz_t x,EFp12 *p,EFp12 *q){
  Fp12 f,Fp12_temp2;
  EFp12 t;
  Fp12_init(&f);
  Fp12_init(&Fp12_temp2);
  EFp12_init(&t);
  mp_bitcnt_t i;
  //f<-1,T<-P
  Fp12_set(&f,&Fp12_unity);
  EFp12_set(&t,q);

  for(i=(int)(mpz_sizeinbase(x,2))-2;i>=0;i--){
    Fp12_pow_ui(&f,&f,2);
    //void EFp12_lineTT(Fp12 *rop,EFp12 *q,EFp12 *t)
    EFp12_lineTT(&Fp12_temp2,p,&t);
    Fp12_mul(&f,&f,&Fp12_temp2);
    EFp12_DBL(&t,&t);
    if(mpz_tstbit(x,i)==1){
      //void EFp12_lineTP(Fp12 *rop,EFp12 *q,EFp12 *p,EFp12 *t)
      EFp12_lineTP(&Fp12_temp2,p,q,&t);
      Fp12_mul(&f,&f,&Fp12_temp2);
      EFp12_ADD(&t,&t,q);
    }
    if(i==0) break;
  }
  Fp12_set(rop,&f);
  Fp12_clear(&f);
  Fp12_clear(&Fp12_temp2);
  EFp12_clear(&t);
}


void wail(){/*a+b=t a*b=p*/
  mpz_t temp1,temp2,p,q;
  int i;
  i=1;
  mpz_init(temp1);
  mpz_init(temp2);
  mpz_init(p);
  mpz_init(q);

  mpz_init(EFp2_r);
  mpz_init(EFp6_r);
  mpz_init(EFp12_r);

  mpz_init(EFp2_t);
  mpz_init(EFp6_t);
  mpz_init(EFp12_t);

  mpz_set(temp1,BLS_t);
  mpz_set_ui(temp2,2);
  while(1){
    if(i==2){
      mpz_set(EFp2_t,temp1);
    }

    if(i==6){
      mpz_set(EFp6_t,temp1);
    }
    if(i==12){
      mpz_set(EFp12_t,temp1);
      break;
    }
    mpz_mul(p,BLS_t,temp1);
    mpz_mul(q,BLS_p,temp2);
    mpz_set(temp2,temp1);
    mpz_sub(temp1,p,q);
    i++;
  }
  /*EFp2_r=EFp2_p + 1 - EFp2_t*/
  mpz_add_ui(temp1,Fp2_p,1);
  mpz_sub(EFp2_r,temp1,EFp2_t);
  mpz_add_ui(temp1,Fp6_p,1);
  mpz_sub(EFp6_r,temp1,EFp6_t);
  mpz_add_ui(temp1,Fp12_p,1);
  mpz_sub(EFp12_r,temp1,EFp12_t);

  mpz_clear(temp1);
  mpz_clear(temp2);
  mpz_clear(p);
  mpz_clear(q);
}

void Fp12_make(){
  Fp6_init(&Fp12_i);
/*
  while(1){
    Fp6_random(&Fp12_i);
    if(Fp6_legendre(&Fp12_i)==-1){
      break;
    }
  }
*/
  Fp_set_ui(&(Fp12_i.x1.x0),1);
  mpz_init(Fp12_p);
  mpz_pow_ui(Fp12_p,Fp6_p,2);
  Fp12_init(&Fp12_unity);
  Fp12_init(&Fp12_zero);
  Fp6_set(&(Fp12_unity.x0),&Fp6_unity);
  Fp6_set(&(Fp12_unity.x1),&Fp6_zero);
  Fp6_set(&(Fp12_zero.x0),&Fp6_zero);
  Fp6_set(&(Fp12_zero.x1),&Fp6_zero);
  Fp12_init(&Fp12_b);
  Fp12_init(&Fp12_a);
  Fp6_set(&(Fp12_b.x0),&Fp6_b);
  Fp6_set(&(Fp12_b.x1),&Fp6_zero);
  Fp6_set(&(Fp12_a.x0),&Fp6_a);
  Fp6_set(&(Fp12_a.x0),&Fp6_zero);
}
int Fp12_legendre(Fp12 *op){
  mpz_t scalar;
  Fp12 temp;
  Fp12_init(&temp);
  mpz_init(scalar);
  mpz_sub_ui(scalar,Fp12_p,1);
  mpz_divexact_ui(scalar,scalar,2);
  Fp12_pow(&temp,op,scalar);
  mpz_clear(scalar);
  if(Fp12_cmp(&temp,&Fp12_unity)==0||Fp12_cmp(&temp,&Fp12_zero)==0){
    Fp12_clear(&temp);
    return 1;
  }else{
    Fp12_clear(&temp);
    return -1;
  }
}
int Fp12_legendre3(Fp12 *op){
  mpz_t scalar;
  Fp12 temp;
  Fp12_init(&temp);
  mpz_init(scalar);
  mpz_sub_ui(scalar,Fp12_p,1);
  mpz_divexact_ui(scalar,scalar,3);
  Fp12_pow(&temp,op,scalar);
  mpz_clear(scalar);
  if(Fp12_cmp(&temp,&Fp12_unity)==0||Fp12_cmp(&temp,&Fp12_zero)==0){
    Fp12_clear(&temp);
    return 1;
  }else{
    Fp12_clear(&temp);
    return -1;
  }
}
void Fp12_pow_ui(Fp12 *rop,Fp12 *op1,unsigned int op2){
  mpz_t temp;
  mpz_init(temp);
  mpz_set_ui(temp,op2);
  Fp12_pow(rop,op1,temp);
  mpz_clear(temp);
}
void Fp12_mul_ui(Fp12 *rop,Fp12 *op1,unsigned int op2){/*スカラー倍算*/
  Fp6_mul_ui(&(rop->x0),&(op1->x0),op2);
  Fp6_mul_ui(&(rop->x1),&(op1->x1),op2);
}
void Fp12_sqrt(Fp12 *rop, Fp12 *a){
  mpz_t q,mpz_temp,mpz_temp2;
  Fp12 n,x,y,b,t,Fp12_temp,Fp12_temp2;
  unsigned int e,r,i,m;
  Fp12_init(&n);
  Fp12_init(&x);
  Fp12_init(&b);
  Fp12_init(&y);
  Fp12_init(&t);
  Fp12_init(&Fp12_temp);
  Fp12_init(&Fp12_temp2);

  mpz_init(q);
  mpz_init(mpz_temp);
  mpz_init(mpz_temp2);

  //1
  while(1){
    Fp12_random(&n);
    if(Fp12_legendre(&n)==-1){
      break;
    }
  }
  //2
  mpz_sub_ui(mpz_temp,Fp12_p,1);
  e = mpz_scan1(mpz_temp,0);
  mpz_fdiv_q_2exp(q,mpz_temp,e);
  //3
  Fp12_pow(&y,&n,q);
  r=e;
  mpz_sub_ui(mpz_temp,q,1);
  mpz_divexact_ui(mpz_temp,mpz_temp,2);
  Fp12_pow(&x,a,mpz_temp);
  //4
  Fp12_pow_ui(&b,&x,2);
  Fp12_mul(&b,&b,a);
  Fp12_mul(&x,a,&x);
  //5
  while(Fp12_cmp(&b,&Fp12_unity)!=0){
    m=0;
    Fp12_set(&Fp12_temp,&b);
    while(Fp12_cmp(&Fp12_temp,&Fp12_unity)!=0){
      Fp12_pow_ui(&Fp12_temp,&Fp12_temp,2);
      m++;
    }
    r=r-m-1;
    mpz_set_ui(mpz_temp,2);
    mpz_pow_ui(mpz_temp,mpz_temp,r);
    Fp12_pow(&t,&y,mpz_temp);
    //ok
    Fp12_pow_ui(&y,&t,2);
    r=m;
    Fp12_mul(&x,&x,&t);
    Fp12_mul(&b,&b,&y);
  }
  Fp12_set(rop,&x);
  mpz_clear(mpz_temp);
  Fp12_clear(&Fp12_temp2);
  Fp12_clear(&Fp12_temp);
  Fp12_clear(&n);
  mpz_clear(q);
  Fp12_clear(&x);
  Fp12_clear(&b);
  Fp12_clear(&y);
  Fp12_clear(&t);
}
void Fp12_pow(Fp12 *rop, Fp12 *op1, mpz_t scalar){
  mp_bitcnt_t i;
  Fp12 temp;
  Fp12_init(&temp);
  /*temp=unity*/
  Fp12_set(&temp,&Fp12_unity);

  i=0;
  for(i=(mpz_sizeinbase(scalar,2)-1);;i--){
    Fp12_mul(&temp,&temp,&temp);
    if(mpz_tstbit(scalar, i)==1){
      Fp12_mul(&temp,&temp,op1);
    }
    if(i==0){
      break;
    }
  }
  Fp12_set(rop,&temp);
  Fp12_clear(&temp);
}
void Fp12_init(Fp12 *op){
  Fp6_init(&(op->x0));
  Fp6_init(&(op->x1));
}
void Fp12_clear(Fp12 *op){
  Fp6_clear(&(op->x0));
  Fp6_clear(&(op->x1));
}
void Fp12_random(Fp12 *rop){
  Fp6_random(&(rop->x0));
  Fp6_random(&(rop->x1));
}
void Fp12_add(Fp12 *rop, Fp12 *op1, Fp12 *op2){
  Fp6_add(&(rop->x0),&(op1->x0),&(op2->x0));
  Fp6_add(&(rop->x1),&(op1->x1),&(op2->x1));
}
void Fp12_sub(Fp12 *rop, Fp12 *op1, Fp12 *op2){
  Fp6_sub(&(rop->x0),&(op1->x0),&(op2->x0));
  Fp6_sub(&(rop->x1),&(op1->x1),&(op2->x1));
}
void Fp12_mul(Fp12 *rop, Fp12 *op1, Fp12 *op2){
  Fp12 temp,temp2;
  Fp12_init(&temp);
  Fp12_init(&temp2);

  Fp6_mul(&(temp.x0),&(op1->x0),&(op2->x0));
  Fp6_mul(&(temp2.x0),&(op1->x1),&(op2->x1));
  Fp6_mul(&(temp2.x0),&(temp2.x0),&Fp12_i);

  Fp6_mul(&(temp.x1),&(op1->x0),&(op2->x1));
  Fp6_mul(&(temp2.x1),&(op1->x1),&(op2->x0));

  Fp12_add(rop,&temp,&temp2);
  Fp12_clear(&temp);
  Fp12_clear(&temp2);
}
void Fp12_conju(Fp12 *rop,Fp12 *op){
  Fp6_set(&(rop->x0),&(op->x0));
  Fp6_sub(&(rop->x1),&Fp6_zero,&(op->x1));
}
void Fp12_inv(Fp12 *rop, Fp12 *op){
  Fp12 temp1,temp2;
  Fp12_init(&temp1);
  Fp12_init(&temp2);
  Fp12_conju(&temp1,op);
  //temp2=(a+bi)*(a-bi)
  Fp12_mul(&temp2,op,&temp1);
  Fp6_div(&(rop->x0),&(temp1.x0),&(temp2.x0));
  Fp6_div(&(rop->x1),&(temp1.x1),&(temp2.x0));
  Fp12_clear(&temp1);
  Fp12_clear(&temp2);
}
void Fp12_div(Fp12 *rop,Fp12 *op1,Fp12 *op2){
  Fp12 inv;
  Fp12_init(&inv);
  Fp12_inv(&inv,op2);
  Fp12_mul(rop,op1,&inv);
  Fp12_clear(&inv);
}
int Fp12_cmp(Fp12 *op1,Fp12 *op2){
  if(Fp6_cmp(&(op1->x0),&(op2->x0))==0 && Fp6_cmp(&(op1->x1),&(op2->x1))==0){
    return 0;
  }else{
    return 1;
  }
}
void Fp12_print(Fp12 *rop){
  printf("(");
  mpz_out_str(stdout,10,rop->x0.x0.x0.x0);
  printf(",");
  mpz_out_str(stdout,10,rop->x0.x0.x1.x0);
  printf(",");
  mpz_out_str(stdout,10,rop->x0.x1.x0.x0);
  printf(",");
  mpz_out_str(stdout,10,rop->x0.x1.x1.x0);
  printf(",");
  mpz_out_str(stdout,10,rop->x0.x2.x0.x0);
  printf(",");
  mpz_out_str(stdout,10,rop->x0.x2.x1.x0);
  printf(",");

  mpz_out_str(stdout,10,rop->x1.x0.x0.x0);
  printf(",");
  mpz_out_str(stdout,10,rop->x1.x0.x1.x0);
  printf(",");
  mpz_out_str(stdout,10,rop->x1.x1.x0.x0);
  printf(",");
  mpz_out_str(stdout,10,rop->x1.x1.x1.x0);
  printf(",");
  mpz_out_str(stdout,10,rop->x1.x2.x0.x0);
  printf(",");
  mpz_out_str(stdout,10,rop->x1.x2.x1.x0);
  printf(")\n");
}
void Fp12_set(Fp12 *rop,Fp12 *op){
  Fp6_set(&(rop->x0),&(op->x0));
  Fp6_set(&(rop->x1),&(op->x1));
}

int main(){
  Fp_testinit();
  BLS_make();
  //Fp12_Frobenius_init();
  Fp2_make();
  Fp6_make();
  Fp12_make();
  wail();
  test2();
  mpz_t mpz_s;
  unsigned int i,s;
  Fp Fp_1,Fp_2;
  Fp2 Fp2_1,Fp2_2;
  Fp6 Fp6_1,Fp6_2,Fp6_3;
  Fp12 Fp12_1,Fp12_2,Fp12_3;

  EFp2 EFp2_1,EFp2_2,EFp2_3;
  EFp6 EFp6_1,EFp6_2,EFp6_3;
  EFp12 EFp12_1,EFp12_2,EFp12_3;
  mpz_init(mpz_s);
  Fp2_init(&Fp2_1);
  Fp2_init(&Fp2_2);

  EFp2_init(&EFp2_1);
  EFp2_init(&EFp2_2);
  EFp2_init(&EFp2_3);

  EFp6_init(&EFp6_1);
  EFp6_init(&EFp6_2);
  EFp6_init(&EFp6_3);
  Fp_init(&Fp_1);
  Fp_init(&Fp_2);
  Fp6_init(&Fp6_1);
  Fp6_init(&Fp6_2);
  Fp6_init(&Fp6_3);


  Fp12_init(&Fp12_1);
  Fp12_init(&Fp12_2);
  Fp12_init(&Fp12_3);


  EFp12_init(&EFp12_1);
  EFp12_init(&EFp12_2);
  EFp12_init(&EFp12_3);

  EFp6_random(&EFp6_1);
  /*
  Fp6_random(&Fp6_1);
  Fp6_random(&Fp6_2);
  Fp6_random(&Fp6_3);
  Fp12_random(&Fp12_1);
  Fp12_random(&Fp12_2);
  Fp12_random(&Fp12_3);
  EFp12_random(&EFp12_1);
  EFp12_random(&EFp12_2);
  EFp12_random(&EFp12_3);
*/

/*
  EFp12_init(&EFp12_1);
  EFp12_init(&EFp12_2);
  EFp12_init(&EFp12_3);
  EFp12_random(&EFp12_1);
  EFp12_random(&EFp12_2);
  EFp12_random(&EFp12_3);
*/
  /*
  Fp6_add(&Fp6_3,&Fp6_1,&Fp6_2);
  Fp6_sub(&Fp6_3,&Fp6_3,&Fp6_2);
  if(Fp6_cmp(&Fp6_3,&Fp6_1)==0){
    printf("ok!");
  }else{
    printf("unko!");
  }
  */
/*
  Fp6_3
  if(Fp6_cmp(&Fp6_3,&Fp6_1)==0){
    printf("ok!");
  }else{
    printf("unko!");
  }
  */
  /*test(conju,mul)*/
  /*
  Fp6_print(&Fp6_3);
  Fp6_print(&Fp6_1);
  printf("===\n");
  Fp6_mul(&Fp6_2,&Fp6_3,&Fp6_1);
  Fp6_div(&Fp6_2,&Fp6_2,&Fp6_3);
  Fp6_print(&Fp6_2);
  Fp6_print(&Fp6_1);
  if(Fp6_cmp(&Fp6_2,&Fp6_1)==0){
    printf("ok!");
  }else{
    printf("unko!");
  }
  */
  /*
  printf("Fp6_i=");
  Fp2_print(&Fp6_i);
  printf("========\n");
  Fp6_print(&Fp6_3);
  printf("========\n");
  Fp6_print(&Fp6_2);
  printf("========\n");
  */
  /*
  Fp6_mul(&Fp6_3,&Fp6_3,&Fp6_1);
  Fp6_mul(&Fp6_3,&Fp6_3,&Fp6_2);
  Fp6_print(&Fp6_3);
  */

  //check

/*
  s=90;
  mpz_t a;
  mpz_init(a);
  mpz_set_ui(a,s);
  Fp12_random(&Fp12_1);
  Fp12_set(&Fp12_2,&Fp12_unity);


  for(i=0;i<s;i++){
    Fp12_mul(&Fp12_2,&Fp12_2,&Fp12_1);
  }

  Fp12_pow_ui(&Fp12_1,&Fp12_1,s);
  if(Fp12_cmp(&Fp12_1,&Fp12_2)==0){
    printf("success\n");
  }else{
    printf("fail\n");
  }
  */
/*
  printf("\nkadai3 begin\n");
  while(1){
    Fp12_random(&Fp12_1);
    if(Fp12_legendre(&Fp12_1)==1) break;
  }
  printf("legendre ok!\n");
  Fp12_print(&Fp12_1);
  printf("aaaa\n");
  Fp12_sqrt(&Fp12_2,&Fp12_1);
  printf("aaaa\n");
  Fp12_print(&Fp12_2);
  Fp12_pow_ui(&Fp12_2,&Fp12_2,2);
  if(Fp12_cmp(&Fp12_1,&Fp12_2)==0){
    printf("(sqrt(Fp2_1))^2=Fp2_1:success\n");
  }else{
    printf("fail\n");
  }

*/
  /*kadai4 check*/
/*
  printf("\nkadai4 begin\n");


  mpz_set_ui(mpz_s,4);
  EFp12_random(&EFp12_1);
  EFp12_2.infinity=1;
  EFp12_set(&EFp12_3,&EFp12_1);

  for(i=0;i<4;i++){
    EFp12_ADD(&EFp12_2,&EFp12_2,&EFp12_1);
  }
  EFp12_print(&EFp12_2);
  for(i=0;i<2;i++){
    EFp12_DBL(&EFp12_3,&EFp12_3);
  }
  EFp12_print(&EFp12_3);
  EFp12_SCM(&EFp12_1,&EFp12_1,mpz_s);
  EFp12_print(&EFp12_1);
  //EFp12_random(&EFp12_2);
  if(EFp12_cmp(&EFp12_1,&EFp12_2)==0 && EFp12_cmp(&EFp12_1,&EFp12_3)==0){
    printf("p*4=p+p+p+p=p*2*2:success!!\n");
  }else{
    printf("false\n");
  }
*/
/*
  Fp6_random(&Fp6_1);
  Fp6_random(&Fp6_2);
  Fp6_inv(&Fp6_1,&Fp6_2);
  */
  /*
  EFp2_random(&EFp2_1);
  EFp2_random(&EFp2_2);
  clockstart();
  EFp2_SCM(&EFp2_2,&EFp2_1,EFp6_r);
  clockend();
  printf("------------\n");
  EFp2_print(&EFp2_2);

  EFp6_random(&EFp6_1);
  clockstart();
  EFp6_SCM(&EFp6_2,&EFp6_1,EFp6_r);
  clockend();
  EFp6_print(&EFp6_2);

  printf("random ok\n");
  EFp12_random(&EFp12_1);
  clockstart();
  EFp12_SCM(&EFp12_2,&EFp12_1,EFp12_r);
  clockend();
*/
  /*
  Fp12_random(&Fp12_1);
  Fp12_pow(&Fp12_3,&Fp12_1,BLS_p);
  Fp12_Frobenius(&Fp12_2,&Fp12_1);



  if(Fp_cmp(&(Fp12_3.x0.x0.x0),&(Fp12_2.x0.x0.x0))!=0) printf("w0:ng\n");
  if(Fp_cmp(&(Fp12_3.x0.x0.x1),&(Fp12_2.x0.x0.x1))!=0) printf("w1:ng\n");
  if(Fp_cmp(&(Fp12_3.x0.x1.x0),&(Fp12_2.x0.x1.x0))!=0) printf("w2:ng\n");
  if(Fp_cmp(&(Fp12_3.x0.x1.x1),&(Fp12_2.x0.x1.x1))!=0) printf("w3:ng\n");
  if(Fp_cmp(&(Fp12_3.x0.x2.x0),&(Fp12_2.x0.x2.x0))!=0) printf("w4:ng\n");
  if(Fp_cmp(&(Fp12_3.x0.x2.x1),&(Fp12_2.x0.x2.x1))!=0) printf("w5:ng\n");
  if(Fp_cmp(&(Fp12_3.x1.x0.x0),&(Fp12_2.x1.x0.x0))!=0) printf("w6:ng\n");
  if(Fp_cmp(&(Fp12_3.x1.x0.x1),&(Fp12_2.x1.x0.x1))!=0) printf("w7:ng\n");
  if(Fp_cmp(&(Fp12_3.x1.x1.x0),&(Fp12_2.x1.x1.x0))!=0) printf("w8:ng\n");
  if(Fp_cmp(&(Fp12_3.x1.x1.x1),&(Fp12_2.x1.x1.x1))!=0) printf("w9:ng\n");
  if(Fp_cmp(&(Fp12_3.x1.x2.x0),&(Fp12_2.x1.x2.x0))!=0) printf("w10:ng\n");
  if(Fp_cmp(&(Fp12_3.x1.x2.x1),&(Fp12_2.x1.x2.x1))!=0) printf("w11:ng\n");
  */

  return 0;
}


void mpz_print(mpz_t rop){
  mpz_out_str(stdout,10,rop);
  printf("\n");
}
/*Fp*/
void Fp_random(Fp *rop){
  mpz_urandomm(rop->x0, state, BLS_p);
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
void Fp_neg(Fp *rop,Fp *op){
  mpz_neg(rop->x0,op->x0);
  mpz_mod(rop->x0, rop->x0, BLS_p);
}
void Fp_pow(Fp *rop, Fp *op1, mpz_t scalar){
  mp_bitcnt_t i;
  Fp temp;
  Fp_init(&temp);
  Fp_set_ui(&temp,1);
  i=0;
  for(i=(mpz_sizeinbase(scalar,2)-1);;i--){
    Fp_mul(&temp,&temp,&temp);
    if(mpz_tstbit(scalar, i)==1){
      Fp_mul(&temp,&temp,op1);
    }
    if(i==0){
      break;
    }
  }
  Fp_set(rop,&temp);
  Fp_clear(&temp);
}
int Fp_legendre(Fp *op){
  mpz_t scalar;
  Fp temp;
  Fp_init(&temp);
  mpz_init(scalar);
  mpz_sub_ui(scalar,BLS_p,1);
  mpz_divexact_ui(scalar,scalar,2);
  Fp_pow(&temp,op,scalar);
  mpz_clear(scalar);
  if(Fp_cmp(&temp,&Fp_unity)==0||Fp_cmp(&temp,&Fp_zero)==0){
    Fp_clear(&temp);
    return 1;
  }else{
    Fp_clear(&temp);
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
  mpz_sub_ui(mpz_temp,BLS_p,1);
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
  mpz_powm_ui(rop->x0, base->x0, exp, BLS_p);
}
void Fp_div(Fp *rop, const Fp *op1, const Fp *op2){
  Fp inv;
  Fp_init(&inv);
  Fp_inv(&inv, op2);
  Fp_mul(rop, op1, &inv);
  Fp_clear(&inv);
}
void Fp_div_ui(Fp *rop, const Fp *op1, int op){
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

/*EFp*/
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
  EFp_clear(&temp);
}
void EFp_DBL(EFp *rop, EFp *op){
  if(op->infinity==1){
    rop->infinity = 1;
  }else if(Fp_cmp_ui(&(op->y),0)==0){
    rop->infinity = 1;
  }else{
    EFp temp;
    EFp_init(&temp);
    Fp slope;
    Fp_init(&slope);
    Fp slope2;
    Fp_init(&slope2);

    rop->infinity = 0;

    //傾き計算
    Fp_pow_ui(&(temp.x), &(op->x), 2);
    Fp_mul_ui(&(temp.x), &(temp.x), 3);

    Fp_mul_ui(&(temp.y), &(op->y), 2);
    Fp_div(&slope, &(temp.x), &(temp.y));

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
    Fp_clear(&slope);
    Fp_clear(&slope2);
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
    Fp_div(&slope, &(temp.y), &(temp.x));
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
int EFp_cmp(EFp *op1, EFp *op2){
  int temp,rop;
  temp = op1->infinity + op2->infinity;
  if(temp==2){
    rop=0;
  }else if(temp==1){
    rop=-1;
  }else{
    if(Fp_cmp(&(op1->x),&(op2->x))==0 && Fp_cmp(&(op1->y),&(op2->y))==0){
      rop=0;
    }else{
      rop=-1;
    }
  }
  return rop;
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
    printf("0");
  }
}
void EFp_set(EFp *rop, EFp *op){
  Fp_set(&(rop->x), &(op->x));
  Fp_set(&(rop->y), &(op->y));
  rop->infinity = op->infinity;
}

/*kadai*/
void fermat(){
  Fp a,p;
  Fp_init(&a);
  Fp_init(&p);

  Fp_random(&a);
  mpz_set(p.x0,BLS_p);
  Fp_sub_ui(&p,&p,1);
  Fp_pow(&a,&a,p.x0);
  printf("a^(p-1)=");
  Fp_print(&a);

  Fp_clear(&a);
  Fp_clear(&p);
}
void p_make(EFp *rop){
  EFp v;
  mpz_t temp;
  mpz_init(temp);
  EFp_init(&v);
  while(1){
    /*いかに計算式*/

    mpz_cdiv_q(temp,BLS_s,BLS_r);
    EFp_random(&v);
    EFp_SCM(&v,&v,temp);
    //EFp_print(&v);
    /*いかはテスト*/
    EFp_SCM(&v,&v,BLS_r);
    if(v.infinity==1) {
      break;
    }else{
      printf("not infinity\n");
    }
  }

  //printf("end!!!!");
  EFp_set(rop,&v);
  EFp_clear(&v);
  mpz_clear(temp);
}
void BLS_search(){
  EFp v;
  Fp temp;
  EFp_init(&v);
  Fp_init(&temp);
  Fp_set_ui(&BLS_b,0);
  while(1){
    Fp_add_ui(&BLS_b,&BLS_b,1);
    printf("start:BLS_b=");
    Fp_print(&BLS_b);
    EFp_random(&v);
    EFp_SCM(&v,&v,BLS_s);
    if(v.infinity==1){
      break;
    }
  }
  printf("end!!!!\n\n\nb = ");
  printf("BLS_b = ");
  Fp_print(&BLS_b);
  printf("BLS_s = ");
  mpz_print(BLS_s);
}
void Fp_testinit(){
  mpz_init(BLS_x);
  mpz_init(BLS_p);
  mpz_init(BLS_r);
  mpz_init(BLS_t);
  mpz_init(BLS_s);

  /*
  mpz_set_str(BLS_p, p_value, 10);
  mpz_set_str(BLS_r, r_value, 10);
  mpz_set_str(BLS_t, t_value, 10);
  mpz_set_str(BLS_s, s_value, 10);
  */
  Fp_init(&BLS_b);
  Fp_set_ui(&BLS_b, b_value);
  Fp_init(&BLS_a);
  Fp_set_ui(&BLS_a, a_value);

  Fp_init(&Fp_unity);
  Fp_set_ui(&Fp_unity,1);

  Fp_init(&Fp_zero);
  Fp_set_ui(&Fp_zero,0);

  gmp_randinit_mt(state);
  gmp_randseed_ui(state, (unsigned int) time(NULL));
}
void BLS_make(){
  mpz_t temp;
  mpz_init(temp);
  mpz_set_str(BLS_x,x_value,10);

  //BLS_p,BLS_rの設定
  BLS_pr(BLS_p,BLS_r);
  if( mpz_probab_prime_p(BLS_p,50)==0 && mpz_probab_prime_p(BLS_r,50)==0){
    printf("x:error\n");
    //return 0;
  }
  //BLS_tの設定
  mpz_add_ui(BLS_t,BLS_x,1);
  //BLS_sの設定
  mpz_set(BLS_s,BLS_p);
  mpz_add_ui(BLS_s,BLS_s,1);
  mpz_sub(BLS_s,BLS_s,BLS_t);
  /*
  printf("x = ");
  mpz_print(x);
  printf("p = ");
  mpz_print(BLS_p);
  printf("r = ");
  mpz_print(BLS_r);
  printf("t = ");
  mpz_print(BLS_t);
  printf("s = ");
  mpz_print(BLS_s);
  */
  mpz_clear(temp);
}
void BLS_make_random(){
  mpz_t x,test,length;
  mpz_init(test);
  mpz_init(length);
  mpz_init(x);

  //mpz_set_ui(x, 2);
  //mpz_pow_ui(x,x,70);

  mpz_set_ui(length, 2);
  mpz_pow_ui(length,length,70);

  while(1){
    mpz_urandomm(x,state,length);

    BLS_pr(BLS_p,BLS_r);
    if( mpz_probab_prime_p(BLS_p,50)!=0 && mpz_probab_prime_p(BLS_r,50)!=0){
      break;

      /*
      mpz_mod_ui(test,BLS_s,4);
      if(mpz_cmp_ui(test,0)==0){
        break;
      }
      */
    }
    //mpz_add_ui(x,x,1);
  }
  mpz_add_ui(BLS_t,x,1);

  mpz_set(BLS_s,BLS_p);
  mpz_add_ui(BLS_s,BLS_s,1);
  mpz_sub(BLS_s,BLS_s,BLS_t);



  printf("x = ");
  mpz_print(x);
  printf("p = ");
  mpz_print(BLS_p);
  printf("r = ");
  mpz_print(BLS_r);
  printf("t = ");
  mpz_print(BLS_t);
  printf("s = ");
  mpz_print(BLS_s);
  mpz_clear(x);
}
void BLS_pr(mpz_t p, mpz_t r){
  mpz_t temp,temp2;

  mpz_init(temp);
  mpz_init(temp2);
  int i;

  mpz_pow_ui(temp, BLS_x, 4);
  mpz_pow_ui(temp2, BLS_x, 2);
  mpz_sub(r, temp, temp2);
  mpz_add_ui(r, r, 1);

  mpz_sub_ui(temp, BLS_x, 1);
  mpz_pow_ui(temp, temp, 2);
  mpz_mul(temp, r, temp);
  mpz_set_ui(temp2, 3);
  if(mpz_divisible_p(temp, temp2)!=0){
    mpz_cdiv_q(p, temp, temp2);
    mpz_add(p, p, BLS_x);
  }else{
    mpz_set_ui(p,1);
  }
  mpz_clear(temp);
  mpz_clear(temp2);
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
void EFp2_random(EFp2 *rop){
  EFp2 temp;
  Fp2 ax;
  EFp2_init(&temp);
  Fp2_init(&ax);
  while(1){
    Fp2_random(&(rop->x));
    Fp2_pow_ui(&(temp.y),&(rop->x),3);
    Fp2_add(&(temp.y),&(temp.y),&Fp2_b);
    Fp2_mul(&ax,&(rop->x),&Fp2_a);
    Fp2_add(&(temp.y),&(temp.y),&ax);
    if(Fp2_legendre(&(temp.y))==1){
      break;
    }
  }
  Fp2_sqrt(&(rop->y),&(temp.y));
  rop->infinity = 0;
  EFp2_clear(&temp);
  Fp2_clear(&ax);
}
void EFp2_print(EFp2 *rop){
  if(rop->infinity==0){
    printf("(");
    Fp2_print(&(rop->x));
    printf(",");
    Fp2_print(&(rop->y));
    printf(")\n");
  }else{
    printf("infinity = %d\n",rop->infinity);
  }
}
void EFp2_set(EFp2 *rop, EFp2 *op){
  Fp2_set(&(rop->x), &(op->x));
  Fp2_set(&(rop->y), &(op->y));
  rop->infinity = op->infinity;
}
void EFp2_DBL(EFp2 *rop, EFp2 *op){
  if(op->infinity==1){
    rop->infinity = 1;
  }else if(Fp2_cmp(&(op->y),&Fp2_zero)==0){
    rop->infinity = 1;
  }else{
    EFp2 temp;
    EFp2_init(&temp);
    Fp2 slope;
    Fp2_init(&slope);
    Fp2 slope2;
    Fp2_init(&slope2);

    rop->infinity = 0;

    //傾き計算
    Fp2_pow_ui(&(temp.x), &(op->x), 2);
    Fp2_mul_ui(&(temp.x), &(temp.x), 3);

    Fp2_mul_ui(&(temp.y), &(op->y), 2);
    Fp2_div(&slope, &(temp.x), &(temp.y));

    //傾きの2乗計算
    Fp2_pow_ui(&slope2, &slope, 2);

    //xrの計算
    Fp2_mul_ui(&(temp.x), &(op->x), 2);
    Fp2_sub(&(temp.x), &slope2, &(temp.x));

    //yrの計算
    Fp2_sub(&(temp.y), &(op->x), &(temp.x));
    Fp2_mul(&(temp.y), &slope, &(temp.y));
    Fp2_sub(&(temp.y), &(temp.y), &(op->y));

    EFp2_set(rop,&temp);
    EFp2_clear(&temp);
    Fp2_clear(&slope);
    Fp2_clear(&slope2);
  }
}
void EFp2_ADD(EFp2 *rop, EFp2 *op1, EFp2 *op2){
  if((op1->infinity==1)&&(op2->infinity==1)){
    rop->infinity = 1;
  }else if(op1->infinity==1){
    EFp2_set(rop, op2);
    rop->infinity = 0;
  }else if(op2->infinity==1){
    EFp2_set(rop, op1);
    rop->infinity = 0;
  }else if(Fp2_cmp(&(op1->x), &(op2->x))==0){
    if(Fp2_cmp(&(op1->y), &(op2->y))==0){
      EFp2_DBL(rop, op1);
    }else{
      rop->infinity = 1;
    }
  }else{
    rop->infinity = 0;

    EFp2 temp;
    EFp2_init(&temp);
    Fp2 slope;
    Fp2_init(&slope);
    Fp2 slope2;
    Fp2_init(&slope2);

    //傾き計算
    Fp2_sub(&(temp.y), &(op1->y), &(op2->y));
    Fp2_sub(&(temp.x), &(op1->x), &(op2->x));
    Fp2_div(&slope, &(temp.y), &(temp.x));
    //傾きの2乗計
    Fp2_pow_ui(&slope2, &slope, 2);
    //xrの計算
    Fp2_add(&(temp.x), &(op1->x), &(op2->x));
    Fp2_sub(&(temp.x), &slope2, &(temp.x));
    //yrの計算
    Fp2_sub(&(temp.y), &(op1->x), &(temp.x));
    Fp2_mul(&(temp.y), &(temp.y), &slope);
    Fp2_sub(&(temp.y), &(temp.y), &(op1->y));

    EFp2_set(rop,&temp);
    Fp2_clear(&slope);
    Fp2_clear(&slope2);
    EFp2_clear(&temp);
  }
}
void EFp2_SCM(EFp2 *rop, EFp2 *op, mpz_t s){
  mp_bitcnt_t i;
  EFp2 temp,sum;
  EFp2_init(&temp);
  EFp2_init(&sum);
  EFp2_set(&temp,op);

  sum.infinity=1;
  for(i=(mpz_sizeinbase(s,2)-1);;i--){
    EFp2_DBL(&sum,&sum);
    if(mpz_tstbit(s,i)==1){
      EFp2_ADD(&sum,&sum,&temp);
    }
    if(i==0) break;
  }
  EFp2_set(rop,&sum);
  EFp2_clear(&sum);
  EFp2_clear(&temp);
}
int EFp2_cmp(EFp2 *op1, EFp2 *op2){
  int temp,rop;
  temp = op1->infinity + op2->infinity;
  if(temp==2){
    rop=0;
  }else if(temp==1){
    rop=-1;
  }else{
    if(Fp2_cmp(&(op1->x),&(op2->x))==0 && Fp2_cmp(&(op1->y),&(op2->y))==0){
      rop=0;
    }else{
      rop=-1;
    }
  }
  return rop;
}

void Fp2_make(){
  Fp_init(&Fp2_i);
  Fp_init(&Fp2_i2);
/*
  while(1){
    Fp_random(&Fp2_i);
    if(Fp_legendre(&Fp2_i)==-1){
      break;
    }
  }
*/
  Fp_sub_ui(&Fp2_i,&Fp_zero,1);
  if(Fp_legendre(&Fp2_i)==1){
    printf("FP2_error\n");
  }

  mpz_init(Fp2_p);
  mpz_pow_ui(Fp2_p,BLS_p,2);
  Fp2_init(&Fp2_unity);
  Fp2_init(&Fp2_zero);
  Fp_set_ui(&(Fp2_unity.x0),1);
  Fp_set_ui(&(Fp2_unity.x1),0);
  Fp_set_ui(&(Fp2_zero.x0),0);
  Fp_set_ui(&(Fp2_zero.x1),0);
  Fp2_init(&Fp2_b);
  Fp2_init(&Fp2_a);
  Fp_set(&(Fp2_b.x0),&BLS_b);
  Fp_set(&(Fp2_b.x1),&Fp_zero);
  Fp_set(&(Fp2_a.x0),&BLS_a);
  Fp_set(&(Fp2_a.x1),&Fp_zero);
}

int Fp2_legendre(Fp2 *op){
  mpz_t scalar;
  Fp2 temp;
  Fp2_init(&temp);
  mpz_init(scalar);
  mpz_sub_ui(scalar,Fp2_p,1);
  mpz_divexact_ui(scalar,scalar,2);
  Fp2_pow(&temp,op,scalar);
  mpz_clear(scalar);
  if(Fp2_cmp(&temp,&Fp2_unity)==0||Fp2_cmp(&temp,&Fp2_zero)==0){
    Fp2_clear(&temp);
    return 1;
  }else{
    Fp2_clear(&temp);
    return -1;
  }
}
int Fp2_legendre3(Fp2 *op){
  mpz_t scalar;
  Fp2 temp;
  Fp2_init(&temp);
  mpz_init(scalar);
  mpz_sub_ui(scalar,Fp2_p,1);
  mpz_divexact_ui(scalar,scalar,3);
  Fp2_pow(&temp,op,scalar);
  mpz_clear(scalar);
  if(Fp2_cmp(&temp,&Fp2_unity)==0||Fp2_cmp(&temp,&Fp2_zero)==0){
    Fp2_clear(&temp);
    return 1;
  }else{
    Fp2_clear(&temp);
    return -1;
  }
}
void Fp2_pow_ui(Fp2 *rop,Fp2 *op1,unsigned int op2){
  mpz_t temp;
  mpz_init(temp);
  mpz_set_ui(temp,op2);
  Fp2_pow(rop,op1,temp);
  mpz_clear(temp);
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
  mp_bitcnt_t i;
  Fp2 temp;
  Fp2_init(&temp);
  /*temp=unity*/
  Fp2_set(&temp,&Fp2_unity);

  i=0;
  for(i=(mpz_sizeinbase(scalar,2)-1);;i--){
    Fp2_mul(&temp,&temp,&temp);
    if(mpz_tstbit(scalar, i)==1){
      Fp2_mul(&temp,&temp,op1);
    }
    if(i==0){
      break;
    }
  }
  Fp2_set(rop,&temp);
  Fp2_clear(&temp);
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
  Fp2 temp,temp2;
  Fp2_init(&temp);
  Fp2_init(&temp2);

  Fp_mul(&(temp.x0),&(op1->x0),&(op2->x0));
  Fp_mul(&(temp2.x0),&(op1->x1),&(op2->x1));
  Fp_mul(&(temp2.x0),&(temp2.x0),&Fp2_i);

  Fp_mul(&(temp.x1),&(op1->x0),&(op2->x1));
  Fp_mul(&(temp2.x1),&(op1->x1),&(op2->x0));

  Fp2_add(rop,&temp,&temp2);
  Fp2_clear(&temp);
  Fp2_clear(&temp2);
}
void Fp2_conju(Fp2 *rop,Fp2 *op){
  Fp_set(&(rop->x0),&(op->x0));
  Fp_neg(&(rop->x1),&(op->x1));
}
void Fp2_inv(Fp2 *rop, Fp2 *op){
  Fp2 temp1,temp2;
  Fp2_init(&temp1);
  Fp2_init(&temp2);
  Fp2_conju(&temp1,op);
  //temp2=(a+bi)*(a-bi)
  Fp2_mul(&temp2,op,&temp1);
  Fp_div(&(rop->x0),&(temp1.x0),&(temp2.x0));
  Fp_div(&(rop->x1),&(temp1.x1),&(temp2.x0));
  Fp2_clear(&temp1);
  Fp2_clear(&temp2);
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

/*Fp6*/
void Fp6_conju(Fp6 *rop1,Fp6 *rop2,Fp6 *op){
  Fp2_set(&(rop1->x0),&(op->x0));
  Fp2_set(&(rop2->x0),&(op->x0));

  Fp2_mul(&(rop1->x1),&(op->x1),&Fp6_i2);
  Fp2_mul(&(rop2->x2),&(op->x2),&Fp6_i2);

  Fp2_mul(&(rop1->x2),&(op->x2),&Fp6_i3);
  Fp2_mul(&(rop2->x1),&(op->x1),&Fp6_i3);

}
void Fp6_make(){
  mpz_t mpz_temp;
  Fp2 Fp2_temp,Fp2_temp2,Fp2_temp3;
  mpz_init(mpz_temp);
  Fp2_init(&Fp2_temp);
  Fp2_init(&Fp2_temp2);
  Fp2_init(&Fp2_temp3);
  Fp2_init(&Fp6_i);
  Fp2_init(&Fp6_i2);
  Fp2_init(&Fp6_i3);
/*
  while(1){
    Fp2_random(&Fp6_i);
    if(Fp2_legendre3(&Fp6_i)==-1){
      break;
    }
  }
*/
  Fp_set_ui(&(Fp6_i.x0),1);
  Fp_set_ui(&(Fp6_i.x1),1);

  mpz_init(Fp6_p);
  mpz_pow_ui(Fp6_p,Fp2_p,3);
  Fp6_init(&Fp6_unity);
  Fp6_init(&Fp6_zero);
  Fp2_set(&(Fp6_unity.x0),&Fp2_unity);
  Fp2_set(&(Fp6_unity.x1),&Fp2_zero);
  Fp2_set(&(Fp6_unity.x2),&Fp2_zero);
  Fp2_set(&(Fp6_zero.x0),&Fp2_zero);
  Fp2_set(&(Fp6_zero.x1),&Fp2_zero);
  Fp2_set(&(Fp6_zero.x2),&Fp2_zero);

  Fp6_init(&Fp6_b);
  Fp6_init(&Fp6_a);
  Fp2_set(&(Fp6_b.x0),&Fp2_b);
  Fp2_set(&(Fp6_b.x1),&Fp2_zero);
  Fp2_set(&(Fp6_b.x2),&Fp2_zero);
  Fp2_set(&(Fp6_a.x0),&Fp2_a);
  Fp2_set(&(Fp6_a.x1),&Fp2_zero);
  Fp2_set(&(Fp6_b.x2),&Fp2_zero);

  //追加
  mpz_sub_ui(mpz_temp,Fp2_p,1);
  mpz_divexact_ui(mpz_temp,mpz_temp,3);
  while(1){
    while(1){
      Fp2_random(&Fp2_temp);
      if(Fp2_cmp(&Fp2_temp,&Fp2_zero)!=0) break;
    }
    Fp2_pow(&Fp2_temp2,&Fp6_i,mpz_temp);
    Fp2_pow_ui(&Fp2_temp3,&Fp2_temp2,2);
    if(Fp2_cmp(&Fp2_temp,&Fp2_unity)!=0 && Fp2_cmp(&Fp2_temp2,&Fp2_unity)!=0 && Fp2_cmp(&Fp2_temp3,&Fp2_unity)!=0) break;
  }
  /*check*/
  Fp2_set(&Fp6_i2,&Fp2_temp2);
  Fp2_set(&Fp6_i3,&Fp2_temp3);
  Fp2_clear(&Fp2_temp);
  Fp2_clear(&Fp2_temp2);
  Fp2_clear(&Fp2_temp3);
  mpz_clear(mpz_temp);
}
int Fp6_legendre(Fp6 *op){
  mpz_t scalar;
  Fp6 temp;
  Fp6_init(&temp);
  mpz_init(scalar);
  mpz_sub_ui(scalar,Fp6_p,1);
  mpz_divexact_ui(scalar,scalar,2);
  Fp6_pow(&temp,op,scalar);
  mpz_clear(scalar);
  if(Fp6_cmp(&temp,&Fp6_unity)==0||Fp6_cmp(&temp,&Fp6_zero)==0){
    Fp6_clear(&temp);
    return 1;
  }else{
    Fp6_clear(&temp);
    return -1;
  }
}
void Fp6_pow_ui(Fp6 *rop,Fp6 *op1,unsigned int op2){
  mpz_t temp;
  mpz_init(temp);
  mpz_set_ui(temp,op2);
  Fp6_pow(rop,op1,temp);
  mpz_clear(temp);
}
void Fp6_mul_ui(Fp6 *rop,Fp6 *op1,unsigned int op2){/*スカラー倍算*/
  Fp2_mul_ui(&(rop->x0),&(op1->x0),op2);
  Fp2_mul_ui(&(rop->x1),&(op1->x1),op2);
  Fp2_mul_ui(&(rop->x2),&(op1->x2),op2);
}
void Fp6_sqrt(Fp6 *rop, Fp6 *a){
  mpz_t q,mpz_temp,mpz_temp2;
  Fp6 n,x,y,b,t,Fp6_temp,Fp6_temp2;
  unsigned int e,r,i,m;
  Fp6_init(&n);
  Fp6_init(&x);
  Fp6_init(&b);
  Fp6_init(&y);
  Fp6_init(&t);
  Fp6_init(&Fp6_temp);
  Fp6_init(&Fp6_temp2);

  mpz_init(q);
  mpz_init(mpz_temp);
  mpz_init(mpz_temp2);

  //1
  while(1){
    Fp6_random(&n);
    if(Fp6_legendre(&n)==-1){
      break;
    }
  }
  //2
  mpz_sub_ui(mpz_temp,Fp6_p,1);
  e = mpz_scan1(mpz_temp,0);
  mpz_fdiv_q_2exp(q,mpz_temp,e);
  //3
  Fp6_pow(&y,&n,q);
  r=e;
  mpz_sub_ui(mpz_temp,q,1);
  mpz_divexact_ui(mpz_temp,mpz_temp,2);
  Fp6_pow(&x,a,mpz_temp);
  //4
  Fp6_pow_ui(&b,&x,2);
  Fp6_mul(&b,&b,a);
  Fp6_mul(&x,a,&x);
  //5
  while(Fp6_cmp(&b,&Fp6_unity)!=0){
    m=0;
    Fp6_set(&Fp6_temp,&b);
    while(Fp6_cmp(&Fp6_temp,&Fp6_unity)!=0){
      Fp6_pow_ui(&Fp6_temp,&Fp6_temp,2);
      m++;
    }
    r=r-m-1;
    mpz_set_ui(mpz_temp,2);
    mpz_pow_ui(mpz_temp,mpz_temp,r);
    Fp6_pow(&t,&y,mpz_temp);
    //ok
    Fp6_pow_ui(&y,&t,2);
    r=m;
    Fp6_mul(&x,&x,&t);
    Fp6_mul(&b,&b,&y);
  }
  Fp6_set(rop,&x);
  mpz_clear(mpz_temp);
  Fp6_clear(&Fp6_temp2);
  Fp6_clear(&Fp6_temp);
  Fp6_clear(&n);
  mpz_clear(q);
  Fp6_clear(&x);
  Fp6_clear(&b);
  Fp6_clear(&y);
  Fp6_clear(&t);
}
void Fp6_pow(Fp6 *rop, Fp6 *op1, mpz_t scalar){
  mp_bitcnt_t i;
  Fp6 temp;
  Fp6_init(&temp);
  /*temp=unity*/
  Fp6_set(&temp,&Fp6_unity);

  i=0;
  for(i=(mpz_sizeinbase(scalar,2)-1);;i--){
    Fp6_mul(&temp,&temp,&temp);
    if(mpz_tstbit(scalar, i)==1){
      Fp6_mul(&temp,&temp,op1);
    }
    if(i==0){
      break;
    }
  }
  Fp6_set(rop,&temp);
  Fp6_clear(&temp);
}
void Fp6_init(Fp6 *op){
  Fp2_init(&(op->x0));
  Fp2_init(&(op->x1));
  Fp2_init(&(op->x2));
}
void Fp6_clear(Fp6 *op){
  Fp2_clear(&(op->x0));
  Fp2_clear(&(op->x1));
  Fp2_clear(&(op->x2));
}
void Fp6_random(Fp6 *rop){
  Fp2_random(&(rop->x0));
  Fp2_random(&(rop->x1));
  Fp2_random(&(rop->x2));
}
void Fp6_add(Fp6 *rop, Fp6 *op1, Fp6 *op2){
  Fp2_add(&(rop->x0),&(op1->x0),&(op2->x0));
  Fp2_add(&(rop->x1),&(op1->x1),&(op2->x1));
  Fp2_add(&(rop->x2),&(op1->x2),&(op2->x2));
}
void Fp6_sub(Fp6 *rop, Fp6 *op1, Fp6 *op2){
  Fp2_sub(&(rop->x0),&(op1->x0),&(op2->x0));
  Fp2_sub(&(rop->x1),&(op1->x1),&(op2->x1));
  Fp2_sub(&(rop->x2),&(op1->x2),&(op2->x2));
}
void Fp6_mul(Fp6 *rop, Fp6 *op1, Fp6 *op2){
  Fp6 temp1,temp2,temp3;
  Fp6_init(&temp1);
  Fp6_init(&temp2);
  Fp6_init(&temp3);
  Fp2_mul(&(temp1.x0),&(op1->x0),&(op2->x0));
  Fp2_mul(&(temp1.x1),&(op1->x0),&(op2->x1));
  Fp2_mul(&(temp1.x2),&(op1->x0),&(op2->x2));

  Fp2_mul(&(temp2.x0),&(op1->x1),&(op2->x2));
  Fp2_mul(&(temp2.x0),&(temp2.x0),&Fp6_i);
  Fp2_mul(&(temp2.x1),&(op1->x1),&(op2->x0));
  Fp2_mul(&(temp2.x2),&(op1->x1),&(op2->x1));

  Fp2_mul(&(temp3.x0),&(op1->x2),&(op2->x1));
  Fp2_mul(&(temp3.x0),&(temp3.x0),&Fp6_i);
  Fp2_mul(&(temp3.x1),&(op1->x2),&(op2->x2));
  Fp2_mul(&(temp3.x1),&(temp3.x1),&Fp6_i);
  Fp2_mul(&(temp3.x2),&(op1->x2),&(op2->x0));

  Fp6_add(&temp1,&temp1,&temp2);
  Fp6_add(rop,&temp1,&temp3);
  Fp6_clear(&temp1);
  Fp6_clear(&temp2);
  Fp6_clear(&temp3);
}

void Fp6_inv(Fp6 *rop, Fp6 *op){
  Fp6 temp1,temp2;
  Fp6_init(&temp1);
  Fp6_init(&temp2);
  Fp6_conju(&temp1,&temp2,op);
  //temp1=(a+bi)*(a+bi)'*(a+bi)''
  //temp2=(a+bi)'*(a+bi)''
  Fp6_mul(&temp1,&temp1,&temp2);
  Fp6_mul(&temp2,op,&temp1);

  Fp2_inv(&(temp2.x0),&(temp2.x0));
  Fp2_mul(&(rop->x0),&(temp1.x0),&(temp2.x0));
  Fp2_mul(&(rop->x1),&(temp1.x1),&(temp2.x0));
  Fp2_mul(&(rop->x2),&(temp1.x2),&(temp2.x0));

  Fp6_clear(&temp1);
  Fp6_clear(&temp2);
}
void Fp6_div(Fp6 *rop,Fp6 *op1,Fp6 *op2){
  Fp6 inv;
  Fp6_init(&inv);
  Fp6_inv(&inv,op2);
  Fp6_mul(rop,op1,&inv);
  Fp6_clear(&inv);
}
int Fp6_cmp(Fp6 *op1,Fp6 *op2){
  if(Fp2_cmp(&(op1->x0),&(op2->x0))==0 && Fp2_cmp(&(op1->x1),&(op2->x1))==0 && Fp2_cmp(&(op1->x2),&(op2->x2))==0){
    return 0;
  }else{
    return 1;
  }
}
void Fp6_print(Fp6 *rop){
  printf("(");
  mpz_out_str(stdout,10,rop->x0.x0.x0);
  printf(",");
  mpz_out_str(stdout,10,rop->x0.x1.x0);
  printf(",");
  mpz_out_str(stdout,10,rop->x1.x0.x0);
  printf(",");
  mpz_out_str(stdout,10,rop->x1.x1.x0);
  printf(",");
  mpz_out_str(stdout,10,rop->x2.x0.x0);
  printf(",");
  mpz_out_str(stdout,10,rop->x2.x1.x0);
  printf(")\n");
}
void Fp6_set(Fp6 *rop,Fp6 *op){
  Fp2_set(&(rop->x0),&(op->x0));
  Fp2_set(&(rop->x1),&(op->x1));
  Fp2_set(&(rop->x2),&(op->x2));
}
void EFp6_init(EFp6 *rop){
  Fp6_init(&(rop->x));
  Fp6_init(&(rop->y));
  rop->infinity = 0;
}
void EFp6_clear(EFp6 *rop){
  Fp6_clear(&(rop->x));
  Fp6_clear(&(rop->y));
}
void EFp6_random(EFp6 *rop){
  EFp6 temp;
  Fp6 ax;
  EFp6_init(&temp);
  Fp6_init(&ax);
  while(1){
    Fp6_random(&(rop->x));
    Fp6_pow_ui(&(temp.y),&(rop->x),3);
    Fp6_add(&(temp.y),&(temp.y),&Fp6_b);
    Fp6_mul(&ax,&(rop->x),&Fp6_a);
    Fp6_add(&(temp.y),&(temp.y),&ax);
    if(Fp6_legendre(&(temp.y))==1){
      break;
    }
  }
  Fp6_sqrt(&(rop->y),&(temp.y));
  rop->infinity = 0;
  EFp6_clear(&temp);
  Fp6_clear(&ax);
}
void EFp6_print(EFp6 *rop){
  if(rop->infinity==0){
    printf("(");
    Fp6_print(&(rop->x));
    printf(",");
    Fp6_print(&(rop->y));
    printf(")\n");
  }else{
    printf("infinity = %d\n",rop->infinity);
  }
}
void EFp6_set(EFp6 *rop, EFp6 *op){
  Fp6_set(&(rop->x), &(op->x));
  Fp6_set(&(rop->y), &(op->y));
  rop->infinity = op->infinity;
}
void EFp6_DBL(EFp6 *rop, EFp6 *op){
  if(op->infinity==1){
    rop->infinity = 1;
  }else if(Fp6_cmp(&(op->y),&Fp6_zero)==0){
    rop->infinity = 1;
  }else{
    EFp6 temp;
    EFp6_init(&temp);
    Fp6 slope;
    Fp6_init(&slope);
    Fp6 slope2;
    Fp6_init(&slope2);

    rop->infinity = 0;

    //傾き計算
    Fp6_pow_ui(&(temp.x), &(op->x), 2);
    Fp6_mul_ui(&(temp.x), &(temp.x), 3);

    Fp6_mul_ui(&(temp.y), &(op->y), 2);
    Fp6_div(&slope, &(temp.x), &(temp.y));

    //傾きの2乗計算
    Fp6_pow_ui(&slope2, &slope, 2);

    //xrの計算
    Fp6_mul_ui(&(temp.x), &(op->x), 2);
    Fp6_sub(&(temp.x), &slope2, &(temp.x));

    //yrの計算
    Fp6_sub(&(temp.y), &(op->x), &(temp.x));
    Fp6_mul(&(temp.y), &slope, &(temp.y));
    Fp6_sub(&(temp.y), &(temp.y), &(op->y));

    EFp6_set(rop,&temp);
    EFp6_clear(&temp);
    Fp6_clear(&slope);
    Fp6_clear(&slope2);
  }
}
void EFp6_ADD(EFp6 *rop, EFp6 *op1, EFp6 *op2){
  if((op1->infinity==1)&&(op2->infinity==1)){
    rop->infinity = 1;
  }else if(op1->infinity==1){
    EFp6_set(rop, op2);
    rop->infinity = 0;
  }else if(op2->infinity==1){
    EFp6_set(rop, op1);
    rop->infinity = 0;
  }else if(Fp6_cmp(&(op1->x), &(op2->x))==0){
    if(Fp6_cmp(&(op1->y), &(op2->y))==0){
      EFp6_DBL(rop, op1);
    }else{
      rop->infinity = 1;
    }
  }else{
    rop->infinity = 0;

    EFp6 temp;
    EFp6_init(&temp);
    Fp6 slope;
    Fp6_init(&slope);
    Fp6 slope2;
    Fp6_init(&slope2);

    //傾き計算
    Fp6_sub(&(temp.y), &(op1->y), &(op2->y));
    Fp6_sub(&(temp.x), &(op1->x), &(op2->x));
    Fp6_div(&slope, &(temp.y), &(temp.x));
    //傾きの2乗計
    Fp6_pow_ui(&slope2, &slope, 2);
    //xrの計算
    Fp6_add(&(temp.x), &(op1->x), &(op2->x));
    Fp6_sub(&(temp.x), &slope2, &(temp.x));
    //yrの計算
    Fp6_sub(&(temp.y), &(op1->x), &(temp.x));
    Fp6_mul(&(temp.y), &(temp.y), &slope);
    Fp6_sub(&(temp.y), &(temp.y), &(op1->y));

    EFp6_set(rop,&temp);
    Fp6_clear(&slope);
    Fp6_clear(&slope2);
    EFp6_clear(&temp);
  }
}
void EFp6_SCM(EFp6 *rop, EFp6 *op, mpz_t s){
  mp_bitcnt_t i;
  EFp6 temp,sum;
  EFp6_init(&temp);
  EFp6_init(&sum);
  EFp6_set(&temp,op);

  sum.infinity=1;
  for(i=(mpz_sizeinbase(s,2)-1);;i--){
    EFp6_DBL(&sum,&sum);
    if(mpz_tstbit(s,i)==1){
      EFp6_ADD(&sum,&sum,&temp);
    }
    if(i==0) break;
  }
  EFp6_set(rop,&sum);
  EFp6_clear(&sum);
  EFp6_clear(&temp);
}
int EFp6_cmp(EFp6 *op1, EFp6 *op2){
  int temp,rop;
  temp = op1->infinity + op2->infinity;
  if(temp==2){
    rop=0;
  }else if(temp==1){
    rop=-1;
  }else{
    if(Fp6_cmp(&(op1->x),&(op2->x))==0 && Fp6_cmp(&(op1->y),&(op2->y))==0){
      rop=0;
    }else{
      rop=-1;
    }
  }
  return rop;
}
void EFp12_init(EFp12 *rop){
  Fp12_init(&(rop->x));
  Fp12_init(&(rop->y));
  rop->infinity = 0;
}
void EFp12_clear(EFp12 *rop){
  Fp12_clear(&(rop->x));
  Fp12_clear(&(rop->y));
}
void EFp12_random(EFp12 *rop){
  EFp12 temp;
  Fp12 ax;
  EFp12_init(&temp);
  Fp12_init(&ax);
  while(1){
    Fp12_random(&(rop->x));
    Fp12_pow_ui(&(temp.y),&(rop->x),3);
    Fp12_add(&(temp.y),&(temp.y),&Fp12_b);
    Fp12_mul(&ax,&(rop->x),&Fp12_a);
    Fp12_add(&(temp.y),&(temp.y),&ax);
    if(Fp12_legendre(&(temp.y))==1){
      break;
    }
  }
  Fp12_sqrt(&(rop->y),&(temp.y));
  rop->infinity = 0;
  EFp12_clear(&temp);
  Fp12_clear(&ax);
}
void EFp12_print(EFp12 *rop){
  if(rop->infinity==0){
    printf("(");
    Fp12_print(&(rop->x));
    printf(",");
    Fp12_print(&(rop->y));
    printf(")\n");
  }else{
    printf("infinity = %d\n",rop->infinity);
  }
}
void EFp12_set(EFp12 *rop, EFp12 *op){
  Fp12_set(&(rop->x), &(op->x));
  Fp12_set(&(rop->y), &(op->y));
  rop->infinity = op->infinity;
}
void EFp12_DBL(EFp12 *rop, EFp12 *op){
  if(op->infinity==1){
    rop->infinity = 1;
  }else if(Fp12_cmp(&(op->y),&Fp12_zero)==0){
    rop->infinity = 1;
  }else{
    EFp12 temp;
    EFp12_init(&temp);
    Fp12 slope;
    Fp12_init(&slope);
    Fp12 slope2;
    Fp12_init(&slope2);

    rop->infinity = 0;

    //傾き計算
    Fp12_pow_ui(&(temp.x), &(op->x), 2);
    Fp12_mul_ui(&(temp.x), &(temp.x), 3);

    Fp12_mul_ui(&(temp.y), &(op->y), 2);
    Fp12_div(&slope, &(temp.x), &(temp.y));

    //傾きの2乗計算
    Fp12_pow_ui(&slope2, &slope, 2);

    //xrの計算
    Fp12_mul_ui(&(temp.x), &(op->x), 2);
    Fp12_sub(&(temp.x), &slope2, &(temp.x));

    //yrの計算
    Fp12_sub(&(temp.y), &(op->x), &(temp.x));
    Fp12_mul(&(temp.y), &slope, &(temp.y));
    Fp12_sub(&(temp.y), &(temp.y), &(op->y));

    EFp12_set(rop,&temp);
    EFp12_clear(&temp);
    Fp12_clear(&slope);
    Fp12_clear(&slope2);
  }
}
void Fp12_Frobenius(Fp12 *rop,Fp12 *op){
  //w0
  Fp_set(&(rop->x0.x0.x0),&(op->x0.x0.x0));
  //w1
  Fp_neg(&(rop->x0.x0.x1),&(op->x0.x0.x1));
  //w2
  Fp_mul(&(rop->x0.x1.x0),&(op->x0.x1.x0),&c6);
  //w3
  Fp_mul(&(rop->x0.x1.x1),&(op->x0.x1.x1),&c6);
  Fp_neg(&(rop->x0.x1.x1),&(rop->x0.x1.x1));
  //w4
  Fp_mul(&(rop->x0.x2.x0),&(op->x0.x2.x0),&c3);
  //w5
  Fp_mul(&(rop->x0.x2.x1),&(op->x0.x1.x1),&c3);
  Fp_neg(&(rop->x0.x1.x1),&(rop->x0.x1.x1));

  //w6
  //Fp_mul(&(rop->x1.x0.x0),&(op->x1.x0.x0),&c12);
  Fp_pow(&(rop->x1.x0.x0),&(op->x1.x0.x0),BLS_p);
  //w7
  //Fp_mul(&(rop->x1.x0.x1),&(op->x1.x0.x1),&c12);
  //Fp_neg(&(rop->x1.x0.x1),&(rop->x1.x0.x1));
  Fp_pow(&(rop->x1.x0.x1),&(op->x1.x0.x1),BLS_p);
  //w8
  Fp_mul(&(rop->x1.x1.x0),&(op->x1.x1.x0),&c4);
  //w9
  Fp_mul(&(rop->x1.x1.x1),&(op->x1.x1.x1),&c4);
  Fp_neg(&(rop->x1.x1.x1),&(rop->x1.x1.x1));

  //w10
  //Fp_mul(&(rop->x1.x2.x0),&(op->x1.x2.x0),&c5);
  Fp_pow(&(rop->x1.x2.x0),&(op->x1.x2.x0),BLS_p);
  //w5
  //Fp_mul(&(rop->x1.x2.x1),&(op->x1.x2.x1),&c5);
  //Fp_neg(&(rop->x1.x2.x1),&(rop->x1.x2.x1));
  Fp_pow(&(rop->x1.x2.x1),&(op->x1.x2.x1),BLS_p);
}

void Fp12_Frobenius_init(){
  mpz_t temp1;
  Fp_init(&c4);
  Fp_init(&c12);
  Fp_init(&c3);
  Fp_init(&c6);
  Fp_init(&c5);
  mpz_init(temp1);
  mpz_sub_ui(temp1,BLS_p,1);
  mpz_div_ui(temp1,temp1,6);
  mpz_print(temp1);
  Fp_pow(&c6,&Fp2_i,temp1);
  Fp_pow_ui(&c3,&c6,2);
  mpz_sub_ui(temp1,BLS_p,1);
  mpz_div_ui(temp1,temp1,4);
  Fp_pow_ui(&c4,&c12,3);

  //Fp_pow_ui(&c5,&c4,4);

  Fp_print(&c6);
  Fp_print(&c3);
  Fp_print(&c12);
  Fp_print(&c4);
  Fp_print(&c5);

  mpz_clear(temp1);
}
void EFp12_ADD(EFp12 *rop, EFp12 *op1, EFp12 *op2){
  if((op1->infinity==1)&&(op2->infinity==1)){
    rop->infinity = 1;
  }else if(op1->infinity==1){
    EFp12_set(rop, op2);
    rop->infinity = 0;
  }else if(op2->infinity==1){
    EFp12_set(rop, op1);
    rop->infinity = 0;
  }else if(Fp12_cmp(&(op1->x), &(op2->x))==0){
    if(Fp12_cmp(&(op1->y), &(op2->y))==0){
      EFp12_DBL(rop, op1);
    }else{
      rop->infinity = 1;
    }
  }else{
    rop->infinity = 0;

    EFp12 temp;
    EFp12_init(&temp);
    Fp12 slope;
    Fp12_init(&slope);
    Fp12 slope2;
    Fp12_init(&slope2);

    //傾き計算
    Fp12_sub(&(temp.y), &(op1->y), &(op2->y));
    Fp12_sub(&(temp.x), &(op1->x), &(op2->x));
    Fp12_div(&slope, &(temp.y), &(temp.x));
    //傾きの2乗計
    Fp12_pow_ui(&slope2, &slope, 2);
    //xrの計算
    Fp12_add(&(temp.x), &(op1->x), &(op2->x));
    Fp12_sub(&(temp.x), &slope2, &(temp.x));
    //yrの計算
    Fp12_sub(&(temp.y), &(op1->x), &(temp.x));
    Fp12_mul(&(temp.y), &(temp.y), &slope);
    Fp12_sub(&(temp.y), &(temp.y), &(op1->y));

    EFp12_set(rop,&temp);
    Fp12_clear(&slope);
    Fp12_clear(&slope2);
    EFp12_clear(&temp);
  }
}
void EFp12_inv(EFp12 *rop,EFp12 *op){
  Fp12_set(&(rop->x),&(op->x));
  Fp12_sub(&(rop->y),&Fp12_zero,&(op->y));
}
void EFp12_sub(EFp12 *rop,EFp12 *op1,EFp12 *op2){
  EFp12 inv;
  EFp12_init(&inv);
  EFp12_inv(&inv,op2);
  EFp12_ADD(rop,op1,&inv);
  EFp12_clear(&inv);
}
void EFp12_SCM(EFp12 *rop, EFp12 *op, mpz_t s){
  mp_bitcnt_t i;
  EFp12 temp,sum;
  EFp12_init(&temp);
  EFp12_init(&sum);
  EFp12_set(&temp,op);

  sum.infinity=1;
  for(i=(mpz_sizeinbase(s,2)-1);;i--){
    EFp12_DBL(&sum,&sum);
    if(mpz_tstbit(s,i)==1){
      EFp12_ADD(&sum,&sum,&temp);
    }
    if(i==0) break;
  }
  EFp12_set(rop,&sum);
  EFp12_clear(&sum);
  EFp12_clear(&temp);
}
int EFp12_cmp(EFp12 *op1, EFp12 *op2){
  int temp,rop;
  temp = op1->infinity + op2->infinity;
  if(temp==2){
    rop=0;
  }else if(temp==1){
    rop=-1;
  }else{
    if(Fp12_cmp(&(op1->x),&(op2->x))==0 && Fp12_cmp(&(op1->y),&(op2->y))==0){
      rop=0;
    }else{
      rop=-1;
    }
  }
  return rop;
}
