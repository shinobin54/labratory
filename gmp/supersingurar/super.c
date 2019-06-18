#include <stdio.h>
#include <gmp.h>
#include <time.h>
#include <unistd.h>
#define p_value "139"
/*
#define x_value ""

#define r_value ""
#define t_value ""
#define s_value ""
#define b_value -13
#define a_value -7
#define i2_value -4/*i^2(i <\- Fp)=-4
*/
/*EFp*/


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
mpz_t Fp_p,x,EFp_r,EFp_t,s;
Fp EFp_a,EFp_b,curve_a,curve_b;

void example_init();
void mpz_print(mpz_t rop);
void Fp_init(Fp *rop);
void Fp_printf(Fp *rop,char *str);
void Fp_clear(Fp *rop);
void Fp_set_random(Fp *rop);
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
void Fp_pow(Fp *rop, Fp *op1,mpz_t scalar);void Fp_init(Fp *rop);
void Fp_set(Fp *rop, Fp *op);
void Fp_set_ui(Fp *rop, int op);
void Fp_clear(Fp *rop);
void Fp_add(Fp *rop, const Fp *op1, const Fp *op2);
void Fp_sub(Fp *rop, const Fp *op1, const Fp *op2);
void Fp_mul(Fp *rop, const Fp *op1, const Fp *op2);
void Fp_inv(Fp *rop, const Fp *op1);
void Fp_add_ui(Fp *rop, const Fp *op1, int op2);
void Fp_mul_ui (Fp *rop, const Fp *op1, int ui);
void Fp_sub_ui (Fp *rop, const Fp *op1, int ui);
int Fp_cmp (const Fp *op1, const Fp *op2);
int Fp_cmp_ui (const Fp *op1, int op2);
void Fp_set_neg(Fp *rop,Fp *op);
void Fp_set_neg_ui(Fp *rop,unsigned int UI);
void Fp_pow_gmp(Fp *rop, const Fp *base, const Fp *exp);
void Fp_pow_ui(Fp *rop, const Fp *base, int exp);
void Fp_div(Fp *rop, const Fp *op1, const Fp *op2);
void Fp_div_ui(Fp *rop, const Fp *op1, int op);
int Fp_legendre(Fp *op);
void Fp_sqrt(Fp *rop,Fp *a);
void Fp_pow(Fp *rop, Fp *op1,mpz_t scalar);
int Fp_cmp_zero(Fp *op);
//EFp
void EFp_rational_point(EFp *P);
void EFp_init(EFp *P);
void EFp_set(EFp *ANS,EFp *A);
void EFp_set_ui(EFp *ANS,unsigned long int UI);
void EFp_set_neg(EFp *ANS,EFp *A);
int  EFp_cmp(EFp *A,EFp *B);
void EFp_ECD(EFp *ANS,EFp *P);
void EFp_ECA(EFp *ANS,EFp *P1,EFp *P2);
void EFp_SCM(EFp *ANS,EFp *P,mpz_t scalar);
void EFp_printf(EFp *P,char *str);

void example_init(){
  /*Fp*/
	mpz_init(Fp_p);
	mpz_init(x);
	mpz_init(EFp_r);
  Fp_init(&curve_a);
  Fp_init(&curve_b);
  /*set Fp_p*/
	mpz_set_str(Fp_p,p_value,10);
  /*set EFp_r*/
  mpz_sub_ui(EFp_r,Fp_p,1);
  Fp_set_neg_ui(&curve_a,7);
  Fp_set_neg_ui(&curve_b,13);

  gmp_randinit_mt(state);
  gmp_randseed_ui(state, (unsigned int) time(NULL));
}
void mpz_print(mpz_t rop){
  mpz_out_str(stdout,10,rop);
  printf("\n");
}
/*Fp*/
void Fp_set_random(Fp *rop){
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
void Fp_add(Fp *rop, const Fp *op1, const Fp *op2){
  mpz_add(rop->x0, op1->x0, op2->x0);
  mpz_mod(rop->x0, rop->x0, Fp_p);
}
void Fp_add_ui(Fp *rop, const Fp *op1, int op2){
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
void Fp_mul(Fp *rop, const Fp *op1, const Fp *op2){
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
void Fp_set_neg(Fp *rop,Fp *op){
  mpz_neg(rop->x0,op->x0);
  mpz_mod(rop->x0, rop->x0, Fp_p);
}
void Fp_set_neg_ui(Fp *rop,unsigned int UI){
  mpz_sub_ui(rop->x0,Fp_p,UI);
  mpz_mod(rop->x0, rop->x0, Fp_p);
}
void Fp_pow(Fp *rop, Fp *op1, mpz_t scalar){
  int i;
  Fp temp;
  Fp_init(&temp);
  Fp_set_ui(&temp,1);

  for(i=(int)(mpz_sizeinbase(scalar,2)-1);i>=0;i--){
    Fp_mul(&temp,&temp,&temp);
    if(mpz_tstbit(scalar, i)==1){
      Fp_mul(&temp,&temp,op1);
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
  mpz_sub_ui(scalar,Fp_p,1);
  mpz_divexact_ui(scalar,scalar,2);
  Fp_pow(&temp,op,scalar);
  mpz_clear(scalar);
  if(Fp_cmp_ui(&temp,1)==0||Fp_cmp_ui(&temp,0)==0){
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
    Fp_set_random(&n);
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
void Fp_printf(Fp *rop,char *str){
  printf("%s",str);
  mpz_out_str(stdout,10,rop->x0);
  printf("\n");
}
int Fp_cmp_zero(Fp *op){
  return Fp_cmp_ui(op,0)!=0;
}
/*EFp*/


void EFp_rational_point(EFp *P){
    Fp tmp1,tmp2,tmp_x;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    Fp_init(&tmp_x);

    while(1){
        Fp_set_random(&P->x);
        Fp_mul(&tmp1,&P->x,&P->x);
        Fp_mul(&tmp2,&tmp1,&P->x);
        Fp_mul(&tmp2,&tmp2,&curve_a);
        Fp_add(&tmp_x,&tmp2,&curve_b);
        if(Fp_legendre(&tmp_x)==1){
            Fp_sqrt(&P->y,&tmp_x);
            break;
        }
    }
}
void EFp_init(EFp *P){
    Fp_init(&P->x);
    Fp_init(&P->y);
    P->infinity=0;
}
void EFp_set(EFp *ANS,EFp *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}
void EFp_set_ui(EFp *ANS,unsigned long int UI){
    Fp_set_ui(&ANS->x,UI);
    Fp_set_ui(&ANS->y,UI);
    ANS->infinity=0;
}
void EFp_set_neg(EFp *ANS,EFp *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set_neg(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}
int  EFp_cmp(EFp *A,EFp *B){
    if(Fp_cmp(&A->x,&B->x)==0 && Fp_cmp(&A->y,&B->y)==0){
        return 0;
    }else if(A->infinity==1&&B->infinity==1){
	return 0;
    }else{
    return 1;
    }
}
void EFp_ECD(EFp *ANS,EFp *P){
    static EFp tmp1_EFp;
    static Fp tmp1_Fp,tmp2_Fp,tmp3_Fp;
    if(Fp_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    EFp_set(&tmp1_EFp,P);

    Fp_add(&tmp1_Fp,&tmp1_EFp.y,&tmp1_EFp.y);
    Fp_inv(&tmp1_Fp,&tmp1_Fp);


    Fp_mul(&tmp2_Fp,&tmp1_EFp.x,&tmp1_EFp.x);
    Fp_add(&tmp3_Fp,&tmp2_Fp,&tmp2_Fp);
    Fp_add(&tmp2_Fp,&tmp2_Fp,&tmp3_Fp);

    Fp_mul(&tmp3_Fp,&tmp1_Fp,&tmp2_Fp);
    Fp_mul(&tmp1_Fp,&tmp3_Fp,&tmp3_Fp);

    Fp_add(&tmp2_Fp,&tmp1_EFp.x,&tmp1_EFp.x);
    Fp_sub(&ANS->x,&tmp1_Fp,&tmp2_Fp);

    Fp_sub(&tmp1_Fp,&tmp1_EFp.x,&ANS->x);
    Fp_mul(&tmp2_Fp,&tmp3_Fp,&tmp1_Fp);
    Fp_sub(&ANS->y,&tmp2_Fp,&tmp1_EFp.y);
}
void EFp_ECA(EFp *ANS,EFp *P1,EFp *P2){
    static EFp tmp1_EFp,tmp2_EFp;
    static Fp tmp1_Fp,tmp2_Fp,tmp3_Fp;
    if(P1->infinity==1){
        EFp_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp_set(ANS,P1);
        return;
    }else if(Fp_cmp(&P1->x,&P2->x)==0){
        if(Fp_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp_ECD(ANS,P1);
            return;
        }
    }

    EFp_set(&tmp1_EFp,P1);
    EFp_set(&tmp2_EFp,P2);

    Fp_sub(&tmp1_Fp,&tmp2_EFp.x,&tmp1_EFp.x);
    Fp_inv(&tmp1_Fp,&tmp1_Fp);
    Fp_sub(&tmp2_Fp,&tmp2_EFp.y,&tmp1_EFp.y);
    Fp_mul(&tmp3_Fp,&tmp1_Fp,&tmp2_Fp);
    Fp_mul(&tmp1_Fp,&tmp3_Fp,&tmp3_Fp);


    Fp_sub(&tmp2_Fp,&tmp1_Fp,&tmp1_EFp.x);
    Fp_sub(&ANS->x,&tmp2_Fp,&tmp2_EFp.x);

    Fp_sub(&tmp1_Fp,&tmp1_EFp.x,&ANS->x);
    Fp_mul(&tmp2_Fp,&tmp3_Fp,&tmp1_Fp);
    Fp_sub(&ANS->y,&tmp2_Fp,&tmp1_EFp.y);
}
void EFp_SCM(EFp *ANS,EFp *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp_set(ANS,P);
        return;
    }

    EFp Tmp_P,Next_P;
    EFp_init(&Tmp_P);
    EFp_set(&Tmp_P,P);
    EFp_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);

    EFp_set(&Next_P,&Tmp_P);
    for(i=1;i<length; i++){
        EFp_ECD(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp_ECA(&Next_P,&Next_P,&Tmp_P);
        }
    }
    EFp_set(ANS,&Next_P);
}
void EFp_printf(EFp *P,char *str){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp_printf(&P->x,"");
        printf(",");
        Fp_printf(&P->y,"");
        printf(")");
    }else{
        printf("0");
    }
}
int main(){
  EFp a,b,c;
  mpz_t t;
  EFp_init(&a);
  EFp_init(&b);
  EFp_init(&c);
  mpz_init(t);
  getchar();
  printf("fff\n");
  EFp_rational_point(&a);
  getchar();
  printf("fff\n");
  EFp_ECA(&b,&a,&a);
  EFp_ECD(&c,&a);

  if(EFp_cmp(&b,&c)==1){
    printf("ok!\n");
  }

	return 0;
}
