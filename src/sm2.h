#ifndef _SM2_H
#define _SM2_H
#include<gmp.h>
#include<unistd.h>
#include<sys/syscall.h>

extern unsigned int* sm3(unsigned int*,unsigned long);

struct point;
struct User;
struct Signature;
struct Curve;

//signature.cpp
void generate_key(User& ,Curve& );
void sign_generate(Signature& ,mpz_t ,mpz_t ,Curve& );
bool sign_verify(Signature& ,mpz_t ,point& ,Curve& );

//convert.cpp
unsigned int gmpToInt(unsigned int*& ,mpz_t ,unsigned int );
void intToGmp(mpz_t ,unsigned int* ,unsigned long );
void pointToGmp(mpz_t ,point& ,mpz_t );
void gmpToPoint(point& ,mpz_t ,Curve& );


//sm2_encryption.cpp
void splice(unsigned int*& ,mpz_t ,unsigned int ,mpz_t ,unsigned int );
void KDF(mpz_t ,unsigned int* ,unsigned int ,unsigned int );
void sm2_encryption(mpz_t ,mpz_t ,point& ,unsigned int ,Curve& );
void sm2_decryption(mpz_t ,mpz_t ,mpz_t ,Curve& );


struct point{
	mpz_t x;
	mpz_t y;
	point(){
		mpz_init(x);
		mpz_init(y);
	}
	~point(){
		mpz_clear(x);
		mpz_clear(y);
	}
};
struct User{
	mpz_t dA;
	point pA;
	User(){
		mpz_init(dA);
	}
	~User(){
		mpz_clear(dA);
	}
};
struct Signature{
	mpz_t r;
	mpz_t s;
	Signature(){
		mpz_init(r);
		mpz_init(s);
	}
	~Signature(){
		mpz_clear(r);
		mpz_clear(s);
	}
};
struct Curve{//Fp-256
	mpz_t p_prime;
	mpz_t a_coeffi;
	mpz_t b_coeffi;
	point G;
	mpz_t n_ranks;

	mpz_t lambda;//计算过程的临时变量
	mpz_t temp;
	mpz_t tempx;
	mpz_t u_sqrt;//sqrt
	Curve(){
		mpz_init_set_str(p_prime,"8542d69e4c044f18e8b92435bf6ff7de457283915c45517d722edb8b08f1dfc3",16);
		mpz_init_set_str(a_coeffi,"787968b4fa32c3fd2417842e73bbfeff2f3c848b6831d7e0ec65228b3937e498",16);
		mpz_init_set_str(b_coeffi,"63e4c6d3b23b0c849cf84241484bfe48f61d59a5b16ba06e6e12d1da27c5249a",16);
		mpz_init_set_str(G.x,"421debd61b62eab6746434ebc3cc315e32220b3badd50bdc4c4e6c147fedd43d",16);
		mpz_init_set_str(G.y,"0680512bcbb42c07d47349d2153b70c4e5d7fdfcbfa36ea1a85841b9e46e09a2",16);
		mpz_init_set_str(n_ranks,"8542d69e4c044f18e8b92435bf6ff7dd297720630485628d5ae74ee7c32e79b7",16);
	
		mpz_init(lambda);
		mpz_init(temp);
		mpz_init(tempx);
		mpz_init_set_str(u_sqrt,"2150b5a7930113c63a2e490d6fdbfdf7915ca0e45711545f5c8bb6e2c23c77f1",16);//+1
	}
	~Curve(){
		mpz_clear(p_prime);
		mpz_clear(a_coeffi);
		mpz_clear(b_coeffi);
		mpz_clear(G.x);
		mpz_clear(G.y);
		mpz_clear(n_ranks);

		mpz_clear(lambda);
		mpz_clear(temp);
		mpz_clear(tempx);
		mpz_clear(u_sqrt);
	}
	void add(point& result,point& p1,point& p2){
		int x1=mpz_cmp_ui(p1.x,0);
		int x2=mpz_cmp_ui(p2.x,0);
		int t1=mpz_cmp(p1.x,p2.x);
		int t2=mpz_cmp(p1.y,p2.y);
		int t3=mpz_cmpabs(p1.y,p2.y);

		if((!x1&&!x2)||(!t1&&t2&&!t3)){//都为无穷点或互逆，结果无穷点
			mpz_set_ui(result.x,0);
			mpz_set_ui(result.y,0);
			return;
		}

		//二者其一为无穷点
		if(!x1){
			mpz_set(result.x,p2.x);
			mpz_set(result.y,p2.y);
			return;
		}
		if(!x2){
			mpz_set(result.x,p1.x);
			mpz_set(result.y,p1.y);
			return;
		}

		if(!t1&&!t2){ //两点相同
			mpz_mul(temp,p1.x,p1.x);// temp = p1.x - p2.x
			mpz_mul_ui(temp,temp,3);//  temp *= 3
			mpz_add(temp,temp,a_coeffi);  
			mpz_mod(temp,temp,p_prime);
			mpz_mul_ui(lambda,p1.y,2);
			mpz_mod(lambda,lambda,p_prime);
			mpz_invert(lambda,lambda,p_prime);
			mpz_mul(lambda,temp,lambda);
		}
		else{ //两个非互逆不同点
			mpz_sub(temp,p2.y,p1.y);
			mpz_sub(lambda,p2.x,p1.x);
			mpz_mod(temp,temp,p_prime);
			mpz_mod(lambda,lambda,p_prime);
			mpz_invert(lambda,lambda,p_prime);
			mpz_mul(lambda,temp,lambda);
		}
		mpz_mod(lambda,lambda,p_prime);
		mpz_mul(temp,lambda,lambda);
		mpz_mod(temp,temp,p_prime);
		mpz_sub(temp,temp,p1.x);
		mpz_sub(temp,temp,p2.x);

		mpz_set(tempx,p1.x);
		mpz_mod(result.x,temp,p_prime);

		mpz_sub(temp,tempx,result.x);
		mpz_mul(temp,temp,lambda);
		mpz_sub(temp,temp,p1.y);
		mpz_mod(result.y,temp,p_prime);
	}
	void mul(point& result,point& p,mpz_t k){
		if(!mpz_cmp_ui(k,2))return add(result,p,p);
		mpz_xor(result.x,result.x,result.x);
		mpz_xor(result.y,result.y,result.y);
		unsigned long index=mpz_sizeinbase(k,2);
		while(index--){
			add(result,result,result);
			if(mpz_tstbit(k,index))
				add(result,result,p);
		}
	}
	bool compute(mpz_t result,mpz_t x){
		mpz_pow_ui(temp,x,3);
		mpz_mul(tempx,x,a_coeffi);
		mpz_add(temp,temp,tempx);
		mpz_add(temp,temp,b_coeffi);
		mpz_mod(temp,temp,p_prime);
		if(mpz_cmp_ui(temp,0)){
			mpz_set_ui(result,0);
			return true;
		}
		mpz_powm(result,temp,u_sqrt,p_prime);
		mpz_powm_ui(result,result,2,p_prime);
		if(!mpz_cmp(result,temp))return true;
		else return false;
	}
};

#endif //_SM2_H
