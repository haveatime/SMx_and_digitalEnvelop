#include"sm2.h"

extern unsigned int* sm3(unsigned int*,unsigned long);

//生成公钥私钥
void generate_key(User& auser,Curve& acurve){
	gmp_randstate_t state;
	gmp_randinit_default(state);
	unsigned int seed;
	syscall(SYS_getrandom,&seed,sizeof(unsigned long int),1);
	gmp_randseed_ui(state,seed);
	do{
		mpz_urandomm(auser.dA,state,acurve.n_ranks);
	}while(!mpz_cmp_ui(auser.dA,0));
	gmp_randclear(state);
	
	acurve.mul(auser.pA,acurve.G,auser.dA);
}

//数字签名生成
void sign_generate(Signature& result,mpz_t message,mpz_t dA,Curve& acurve){
	gmp_randstate_t state;
	gmp_randinit_default(state);
	unsigned int seed;
	syscall(SYS_getrandom,&seed,sizeof(unsigned long int),1);
	gmp_randseed_ui(state,seed);

	mpz_t k,temp;
	point apoint;
	unsigned int sz;
	unsigned int *mess,*hashResult;
	mpz_init(k);
	mpz_init(temp);
	while(true){
		do{
			mpz_urandomm(k,state,acurve.n_ranks);
		}while(!mpz_cmp_ui(k,0));//k属于[1,n-1]
		acurve.mul(apoint,acurve.G,k);
		sz=gmpToInt(mess,message,mpz_sizeinbase(message,16)*4);
		hashResult=sm3(mess,sz);	
		intToGmp(result.r,hashResult,256);
		mpz_add(result.r,result.r,apoint.x);
		mpz_mod(result.r,result.r,acurve.n_ranks);
		
		//r==0 或 r+k==n 则返回重新生成k
		mpz_add(temp,result.r,k);
		if(!mpz_cmp_ui(result.r,0)||!mpz_cmp(acurve.n_ranks,temp))continue;
		
		mpz_mul(temp,result.r,dA);
		mpz_sub(temp,k,temp);
		mpz_add_ui(result.s,dA,1);
		mpz_invert(result.s,result.s,acurve.n_ranks);
		mpz_mul(result.s,result.s,temp);
		mpz_mod(result.s,result.s,acurve.n_ranks);
		
		if(!mpz_cmp_ui(result.s,0))continue;//s==0 重新生成k
		break;
	}
	delete mess;
	delete hashResult;
	gmp_randclear(state);
	mpz_clear(temp);
	mpz_clear(k);
}

//签名验证
bool sign_verify(Signature& result,mpz_t message,point& pA,Curve& acurve){
	unsigned int* mess;
	unsigned int sz=gmpToInt(mess,message,mpz_sizeinbase(message,16)*4);
	unsigned int* hashResult=sm3(mess,sz);

	mpz_t e,t;
	mpz_init(e);
	mpz_init(t);
	intToGmp(e,hashResult,256);
	
	mpz_add(t,result.r,result.s);
	mpz_mod(t,t,acurve.n_ranks);
	point pointa,pointb;
	acurve.mul(pointa,acurve.G,result.s);
	acurve.mul(pointb,pA,t);
	acurve.add(pointa,pointa,pointb);
	mpz_add(e,e,pointa.x);
	mpz_mod(e,e,acurve.n_ranks);
	bool flag;
	if(!mpz_cmp_ui(t,0)||mpz_cmp(result.r,e))flag=false;//t==0 或 r！=e 验证不通过
	else flag=true;
	delete mess;
	delete hashResult;
	mpz_clear(t);
	mpz_clear(e);
	return flag;
}
