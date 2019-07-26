//测试程序
#include<iostream>
#include"digitalEnvelop.h"
using namespace std;

//用于测试的压缩，解压函数
//在此不做任何压缩，解压
unsigned int compress(mpz_t result,mpz_t mess,unsigned int bits){
	mpz_set(result,mess);
	return bits;
}
unsigned int decompress(mpz_t result,mpz_t mess,unsigned int bits){
	mpz_set(result,mess);
	return bits;
}

int main(){
	gmp_randstate_t state;
	gmp_randinit_default(state);
	unsigned long int seed;
	syscall(SYS_getrandom,&seed,sizeof(unsigned long),1);
	gmp_randseed_ui(state,seed);

	mpz_t message;
	mpz_init(message);
	mpz_urandomb(message,state,9999);

	cout<<"待加密的数据(长度："<<mpz_sizeinbase(message,2)<<" bits"<<endl;
	mpz_out_str(NULL,16,message);
	cout<<endl;

	Curve acurve;
	User A,B;
	generate_key(A,acurve);
	generate_key(B,acurve);
	mpz_t result;
	mpz_init(result);
	letter(result,message,mpz_sizeinbase(message,2),B.pA,A.dA,acurve,compress);

	cout<<"加密后的数据："<<endl;
	mpz_out_str(NULL,16,result);
	cout<<endl;

	unsigned int lenn=letter_de(result,result,A.pA,B.dA,acurve,decompress);
	cout<<"解密后数据,长度："<<lenn<<endl;
	mpz_out_str(NULL,16,result);
	cout<<endl;

	mpz_clear(message);
	mpz_clear(result);
}

