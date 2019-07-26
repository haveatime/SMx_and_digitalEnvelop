#include"sm2.h"
#include<cstdlib>

extern unsigned int* sm3(unsigned int*,unsigned long);

//拼接函数 result=a||b
void splice(unsigned int*& result,mpz_t a,unsigned int alen,mpz_t b,unsigned int blen){
	mpz_t temp;
	mpz_init(temp);
	mpz_mul_2exp(temp,a,blen);
	mpz_add(temp,temp,b);
	gmpToInt(result,temp,alen+blen);
	mpz_clear(temp);
}

//密钥派生函数
//输出为klen比特的比特串
//bits：temp数组的比特串长度
void KDF(mpz_t K,unsigned int* temp,unsigned int klen,unsigned int bits){
	unsigned int *arr;
	unsigned int l=(klen+255)>>8,size=(31+bits)>>5,ct=1;
	arr=new unsigned int[size+1];
	for(int i=0;i<size;i++)
		arr[i]=temp[i];
	unsigned int a=bits%32,b,w=0;
	b=32-a;
	unsigned int *hash,result[l*8];
	for(int i=1;i<=l;i++,w+=8){
		//  temp||ct
		arr[size]=ct;
		if(a!=0){
			arr[size-1]|=(ct>>a);
			arr[size]=(ct<<b);
		}

		hash=sm3(arr,bits+32);
		for(int j=0;j<8;j++)
			result[w+j]=hash[j];
		arr[size-1]=temp[size-1];
		ct++;
		delete[] hash;
	}
	while(w<((klen+31)>>5))
		result[w++]=result[w-9];
	intToGmp(K,result,klen);
	delete[] arr;
}

//sm2 公钥加密算法
//result：密文(加密结果）
//message，klen：长度为klen的消息; 
//pA：公钥;
//acurve：椭圆曲线
void sm2_encryption(mpz_t result,mpz_t message,point& pA,unsigned int klen,Curve& acurve){
	gmp_randstate_t state;
	gmp_randinit_default(state);
	unsigned long int seed;
	syscall(SYS_getrandom,&seed,sizeof(unsigned long int),1);
	gmp_randseed_ui(state,seed);

	mpz_t k,t;
	mpz_t C1,C2,C3;
	mpz_init(k);
	mpz_init(t);
	mpz_init(C1);
	mpz_init(C2);
	mpz_init(C3);
	point c1,apoint;
	do{
		do{
			mpz_urandomm(k,state,acurve.n_ranks);
		}while(!mpz_cmp_ui(k,0));//k属于[1,n-1]
		acurve.mul(c1,acurve.G,k);//c1=[k]G
		pointToGmp(C1,c1,acurve.p_prime);
		//省略步骤S=[h]Pb
		acurve.mul(apoint,pA,k);// apoint=[k]pA
		unsigned int* xy;
		splice(xy,apoint.x,mpz_sizeinbase(apoint.x,2),apoint.y,mpz_sizeinbase(apoint.y,2));
		unsigned int size=mpz_sizeinbase(apoint.x,2)+mpz_sizeinbase(apoint.y,2);
		KDF(t,xy,klen,size);
		delete[] xy;
	}while(!mpz_cmp_ui(t,0));

	mpz_xor(C2,message,t);

	// temp = x||M||y
	// C3 = Hash(temp)
	unsigned int* temp;
	unsigned int bits=mpz_sizeinbase(apoint.x,2)+mpz_sizeinbase(message,2);
	splice(temp,apoint.x,mpz_sizeinbase(apoint.x,2),message,mpz_sizeinbase(message,2));

	intToGmp(C3,temp,bits);

	delete[] temp;
	bits=mpz_sizeinbase(C3,2)+mpz_sizeinbase(apoint.y,2);
	splice(temp,C3,mpz_sizeinbase(C3,2),apoint.y,mpz_sizeinbase(apoint.y,2));
	unsigned int *c3=sm3(temp,bits);
	intToGmp(C3,c3,256);
	delete[] temp;

	//result = C1||C2||C3
	mpz_mul_2exp(result,C1,klen+256);
	mpz_mul_2exp(C2,C2,256);
	mpz_add(result,result,C2);
	mpz_add(result,result,C3);
	
	delete[] c3;
	mpz_clear(t);
	mpz_clear(k);
	mpz_clear(C1);
	mpz_clear(C2);
	mpz_clear(C3);
	gmp_randclear(state);
}

//sm2 解密算法
//result：密文
void sm2_decryption(mpz_t result,mpz_t message,mpz_t dB,Curve& acurve){
	mpz_t C1,C2,C3;
	mpz_init(C1);
	mpz_init(C2);
	mpz_init(C3);

	//从message中取出C1，C2，C3
	mpz_setbit(C3,0);
	mpz_setbit(C2,0);
	unsigned long len=mpz_sizeinbase(message,2);
	unsigned long l=(7+mpz_sizeinbase(acurve.p_prime,2))>>3;
	unsigned long klen=len-2*8*l-3-256;
	mpz_tdiv_q_2exp(C1,message,len-2*8*l-3);
	mpz_mul_2exp(C2,C2,len-2*8*l-3);
	mpz_mul_2exp(C3,C3,256);
	mpz_tdiv_r(C2,message,C2);
	mpz_mod(C3,C2,C3);
	mpz_tdiv_q_2exp(C2,C2,256);

	//将C1转为椭圆曲线上点
	point c1,apoint;
	gmpToPoint(c1,C1,acurve);

	acurve.mul(apoint,c1,dB);	
	unsigned int* arr;
	len=mpz_sizeinbase(apoint.x,2)+mpz_sizeinbase(apoint.y,2);
	splice(arr,apoint.x,mpz_sizeinbase(apoint.x,2),apoint.y,mpz_sizeinbase(apoint.y,2));
	KDF(C1,arr,klen,len);
	mpz_xor(result,C1,C2);//求明文
	delete[] arr;

	// arr=  x2||M||y2
	len=mpz_sizeinbase(apoint.x,2)+mpz_sizeinbase(result,2);
	splice(arr,apoint.x,mpz_sizeinbase(apoint.x,2),result,mpz_sizeinbase(result,2));
	intToGmp(C1,arr,len);
	delete[] arr;
	splice(arr,C1,len,apoint.y,mpz_sizeinbase(apoint.y,2));
	len+=mpz_sizeinbase(apoint.y,2);
	unsigned int* hash=sm3(arr,len);
	intToGmp(C1,hash,256);

	if(mpz_cmp(C1,C3))exit(1);
	mpz_clear(C1);
	mpz_clear(C2);
	mpz_clear(C3);
	delete[] hash;
	delete[] arr;
}

