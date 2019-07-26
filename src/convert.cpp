#include"sm2.h"
#include<cstdlib>

const unsigned int base[32]={ 0x1,      0x2,      0x4,      0x8,      0x10,      0x20,      0x40,      0x80,
					  0x100,    0x200,    0x400,    0x800,    0x1000,    0x2000,    0x4000,    0x8000,
					  0x10000,  0x20000,  0x40000,  0x80000,  0x100000,  0x200000,  0x400000,  0x800000,
					  0x1000000,0x2000000,0x4000000,0x8000000,0x10000000,0x20000000,0x40000000,0x80000000};

// mpz_t 类型数 转换为 int数组
// len为mpz_t数的比特长度 
unsigned int gmpToInt(unsigned int*& result,mpz_t mess,unsigned int len){
	unsigned long length=mpz_sizeinbase(mess,2);
	if(len>length)length=len;
	unsigned long temp=length;
	result=new unsigned int[(31+length)>>5];
	int bit=31;
	int index=0;
	result[0]=0;
	while(length--){
		if(bit==-1){
			index++;
			result[index]=0;
			bit=31;
		}
		if(mpz_tstbit(mess,length)) result[index]|=base[bit];
		bit--;
	}
	return temp;
}

// int数组 转为 mpz_t 数，bits为数组的比特长度
void intToGmp(mpz_t result,unsigned int* arr,unsigned long bits){
	unsigned long index=bits/32-1;
	if(bits%32)index++;
	unsigned int k=bits%32;
	unsigned long i;
	if(k!=0){
		i=k;
		k=arr[index]>>(32-k);
		mpz_set_ui(result,k);
	}
	else{
		i=0;
		index++;
		mpz_set_ui(result,0);
	}
	while(index--){
		for(int j=0;j<32;j++)
			if(arr[index]>>j&1)mpz_setbit(result,i+j);
		i+=32;
	}
}

//椭圆曲线上 点到比特串 的转换
void pointToGmp(mpz_t result,point& p,mpz_t Fq){
	mpz_t X1,Y1;
	mpz_init_set(X1,p.x);
	mpz_init_set(Y1,p.y);
	unsigned long t=mpz_sizeinbase(Fq,2);
	unsigned long l=(t+7)>>3;
	int y=mpz_tstbit(p.y,0);//点的y坐标最右一比特
	mpz_mul_2exp(result,X1,8*l);
	mpz_xor(result,result,Y1);
	mpz_setbit(result,16*l+2);
	mpz_clear(X1);
	mpz_clear(Y1);
}

// 比特串到椭圆曲线上点 的转换
void gmpToPoint(point& result,mpz_t bitstring,Curve& acurve){
	unsigned long t=mpz_sizeinbase(acurve.p_prime,2);
	unsigned long l=(t+7)>>3;
	unsigned length=mpz_sizeinbase(bitstring,2);
	mpz_tdiv_r_2exp(result.x,bitstring,2*l*8);
	mpz_tdiv_r_2exp(result.y,bitstring,l*8);
	mpz_tdiv_q_2exp(result.x,result.x,l*8);
	unsigned int begin=2*l*8,end=length-1,pc=0;
	while(end>=begin){
		pc<<=1;
		if(mpz_tstbit(bitstring,end--))pc++;
	}
	bool y;
	if(pc!=4)exit(1);
	mpz_t yy;
	mpz_init(yy);
	bool f=acurve.compute(yy,result.x);
	if(!f)exit(1);
	mpz_powm_ui(yy,yy,2,acurve.p_prime);
	if(!mpz_cmp(yy,result.y))exit(1);
	mpz_clear(yy);
}
