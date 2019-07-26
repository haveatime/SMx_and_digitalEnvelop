#ifndef _DIGITALENVELOP_H
#define _DIGITALENVELOP_H

#include"sm3.h"
#include"sm2.h"
#include"sm4.h"

//数字信封生成，带认证功能
//param:	
//		result  :生成的数字信封
//		message :明文
//		bits    :发送数据message的bit数
//		pkey_B  :接收方公钥
//		skey_A  :发送方私钥
//		acurve  :所使用的椭圆曲线
//      compress:压缩函数，用以压缩明文
//				 函数参数(压缩结果，待压缩数据，待压缩数据bits）
//				 函数返回：压缩结果bits
void letter(mpz_t result,mpz_t message,unsigned int bits,point& pkey_B,mpz_t skey_A,Curve& acurve,unsigned int(*compress)(mpz_t,mpz_t,unsigned int));

//数字信封解开
//param：
//		result    :明文	
//		message   :收到的数字信封
//		pkey_A    :发送方公钥
//		skey_B    :接收方私钥
//		acurve    :所使用的椭圆曲线
//		decompress:解压函数，函数返回解压结果的bits
//				   函数参数(解压结果，待解压数据，待解压数据bits），
//return: 明文的bit数，若数字签名认证失败，则为0
unsigned int letter_de(mpz_t result,mpz_t message,point& pkey_A,mpz_t skey_B,Curve& acurve,unsigned int (*decompress)(mpz_t,mpz_t,unsigned int));



void letter(mpz_t result,mpz_t message,unsigned int bits,point& pkey_B,mpz_t skey_A,Curve& acurve,unsigned int (*compress)(mpz_t,mpz_t,unsigned int)){	
	srand((int)time(0));
	unsigned int key[4];//随机生成sm4密钥
	key[0]=rand();
	key[1]=rand();
	key[2]=rand();
	key[3]=rand();
	
	unsigned int* mess;
	unsigned int mess_size=gmpToInt(mess,message,bits);
	unsigned int* pointer=sm3(mess,mess_size);
	mpz_t temp;
	Signature sign;
	mpz_init(temp);
	intToGmp(temp,pointer,256);
	sign_generate(sign,temp,skey_A,acurve);//生成了消息的数字签名sign
	delete pointer;

	unsigned int rlen=mpz_sizeinbase(sign.r,2),slen=mpz_sizeinbase(sign.s,2);
	splice(mess,sign.r,rlen,sign.s,slen);// r||s
	intToGmp(temp,mess,rlen+slen);
	delete mess;

	mpz_t cpm;
	mpz_init(cpm);
	unsigned int len=compress(cpm,message,bits);//压缩明文
	splice(mess,cpm,len,temp,rlen+slen);// cpm||r||s
	mess_size=len+rlen+slen;
	mess_size=sm4_encryption(pointer,mess,mess_size,key);
	intToGmp(temp,pointer,mess_size);//对称加密：消息使用sm4加密

	mpz_mul_2exp(temp,temp,32);//在其后添加32bits
	mpz_add_ui(temp,temp,rlen<<16);//用于记录数字签名r，s各自长度
	mpz_add_ui(temp,temp,slen);

	delete pointer;
	delete mess;

	intToGmp(result,&key[0],128);	
	sm2_encryption(result,result,pkey_B,128,acurve);
	mpz_mul_2exp(result,result,mess_size+32);
	mpz_xor(result,result,temp);// result||temp
	mpz_clear(temp);
}


unsigned int letter_de(mpz_t result,mpz_t message,point& pkey_A,mpz_t skey_B,Curve& acurve,unsigned int (*decompress)(mpz_t,mpz_t,unsigned int)){
	Signature sign;
	mpz_t C1,C2;
	mpz_init(C1);
	mpz_init(C2);
	
	unsigned int bits;
	unsigned len=(7+mpz_sizeinbase(acurve.p_prime,2))>>3;
	len=mpz_sizeinbase(message,2)-387-len*2*8;
	mpz_tdiv_q_2exp(C1,message,len);
	mpz_tdiv_r_2exp(C2,message,len);

	//提取出签名r，s的长度
	mpz_init_set_ui(sign.r,0xffff);
	mpz_and(sign.s,C2,sign.r);
	mpz_mul_2exp(sign.r,sign.r,16);
	mpz_and(sign.r,C2,sign.r);
	mpz_tdiv_q_2exp(sign.r,sign.r,16);
	mpz_tdiv_q_2exp(C2,C2,32);

	len-=32;
	sm2_decryption(C1,C1,skey_B,acurve);//sm2解密获得密钥
	unsigned int *key,*mess,*c2;
	len=gmpToInt(c2,C2,len);
	gmpToInt(key,C1,128);
	len=sm4_decryption(mess,c2,len,key);//sm4加密
	bits=len;
	intToGmp(result,mess,len);

	//提取签名r，s
	len=mpz_get_ui(sign.s);
	bits-=len;
	mpz_tdiv_r_2exp(sign.s,result,len);
	mpz_tdiv_q_2exp(result,result,len);

	len=mpz_get_ui(sign.r);
	bits-=len;
	mpz_tdiv_r_2exp(sign.r,result,len);
	mpz_tdiv_q_2exp(C2,result,len);
	bits=decompress(result,C2,bits);//解压

	//签名认证
	delete mess;
	unsigned int mess_size=gmpToInt(mess,result,bits);
	unsigned int* hash=sm3(mess,mess_size);
	intToGmp(C1,hash,256);
	if(!sign_verify(sign,C1,pkey_A,acurve))bits=0;//认证失败，返回0

	mpz_clear(C1);
	mpz_clear(C2);
	delete hash;
	delete key;
	return bits;
}

#endif
