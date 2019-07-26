#ifndef _SM3_H
#define _SM3_H
#include<cstddef>
typedef unsigned int uint;
typedef unsigned long ulong;

#define shiftLeft(x,j) ((x<<j)|(x>>(32-j)))
#define P0(x) (x^shiftLeft(x,9)^shiftLeft(x,17))
#define P1(x) (x^shiftLeft(x,15)^shiftLeft(x,23))
#define FF(x,y,z,j) ((j<16)?(x^y^z):((x&y) | (x&z) | (y&z)))
#define GG(x,y,z,j) ((j<16)?(x^y^z):((x&y) | ((~x)&z)))

const uint T[64] = {
	0x79cc4519, 0xf3988a32, 0xe7311465, 0xce6228cb, 0x9cc45197, 0x3988a32f, 0x7311465e, 0xe6228cbc,
	0xcc451979, 0x988a32f3, 0x311465e7, 0x6228cbce, 0xc451979c, 0x88a32f39, 0x11465e73, 0x228cbce6,
	0x9d8a7a87, 0x3b14f50f, 0x7629ea1e, 0xec53d43c, 0xd8a7a879, 0xb14f50f3, 0x629ea1e7, 0xc53d43ce,
	0x8a7a879d, 0x14f50f3b, 0x29ea1e76, 0x53d43cec, 0xa7a879d8, 0x4f50f3b1, 0x9ea1e762, 0x3d43cec5,
	0x7a879d8a, 0xf50f3b14, 0xea1e7629, 0xd43cec53, 0xa879d8a7, 0x50f3b14f, 0xa1e7629e, 0x43cec53d,
	0x879d8a7a, 0xf3b14f5,  0x1e7629ea, 0x3cec53d4, 0x79d8a7a8, 0xf3b14f50, 0xe7629ea1, 0xcec53d43,
	0x9d8a7a87, 0x3b14f50f, 0x7629ea1e, 0xec53d43c, 0xd8a7a879, 0xb14f50f3, 0x629ea1e7, 0xc53d43ce,
	0x8a7a879d, 0x14f50f3b, 0x29ea1e76, 0x53d43cec, 0xa7a879d8, 0x4f50f3b1, 0x9ea1e762, 0x3d43cec5
};

struct Message{
	uint* arr;//存储消息
	ulong countBits;//比特数
	ulong size;//arr数组长度
	uint W[68];//消息扩展以及杂凑结果存于W【0-7】
	Message(uint* anarr,ulong sz){
		if(anarr!=NULL){
			countBits = sz;
			size = (31 + countBits) / 32;
			arr = new uint[size];
			ulong i;
			for (i = 0; i < size; i++)
				arr[i] = anarr[i];
				

		}
		else{
			arr=NULL;
			countBits=0;
			size=0;
		}
	}
	~Message() {
		delete[] arr;
	}
};

void CF(uint word[8], Message& mess,ulong begin) {
	//消息扩展
	for (int i = 0; i < 16; i++)
		mess.W[i] = mess.arr[i+begin*16];
	for (int i = 16; i < 68; i++) {
		uint a = mess.W[i - 16] ^ mess.W[i - 9] ^ shiftLeft(mess.W[i - 3], 15);
		mess.W[i] = P1(a) ^ shiftLeft(mess.W[i - 13], 7) ^ mess.W[i - 6];
	}
	
	//压缩过程
	uint wordtemp[8];
	for (int i = 0; i < 8; i++)
		wordtemp[i] = word[i];
	for (int j = 0; j < 64; j++) {
		uint ss1, ss2, tt1, tt2;
		ss2 = shiftLeft(word[0], 12);
		ss1 = ss2 +T[j]  + word[4];
		ss1=shiftLeft(ss1, 7);
		ss2 ^= ss1;
		tt1 = FF(word[0], word[1], word[2], j) + word[3] + ss2 + (mess.W[j]^mess.W[j+4]);
		tt2 = GG(word[4], word[5], word[6], j) + word[7] + ss1 + mess.W[j];
		word[3] = word[2];
		word[2] = shiftLeft(word[1], 9);
		word[1] = word[0];
		word[0] = tt1;
		word[7] = word[6];
		word[6] = shiftLeft(word[5], 19);
		word[5] = word[4];
		word[4] = P0(tt2);
	}
	for (int i = 0; i < 8; i++)
		word[i] ^= wordtemp[i];
}
void iteration(Message& mess) {
	//消息填充
	ulong newsz = 512 - mess.countBits % 512;
	if (newsz < 65)newsz += 512;
	newsz = (newsz + mess.countBits) / 512;
	uint* newarr = new uint[newsz*16];
	ulong i = 0;
	for (i = 0; i < mess.size; i++)
		newarr[i] = mess.arr[i];
	
	if (mess.countBits % 32)//添加1
		newarr[i - 1] |= (1 << (31 - mess.countBits % 32));
	else newarr[i++] = 1 << 31;

	for (; i < newsz * 16; i++)//k个0
		newarr[i] = 0;

	newarr[newsz * 16 - 1] = (uint)mess.countBits;//最后64位
	newarr[newsz * 16 - 2] = uint(mess.countBits / 0x100000000);
	
	uint* temp = mess.arr;
	mess.arr = newarr;
	delete[] temp;
	mess.countBits = newsz * 512;
	mess.size = newsz*16;

	//迭代过程
	uint word[8] = { 0x7380166f,0x4914b2b9,0x172442d7,0xda8a0600,0xa96f30bc,0x163138aa,0xe38dee4d,0xb0fb0e4e };
	ulong sz = mess.countBits / 512;
	i=0;
	while (sz--)
		CF(word, mess,i++);
	for(int i=0;i<8;i++)
		mess.W[i]=word[i];
}
uint* sm3(uint* arr,ulong sz){
	Message* mess=new Message(arr,sz);
	iteration(*mess);
	uint* k=new uint[8];
	for(int i=0;i<8;i++)
		k[i]=mess->W[i];
	delete mess;
	return k;
}

#endif
