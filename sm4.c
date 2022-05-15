#include <stdio.h>
#include <iostream>
#include <pthread.h>
#include <immintrin.h>
using namespace std;
#define u8 unsigned char
#define u32 unsigned long

#define NUM_P 32

// S盒
const u8 Sbox[256] = {
	0xd6,0x90,0xe9,0xfe,0xcc,0xe1,0x3d,0xb7,0x16,0xb6,0x14,0xc2,0x28,0xfb,0x2c,0x05,
	0x2b,0x67,0x9a,0x76,0x2a,0xbe,0x04,0xc3,0xaa,0x44,0x13,0x26,0x49,0x86,0x06,0x99,
	0x9c,0x42,0x50,0xf4,0x91,0xef,0x98,0x7a,0x33,0x54,0x0b,0x43,0xed,0xcf,0xac,0x62,
	0xe4,0xb3,0x1c,0xa9,0xc9,0x08,0xe8,0x95,0x80,0xdf,0x94,0xfa,0x75,0x8f,0x3f,0xa6,
	0x47,0x07,0xa7,0xfc,0xf3,0x73,0x17,0xba,0x83,0x59,0x3c,0x19,0xe6,0x85,0x4f,0xa8,
	0x68,0x6b,0x81,0xb2,0x71,0x64,0xda,0x8b,0xf8,0xeb,0x0f,0x4b,0x70,0x56,0x9d,0x35,
	0x1e,0x24,0x0e,0x5e,0x63,0x58,0xd1,0xa2,0x25,0x22,0x7c,0x3b,0x01,0x21,0x78,0x87,
	0xd4,0x00,0x46,0x57,0x9f,0xd3,0x27,0x52,0x4c,0x36,0x02,0xe7,0xa0,0xc4,0xc8,0x9e,
	0xea,0xbf,0x8a,0xd2,0x40,0xc7,0x38,0xb5,0xa3,0xf7,0xf2,0xce,0xf9,0x61,0x15,0xa1,
	0xe0,0xae,0x5d,0xa4,0x9b,0x34,0x1a,0x55,0xad,0x93,0x32,0x30,0xf5,0x8c,0xb1,0xe3,
	0x1d,0xf6,0xe2,0x2e,0x82,0x66,0xca,0x60,0xc0,0x29,0x23,0xab,0x0d,0x53,0x4e,0x6f,
	0xd5,0xdb,0x37,0x45,0xde,0xfd,0x8e,0x2f,0x03,0xff,0x6a,0x72,0x6d,0x6c,0x5b,0x51,
	0x8d,0x1b,0xaf,0x92,0xbb,0xdd,0xbc,0x7f,0x11,0xd9,0x5c,0x41,0x1f,0x10,0x5a,0xd8,
	0x0a,0xc1,0x31,0x88,0xa5,0xcd,0x7b,0xbd,0x2d,0x74,0xd0,0x12,0xb8,0xe5,0xb4,0xb0,
	0x89,0x69,0x97,0x4a,0x0c,0x96,0x77,0x7e,0x65,0xb9,0xf1,0x09,0xc5,0x6e,0xc6,0x84,
	0x18,0xf0,0x7d,0xec,0x3a,0xdc,0x4d,0x20,0x79,0xee,0x5f,0x3e,0xd7,0xcb,0x39,0x48
};

// 密钥扩展算法的常数FK 
const u32 FK[4] = {
	0xa3b1bac6, 0x56aa3350, 0x677d9197, 0xb27022dc
};

// 密钥扩展算法的固定参数CK 
const u32 CK[32] = {
	0x00070e15, 0x1c232a31, 0x383f464d, 0x545b6269,
	0x70777e85, 0x8c939aa1, 0xa8afb6bd, 0xc4cbd2d9,
	0xe0e7eef5, 0xfc030a11, 0x181f262d, 0x343b4249,
	0x50575e65, 0x6c737a81, 0x888f969d, 0xa4abb2b9,
	0xc0c7ced5, 0xdce3eaf1, 0xf8ff060d, 0x141b2229,
	0x30373e45, 0x4c535a61, 0x686f767d, 0x848b9299,
	0xa0a7aeb5, 0xbcc3cad1, 0xd8dfe6ed, 0xf4fb0209,
	0x10171e25, 0x2c333a41, 0x484f565d, 0x646b7279
};
u32 *X;  // 明文 
u32 *Y;//密文
u32 MK[4] = { 0x01234567 ,0x89abcdef, 0xfedcba98, 0x76543210 }; // 密钥 
u32 RK[32]; // 轮密钥  
u32 K[4]; // 中间数据 

u32 find_sbox(u32 b); // 查S盒的函数B 
u32 loop_left(u32 a, short length); // 循环左移函数 
u32 L1(u32 a); // 线性变换L
u32 L2(u32 a); // 线性变换L'
u32 T(u32 a, short mode); // 合成变换T
void extend1(u32 MK[], u32 K[]); // 密钥扩展算法第一步
void extend2(u32 RK[], u32 K[]); // 密钥扩展算法第二步
void set_RK(u32 MK[], u32 K[], u32 RK[]); // 轮密钥获取算法
void round(u32 X[], u32 RK[]); // 迭代算法
void reverse(u32 X[], u32 Y[]); // 反转函数 
void encryptSM4(u32 X[], u32 RK[], u32 Y[]); // 加密算法
void decryptSM4(u32 X[], u32 RK[], u32 Y[]); // 解密算法



u32 find_sbox(u32 b) {
	u8 a[4];
	short i;
	a[0] = (u8)(b >> 24);
	a[1] = (u8)(b >> 16);
	a[2] = (u8)(b >> 8);
	a[3] = (u8)b;
	b = Sbox[a[0]]<<24 | Sbox[a[1]] <<16 | Sbox[a[2]] <<8|  Sbox[a[3]];
	return b;
}



u32 loop_left(u32 a, short length) {
	a = a << length | a >> (32 - length);
	return a;
}


u32 L1(u32 a) {
	return a ^ loop_left(a, 2) ^ loop_left(a, 10) ^ loop_left(a, 18) ^ loop_left(a, 24);
}


u32 L2(u32 a) {
	return a ^ loop_left(a, 13) ^ loop_left(a, 23);
}


u32 T(u32 a, short mode) {
	return mode == 1 ? L1(find_sbox(a)) : L2(find_sbox(a));
}
__m256i T_simd(__m256i x) {
	u32* a = (u32*)&x;
	for (int i = 0; i < 8; i++)
	{
		a[i] = L1(find_sbox(a[i]));
	}
	return *(__m256i*)a;
}

void extend1(u32 MK[], u32 K[]) {
	int i;
	for (i = 0; i < 4; i++) {
		K[i] = MK[i] ^ FK[i];
	}
}


void extend2(u32 RK[], u32 K[]) {
	short i;
	for (i = 0; i < 32; i++) {
		K[(i + 4) % 4] = K[i % 4] ^ T(K[(i + 1) % 4] ^ K[(i + 2) % 4] ^ K[(i + 3) % 4] ^ CK[i], 2);
		RK[i] = K[(i + 4) % 4];
	}
}


void set_RK(u32 MK[], u32 K[], u32 RK[]) {
	extend1(MK, K);
	extend2(RK, K);
}


void round(u32 X[], u32 RK[]) {
	short i;
	for (i = 0; i < 32; i++) {
		X[(i + 4) % 4] = X[i % 4] ^ T(X[(i + 1) % 4] ^ X[(i + 2) % 4] ^ X[(i + 3) % 4] ^ RK[i], 1);
	}
}
void round_simd(u32 X[], u32 RK[]) {
	short i;
	__m256i x0, x1, x2, x3, rk;
	for (i = 0; i < 32; i++) {
		x0 = _mm256_set_epi32(X[i % 4], X[i % 4 + 4], X[i % 4 + 8], X[i % 4 + 12], X[i % 4 + 16], X[i % 4 + 20], X[i % 4 + 24], X[i % 4 + 28]);
		x1 = _mm256_set_epi32(X[(i + 1) % 4], X[(i + 1) % 4 + 4], X[(i + 1) % 4 + 8], X[(i + 1) % 4 + 12], X[(i + 1) % 4 + 16], X[(i + 1) % 4 + 20], X[(i + 1) % 4 + 24], X[(i + 1) % 4 + 28]);
		x2 = _mm256_set_epi32(X[(i + 2) % 4], X[(i + 2) % 4 + 4], X[(i + 2) % 4 + 8], X[(i + 2) % 4 + 12], X[(i + 2) % 4 + 16], X[(i + 2) % 4 + 20], X[(i + 2) % 4 + 24], X[(i + 2) % 4 + 28]);
		x3 = _mm256_set_epi32(X[(i + 3) % 4], X[(i + 3) % 4 + 4], X[(i + 3) % 4 + 8], X[(i + 3) % 4 + 12], X[(i + 3) % 4 + 16], X[(i + 3) % 4 + 20], X[(i + 3) % 4 + 24], X[(i + 3) % 4 + 28]);
		rk = _mm256_set1_epi32(RK[i]);
		x1 = _mm256_xor_si256(x1, x2);
		x1 = _mm256_xor_si256(x1, x3);
		x1 = _mm256_xor_si256(x1, rk);
		x1 = T_simd(x1);
		x0 = _mm256_xor_si256(x0, x1);
		//X[(i + 4) % 4] = X[i % 4] ^ T(X[(i + 1) % 4] ^ X[(i + 2) % 4] ^ X[(i + 3) % 4] ^ RK[i], 1);
		u32* a = (u32*)&x0;
		for (int j = 0; j < 8; j++)
		{
			X[(i + 4) % 4 + 4 * j] = a[j];
		}
	}
}


void reverse(u32 X[], u32 Y[]) {
	short i;
	for (i = 0; i < 4; i++) {
		Y[i] = X[4 - 1 - i];
	}
}
void reverse_simd(u32 X[], u32 Y[]) {
	short i;
	for (i = 0; i < 4; i++) {
		Y[i] = X[4 - 1 - i];
		Y[i+4] = X[4 - 1 - i+4];
		Y[i+8] = X[4 - 1 - i+8];
		Y[i+12] = X[4 - 1 - i+12];
		Y[i+16] = X[4 - 1 - i+16];
		Y[i+20] = X[4 - 1 - i+20];
		Y[i+24] = X[4 - 1 - i+24];
		Y[i+28] = X[4 - 1 - i+28];
	}
}


void encryptSM4(u32 X[], u32 RK[], u32 Y[]) {
	round(X, RK);
	reverse(X, Y);
}

void encryptSM4_simd(u32 X[], u32 RK[], u32 Y[]) {

	round_simd(X, RK);
	reverse_simd(X, Y);
}


void decryptSM4(u32 X[], u32 RK[], u32 Y[]) {
	short i;
	u32 reverseRK[32];
	for (i = 0; i < 32; i++) {
		reverseRK[i] = RK[32 - 1 - i];
	}
	round(X, reverseRK);
	reverse(X, Y);
}

