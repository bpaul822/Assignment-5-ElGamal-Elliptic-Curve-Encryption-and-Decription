#ifndef _BIGNUM_H
#define _BIGNUM_H
#include<string>
using namespace std;
const int Size2048 = 309;

struct BigNum
{
    int Num[309] = {};
    bool negative = false;
};

struct DivResult
{
    BigNum Result;
    BigNum Remainder;
};

struct ArrayOfArray
{
    BigNum Result;
    BigNum Count;
};

DivResult DivSmall(BigNum first, BigNum second);
BigNum AddFront(BigNum input, int val);
DivResult DivLarge(BigNum first, BigNum second);
bool IsPrime(BigNum input);
BigNum Inverse(BigNum input, BigNum mod);
BigNum Sub(BigNum firstOriginal, BigNum second);
BigNum gcd(BigNum a, BigNum b);
int Compare(BigNum a, BigNum b);
BigNum StringToArray(string input);
BigNum CopyOf(BigNum input);
bool Equalone(BigNum input);
string value_number(BigNum input);
BigNum PwrMod(BigNum firstOriginal, BigNum secondOriginal, BigNum Mod);
bool EqualZero(BigNum input);
BigNum Add(BigNum first, BigNum second);
BigNum Mul(BigNum first, BigNum second);

#endif 
