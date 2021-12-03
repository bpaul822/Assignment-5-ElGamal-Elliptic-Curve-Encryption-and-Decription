#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include "bignum.h"

using namespace std;

//Function to compute the sum of 2 points on the elliptic curve
BigNum Gx,Gy,a,b,p,Xr,Yr,Nb,k,One,add_big;
string add_str = "1";
void pointsum (BigNum Xp,BigNum Yp ,BigNum Xq,BigNum Yq) {
    One.Num[0] = 1;
    BigNum yp2;
    BigNum lamda;
    BigNum Inv,diffxpxr,diffyqyp;  
    if ((!Compare(Xp,Xq))  & (!Compare(Yp,Yq))){
      //modInverse(2*Yp,p);
      add_str = "2";
      add_big = StringToArray(add_str);
      yp2 = Add(yp2,add_big);
      Inv = Inverse(Mul(yp2,Yp),p);
        while(Inv.negative) {
          Inv  = Add(Inv,p);
        }
      //lamda = (3*Xp*Xp + a)* modInverse(2*Yp,p);
      lamda = Mul((Add(Mul(Mul(StringToArray("3"),Xp),Xp),a)),Inv); 
    }
    else{
      Inv = Sub(Xq,Xp);
      while(Inv.negative){
        Inv  = Add(Inv,p);}
      diffyqyp =  Sub(Yq,Yp);
      while(diffyqyp.negative){
        diffyqyp = Add(diffyqyp,p);}
      lamda = Mul(diffyqyp,(Inverse(Inv,p))); 
      //lamda = (Yq -Yp)* modInverse((Xq-Xp),p);
    }
    
    while(lamda.negative){
      lamda = Add(lamda,p);}

    //Xr = (lamda*lamda) -Xp -Xq;
    Xr =  Sub((Sub((Mul(lamda,lamda)),Xp)),Xq); 
    while(Xr.negative){
      Xr  = Add(Xr,p); }
    Xr = PwrMod(Xr,One,p);
    diffxpxr =  Sub(Xp,Xr);
    
    while(diffxpxr.negative) {
      diffxpxr = Add(diffxpxr,p); }
    Yr = Sub((Mul(lamda,diffxpxr)),Yp);  
    while(Yr.negative) {
      Yr  = Add(Yr,p); }
    Yr = PwrMod(Yr,One,p); 
    while(Xr.negative) {
      Xr  = Add(Xr,p); }
    while(Yr.negative) {
      Yr  = Add(Yr,p); }
}

DivResult Quo;
//Computes x = y + y + y .. N times using Double and Add algorithm in log(N) steps where x =(Xp,Yp), y = (Xq,Yq)
void dad (BigNum Xp,BigNum Yp ,BigNum Xq,BigNum Yq, BigNum N) {
    BigNum flag;
    flag.Num[0]  = 0;
    while (!(N.negative) &  Compare(N,StringToArray("0"))) {  
         if (!Compare(PwrMod(N,One,StringToArray("2")),One)){
            if (!Compare(flag,StringToArray("0"))){
              Xp = Xq; Yp = Yq;    //x = y;
              flag.Num[0]  = 1;
            }
           else {
             pointsum(Xp,Yp,Xq,Yq); // x = x+y;//x is the final result of addition
             Xp = Xr; Yp = Yr;
           }
         }

    pointsum(Xq,Yq,Xq,Yq); //y= y << 1;//double y
    Xq = Xr; Yq = Yr;
	Quo = DivLarge(N,StringToArray("2"));
    N = Quo.Result;
    }
    Xr = Xp; Yr =Yp; 
}

//This func set ECC parameters 
//Pcurve = 2**256 - 2**32 - 2**9 - 2**8 - 2**7 - 2**6 - 2**4 -1 # The proven prime
//#Elliptic curve: y^2 = x^3 + Acurve * x + Bcurve
void ecc_para() {
    Nb.Num[0] = 0;
    add_big = StringToArray("100000");//private key
    Nb = Add(Nb,add_big);
    
    p.Num[0] = 0;
    add_str = "115792089237316195423570985008687907853269984665640564039457584007908834671663"; 
    add_big = StringToArray(add_str);
    p = Add(p,add_big);
    
    a.Num[0] = 0;
    add_str = "0"; //a = 0 
    add_big = StringToArray(add_str);
    a = Add(a,add_big);
    
    b.Num[0] = 0;
    add_str = "7"; //b = 7 
    add_big = StringToArray(add_str);
    b = Add(b,add_big);
     

    Gx.Num[0] = 0;
    add_str = "55066263022277343669578718895168534326250603453777594175500187360389116729240";
    add_big = StringToArray(add_str);
    Gx = Add(Gx,add_big);

    Gy.Num[0] = 0;
    add_str = "32670510020758816978083085130507043184471273380659243275938904335757337482424";
    add_big = StringToArray(add_str);
    Gy = Add(Gy,add_big);
    
    BigNum pt_count;//=1;
    pt_count.Num[0] = 1;//=1;
  
    /*
     do{//Generating all the points in the curve
        //for very large primes this takes lot of time hence disabling
        //pt_count++;
        pt_count = Add(pt_count,One); 
       
        pointsum(Xr,Yr,Gx,Gy);
        //cout <<  "pt_count is " << value_number(pt_count) << endl; 
        
      }while (!(!Compare(Gx,Xr) && !Compare(Gy,Sub(p,Yr)))); 
      // pt_count = Add(pt_count,One); 
       */ 
    cout <<  "\nCurve parameters Ep(a,b) is E " << value_number(p) << "(0" << value_number(a) << "," << value_number(b) << ")" << endl; 
    cout <<  "Generator used is (" << value_number(Gx) << "," << value_number(Gy) << endl; 
    cout <<  "Private Key is " << value_number(Nb) << endl; 
}

int main()
{
    BigNum kGx,kGy,kNbGx,kNbGy,Tx,Ty;
    BigNum Pby,kPby,Pbx,kPbx;///,One;
	One.Num[0] = 1;
    srand(time(0));
  
     //set curve parameters 
      ecc_para(); 

     //Public Key Generation NbG
       dad(Gx,Gy,Gx,Gy,Nb);
       Pbx = Xr; Pby = Yr;
       cout <<  "Public Key is (" << value_number(Pbx) << "," << value_number(Pby) << ")"  << endl; 

       //Pmx = , Pmy = //Chose a plain text PM lying on the elliptic curve
       BigNum Pmx,Pmy; 
       string p_str;
	   cout<<"\nEnter the message Pmx (<p): \n";
	   cin>>p_str;
	   Pmx = StringToArray(p_str);
	   cout<<"\nEnter the message Pmy (<p): \n";
	   cin>>p_str;
	   Pmy = StringToArray(p_str);
       cout <<  "\nPlain text Chosen is (" << value_number(Pmx) << "," << value_number(Pmy) << ")"  << endl; 

      //Encryption cipher text C1 = kG
	   k = PwrMod(StringToArray(to_string(rand())),One,StringToArray("100"));
       dad(Gx,Gy,Gx,Gy,k);
       kGx = Xr; kGy = Yr;
       cout <<  "Cipher text C1 is (" << value_number(kGx) << "," << value_number(kGy) << ")"  << endl; 

      //Encryption Cipher text C2 = Pm + kNbG
       dad(Pbx,Pby,Pbx,Pby,k);
       kPbx = Xr; kPby = Yr;//kNbG Public Key
       pointsum(Pmx,Pmy,kPbx,kPby);//Pm + kNbG
       Tx = Xr; Ty = Yr;
       cout <<  "Cipher text C2 is (" << value_number(Tx) << "," << value_number(Ty) << ")"  << endl; 

       //Decryption C2 -NbC1
       dad(kGx,kGy,kGx,kGy,Nb);//NbC1
       kNbGx = Xr; kNbGy = Yr;
       pointsum (Tx,Ty,kNbGx,Sub(p,kNbGy));//-NBC1 = p - NBC1
       cout <<  "\nDeciphered Plain text is (" << value_number(Xr) << "," << value_number(Yr) << ")"  << endl; 
   return 0;
}

