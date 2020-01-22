#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 8

double Fr[N],Fi[N];

void FFT1JT(int S,int M,double Fr[N],double Fi[N])
{
  double PI=3.14159265358979324;

  //bit reverse
  int K,L,N2,I,J;
  double Ar,Ai,Cr,Ci,Wr,Wi;
  double Q;
  double tmp_A;

  int JL;

  N2=N/2;

  L=N2;

  for(K=1;K<=N-2;K++){
    if(K<L){
      Wr=Fr[L];
      Fr[L]=Fr[K];
      Fr[K]=Wr;
      Wi=Fi[L];
      Fi[L]=Fi[K];
      Fi[K]=Wi;
    }
    I=N2;
    while(I<L+1){
      L=L-I;
      I=I/2;
    }
    L=L+I;
  }

  //DFT

  int L1,L2;
  double SIG;

  for (L=1;L<=M;L++){
    L1=1<<L;
    L2=L1/2;
    Ar=1.0;
    Ai=0.0;
    SIG=1.0;
    if(S==0){
      SIG=-SIG;
    }
    Q=PI/L2;
    Cr=cos(Q);
    Ci=SIG*sin(Q);
    for(I=0;I<L2;I++){
      for(J=I;J<=N-1;J+=L1){
        K=J+L2;
        Wr=Fr[K]*Ar-Fi[K]*Ai;
        Wi=Fi[K]*Ar+Fr[K]*Ai;
        Fr[K]=Fr[J]-Wr;
        Fi[K]=Fi[J]-Wi;
        Fr[J]=Fr[J]+Wr;
        Fi[J]=Fi[J]+Wi;
      }
      tmp_A=Ar;
      Ar=Ar*Cr-Ai*Ci;
      Ai=Ai*Cr+tmp_A*Ci;
    }
  }

  if(S!=0){
    for(I=0;I<N;I++){
      Fr[I]=Fr[I]/N;
      Fi[I]=Fi[I]/N;
    }
  }
}

int main()
{
int T,S;
int M=3;
int i;



for (i=0;i<N;i++){
  if (i<3) Fr[i]=1;Fi[i]=0;
  if (i>=3) Fr[i]=0;Fi[i]=0;
}

printf("FFT start\n");
printf("number of bit M=%d\n",M);


  S=0;//S=0 FFT S!=0 iFFT

  for(i=0;i<N;i++){
    printf("%d   %lf   %lf\n",i,Fr[i],Fi[i]);
  }

  FFT1JT(S,M,Fr,Fi);

  printf("\n");

  for(i=0;i<N;i++){
    printf("%d   %lf   %lf\n",i,Fr[i],Fi[i]);
  }


}
