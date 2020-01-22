#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void FFT1JT(int N,int S,int M,double Fr[2],double Fi[2])
{
//#pragma HLS PIPELINE
  double PI=3.14159265F;

  //bit reverse
  int K,L,I,J;
  double Ar,Ai,Cr,Ci,Wr,Wi;
  double Q;
  double tmp_A;

  int L1,L2;
  double SIG;

  for (L=1;L<=M;L++){
    L1=1<<L;
    L2=L1/2;
    Ar=1.0F;
    Ai=0.0F;
    SIG=1.0F;
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
    for(I=1;I<=N;I++){
      Fr[I-1]=Fr[I-1]/N;
      Fi[I-1]=Fi[I-1]/N;
    }
  }
}

int main()
{
int N=2;

double Fr[N],Fi[N];
int T,S;
int i,j;
int M=1;

printf("FFT start\n");
printf("number of bit M=%d\n",M);

for (i=0;i<N;i++){
  Fr[i]=1.0F;Fi[i]=0.0F;
}

  S=0;//S=0 FFT S!=0 iFFT

  for(i=0;i<N;i++){
    printf("%d   %f   %f\n",i,Fr[i],Fi[i]);
  }
  FFT1JT(N,S,M,Fr,Fi);

  printf("\n");

  for(i=0;i<N;i++){
    printf("%d   %f   %f\n",i,Fr[i],Fi[i]);
  }
  return(0);

}

