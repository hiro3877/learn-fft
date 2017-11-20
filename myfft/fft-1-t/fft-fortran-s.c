#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void FFT1JF(int S,int M,double Fr[],double Fi[])
{
  int N;
  N=1<<M;

  double pi=3.14159265358979324;

  //bit reverse
  int K,L,N2,I,J,i;
  double Ar,Ai,Cr,Ci,Wr,Wi;
  double Q;

  N2=N/2;

  L=N2;

  for(K=1;K<=N-2;K++){
    if(K<L){
      Wr=Fr[L+1-1];
      Fr[L+1-1]=Fr[K+1-1];
      Fr[K+1-1]=Wr;
      Wi=Fi[L+1-1];
      Fi[L+1-1]=Fi[K+1-1];
      Fi[K+1-1]=Wi;
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

  for (L=0;L<M;L++){
    L1=1<<(M-L+1);
	
    L2=L1/2;
    Ar=1;
    Ai=0;
    SIG=1;
    if(S=0){
      SIG=-SIG;
    }
    Q=pi/L2;
    Cr=cos(Q);
    Ci=SIG*sin(Q);
    for(I=0;I<L2;I++){
      for(J=0;J<N;J=J+L1){
        K=J+L2;
        Wr=Fr[J]-Fr[K];
        Wr=Fi[J]-Fi[K];
        Fr[J]=Fr[J]+Fr[K];
        Fi[J]=Fi[J]+Fi[K];
        Fr[K]=Wr*Ar+Wi*Ai;
        Fi[K]=Wr*Ai+Wi*Ar;
      }
      Ar=Ar*Cr+Ai+Ci;
      Ai=Ai*Cr+Ar+Ci;
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

double Fr[512],Fi[512];
int T,S,N;
int M=3;
int i;

N=1 << M;

for (i=0;i<N;i++){
  if (i<3) Fr[i]=1;Fi[i]=0;
  if (i>=3) Fr[i]=0;Fi[i]=0;
}

printf("FFT start\n");
printf("number of bit M=%d\n",M);



if(N>512){
  printf("M is too large\n");
  exit(1);
}

  S=0;//S=0 FFT S=0 iFFT

  for(i=0;i<N;i++){
    printf("%d   %lf   %lf\n",i,Fr[i],Fi[i]);
  }

  FFT1JF(S,M,Fr,Fi);

    for(i=0;i<N;i++){
    printf("%d   %lf   %lf\n",i,Fr[i],Fi[i]);
  }


}
