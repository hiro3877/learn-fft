#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sys/time.h>


#define N 8

double Fr[N][N],Fi[N][N];


void FFT1JT(int S,double Fr[N],double Fi[N])
{
  double PI=3.14159265358979324;

  //bit reverse
  int K,L,N2,I,J;
  double Ar,Ai,Cr,Ci,Wr,Wi;
  double Q;
  double tmp_A;
  double tmp;
  int JL;

  N2=N/2;

  L=N2;

  tmp=Fr[1];
  Fr[1]=Fr[4];
  Fr[4]=tmp;
  tmp=Fi[1];
  Fi[1]=Fi[4];
  Fi[4]=tmp;

  tmp=Fr[3];
  Fr[3]=Fr[6];
  Fr[6]=tmp;
  tmp=Fi[3];
  Fi[3]=Fi[6];
  Fi[6]=tmp;


  //DFT

  int L1,L2;
  double SIG;

  for (L=1;L<=3;L++){
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


void FFT2JS(int S,int M,double Fr[N][N],double Fi[N][N])
{
  int I,J,P,Q;
  double Xr[N],Xi[N];
  for (Q=0;Q<N;Q++){
    for(I=0;I<N;I++){
      Xr[I]=Fr[I][Q];
      Xi[I]=Fi[I][Q];
    }
    FFT1JT(S,Xr,Xi);
    for(I=0;I<N;I++){
      Fr[I][Q]=Xr[I];
      Fi[I][Q]=Xi[I];
    }
  }
  for(P=0;P<N;P++){
    for(J=0;J<N;J++){
      Xr[J]=Fr[P][J];
      Xi[J]=Fi[P][J];
    }
    FFT1JT(S,Xr,Xi);
    for(J=0;J<N;J++){
      Fr[P][J]=Xr[J];
      Fi[P][J]=Xi[J];
    }
  }
}


int main()
{
	
struct timeval s, e;
gettimeofday(&s, NULL);

int T,S;
int M=3;
int i,j;

for (i=0;i<N;i++){
	for (j=0;j<N;j++){

		if(i<3 && j<2) Fr[i][j]=1.0;Fi[i][j]=0.0;
		if(i>=3 || j>=2) Fr[i][j]=0.0;Fi[i][j]=0.0;

	}
}

for (i=0;i<N;i++){
  for (j=0;j<N;j++){
	  printf("%d %d   %lf %lf\n",i,j,Fr[i][j],Fi[i][j]);
  }
}

printf("FFT start\n");


  S=0;//S=0 FFT S!=0 IFFT
  
  gettimeofday(&s, NULL);
  
  //for(i=0;i<1;i++){

  FFT2JS(S,M,Fr,Fi);
  
  //}
  
  gettimeofday(&e, NULL);
  
for (i=0;i<N;i++){
	for (j=0;j<N;j++){
  	  printf("%d %d   %lf %lf\n",i,j,Fr[i][j],Fi[i][j]);
    }
  }
  
  printf("time = %lf\n", (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6);

  

}
