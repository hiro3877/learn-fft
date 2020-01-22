#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 8
#define M 3
#define M2 1
#define N2 4

double Fr[N][N],Fi[N][N];

void FFT1JT(double Fr[],double Fi[])
{

  double PI=3.14159265358979324;

  //bit reverse
  int K,L,I,J;
  double Ar,Ai,Cr,Ci,Wr,Wi;
  double Q;
  double tmp_A;

  int JL;

  
  /*
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
  }*/

  //DFT

  int L1,L2;
  double SIG;

  for (L=1;L<=M;L++){
    L1=1<<L;
    L2=L1/2;
    Ar=1.0;
    Ai=0.0;
    SIG=1.0;
    SIG=-SIG;
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
}

void FFT2JS(double Fr[N][N],double Fi[N][N])
{
  int I,J,P,Q;
  double Xr[N],Xi[N];
  
  
  
  //bit reverse
  int K,L;

  int IR,I1,I2,I3,I4;
  int JR,J1,J2;
  int KF,K1,K2,K3,K4;
  int IP,JP,IM,KI;
  int JM,KJ,KK;
  
  double Tr,Ti;

  for (L=0;L<M2;L++){
	  I1=1<<L;
	  I2=N/I1;
	  I3=I2/2;
	  I4=I3/2;
	  for(I=0;I<I1;I++){
		  JR=I*I2;
		  J1=JR+I1;
		  J2=JR+I3;
		  for (J=0;J<I4;J++){
			  KF=(J/I1)*I1+J;
			  K1=KF+J1;
			  K2=KF+J2;
			  for (K=0;K<N;K++){
				  Tr=Fr[K][K1];
				  Ti=Fi[K][K1];
				  Fr[K][K1]=Fr[K][K2];
				  Fi[K][K1]=Fi[K][K2];
				  Fr[K][K2]=Tr;
				  Fi[K][K2]=Ti;
			  }
        for (K=0;K<N;K++){
          Tr=Fr[K1][K];
          Ti=Fi[K1][K];
          Fr[K1][K]=Fr[K2][K];
          Fi[K1][K]=Fi[K2][K];
          Fr[K2][K]=Tr;
          Fi[K2][K]=Ti;
        }
		  }
	  }
  }

  
  
  
  for (Q=0;Q<N;Q++){
    for(I=0;I<N;I++){
      Xr[I]=Fr[I][Q];
      Xi[I]=Fi[I][Q];
    }
    FFT1JT(Xr,Xi);
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
    FFT1JT(Xr,Xi);
    for(J=0;J<N;J++){
      Fr[P][J]=Xr[J];
      Fi[P][J]=Xi[J];
    }
  }
}

int main()
{

int T;
int i,j;


for (i=0;i<N;i++){
	for (j=0;j<N;j++){

		if(i<3 && j<2) Fr[i][j]=1.0;Fi[i][j]=0.0;
		if(i>=3 || j>=2) Fr[i][j]=0.0;Fi[i][j]=0.0;
	}
}



printf("FFT start\n");
printf("number of bit M=%d\n",M);


//S=0 FFT S!=0 IFFT

for (i=0;i<N;i++){
	for (j=0;j<N;j++){
  	  printf("%d %d   %lf %lf\n",i,j,Fr[i][j],Fi[i][j]);
    }
  }


  FFT2JS(Fr,Fi);

  for (i=0;i<N;i++){
	for (j=0;j<N;j++){
  	  printf("%d %d   %lf %lf\n",i,j,Fr[i][j],Fi[i][j]);
    }
  }


}
