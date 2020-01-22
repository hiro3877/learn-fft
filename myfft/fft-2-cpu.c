#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sys/time.h>

#define N 8
#define M 3
#define M2 1

double Fr[N][N],Fi[N][N];


void FFT2JT(double Fr[N][N],double Fi[N][N])
{
	
	struct timeval s, e;

  char i,j;

  double P=6.283185307179584/N;

  //bit reverse
  char K,L,I,J;

  double Cr,Ci,S1r,S1i,S3r,S3i,S5r,S5i,S7r,S7i,Tr,Ti;
  double Q;
  double Xr[N],Xi[N];

  char IR,I1,I2,I3,I4;
  char JR,J1,J2;
  char KF,K1,K2,K3,K4;
  char IP,JP,IM,KI;
  char JM,KJ,KK;
  
  gettimeofday(&s, NULL);
  
  for(i=0;i<100;i++){

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
 


  //DFT

  char L1,L2;
  char SIG=-1;

  Cr=cos(P);
  Ci=SIG*sin(P);
  Xr[0]=1;
  Xi[0]=0;

  for(I=0;I<N-1;I++){
	  Xr[I+1]=Xr[I]*Cr-Xi[I]*Ci;
	  Xi[I+1]=Xr[I]*Ci+Xi[I]*Cr;
  }

  for(L=0;L<M;L++){
	  K1=1<<(M-L);
	  K2=K1/2;
	  K3=N/K1;
	  K4=K3*2;
	  for(IP=0;IP<K2;IP++){
		  I1=IP*K4;
		  for(JP=0;JP<K2;JP++){
			  J1=JP*K4;
			  for(I=0;I<K3;I++){
				  IR=I+I1;
				  IM=IR+K3;
				  KI=I*K2;
				  for(J=0;J<K3;J++){
					  JR=J+J1;
					  JM=JR+K3;
					  KJ=J*K2;
					  KK=KI+KJ;
					  S1r=Fr[IR][JR];
					  S1i=Fi[IR][JR];
					  S3r=Fr[IR][JM]*Xr[KJ]-Fi[IR][JM]*Xi[KJ];
					  S3i=Fi[IR][JM]*Xr[KJ]+Fr[IR][JM]*Xi[KJ];
					  S5r=Fr[IM][JR]*Xr[KI]-Fi[IM][JR]*Xi[KI];
					  S5i=Fi[IM][JR]*Xr[KI]+Fr[IM][JR]*Xi[KI];
					  S7r=Fr[IM][JM]*Xr[KK]-Fi[IM][JM]*Xi[KK];
					  S7i=Fi[IM][JM]*Xr[KK]+Fr[IM][JM]*Xi[KK];
					  Fr[IR][JR]=S1r+S3r+S5r+S7r;
					  Fi[IR][JR]=S1i+S3i+S5i+S7i;
					  Fr[IR][JM]=S1r-S3r+S5r-S7r;
					  Fi[IR][JM]=S1i-S3i+S5i-S7i;
					  Fr[IM][JR]=S1r+S3r-S5r-S7r;
					  Fi[IM][JR]=S1i+S3i-S5i-S7i;
					  Fr[IM][JM]=S1r-S3r-S5r+S7r;
					  Fi[IM][JM]=S1i-S3i-S5i+S7i;
				  }
			  }
		  }
	  }
  }
  
  }
  
  gettimeofday(&e, NULL);
  printf("time = %lf\n", (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6);
  
  
}



int main()
{
	
	

char T;
char i,j;


for (i=0;i<N;i++){
	for (j=0;j<N;j++){

		if(i<3 && j<2) Fr[i][j]=1.0;Fi[i][j]=0.0;
		if(i>=3 || j>=2) Fr[i][j]=0.0;Fi[i][j]=0.0;
	}
}


printf("FFT start\n");
printf("number of data N*N N=%d\n",N);

  //S=0 FFT S!=0 IFFT

/*for (i=0;i<N;i++){
	for (j=0;j<N;j++){
  	  printf("%d %d   %lf %lf\n",i,j,Fr[i][j],Fi[i][j]);
    }
  }*/
  
  
  

  FFT2JT(Fr,Fi);
  
  
  


/*for (i=0;i<N;i++){
	for (j=0;j<N;j++){
  	  printf("%d %d   %lf %lf\n",i,j,Fr[i][j],Fi[i][j]);
    }
  }*/
  
  
  
    return 0;
	
	
}

