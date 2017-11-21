#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void FFT2JT(int S,int M,double Fr[][128],double Fi[][128])
{
  int N;
  N=1<<M;
  
  int i,j;
  
  /*for(i=0;i<N;i++){
    for(j=0;j<N;j++){
		printf("%d  %d    %lf    %lf\n",i,j,Fr[i][j],Fi[i][j]);
	}
  }*/
	

  double P=6.283185307179584/N;

  //bit reverse
  int K,L,M2,I,J;
  
  double Cr,Ci,S1r,S1i,S3r,S3i,S5r,S5i,S7r,S7i,Tr,Ti;
  double Q;
  double Xr[128],Xi[128];
  
  int IR,I1,I2,I3,I4;
  int JR,J1,J2;
  int KF,K1,K2,K3,K4;
  int IP,JP,IM,KI;
  int JM,KJ,KK;

  M2=M/2;
  
  for (L=1;L<=M2;L++){
	  I1=1<<(L-1);
	  I2=N/I1;
	  I3=I2/2;
	  I4=I3/2;
	  for(I=0;I<=I1-1;I++){
		  JR=I*I2;
		  J1=JR+I1;
		  J2=JR+I3;
		  for (J=0;J<=I4-1;J++){
			  KF=(J/I1)*I1+J;
			  K1=KF+J1+1;
			  K2=KF+J2+1;
			  for (K=1;K<=N;K++){
				  Tr=Fr[K-1][K1-1];
				  Ti=Fi[K-1][K1-1];
				  Fr[K-1][K1-1]=Fr[K-1][K2-1];
				  Fi[K-1][K1-1]=Fi[K-1][K2-1];
				  Fr[K-1][K2-1]=Tr;
				  Fi[K-1][K2-1]=Ti;
			  }
			  for (K=1;K<=N;K++){
				  Tr=Fr[K1-1][K-1];
				  Ti=Fi[K1-1][K-1];
				  Fr[K1-1][K-1]=Fr[K2-1][K-1];
				  Fi[K1-1][K-1]=Fi[K2-1][K-1];
				  Fr[K2-1][K-1]=Tr;
				  Fi[K2-1][K-1]=Ti;
			  }
		  }
	  } 
  }


  //DFT

  int L1,L2;
  double SIG;

  SIG=1.0;
  if(S=0) SIG=-SIG;
  Cr=cos(P);
  Ci=SIG*sin(P);
  Xr[0]=1;
  Xi[0]=0;
  
  for(I=1;I<=N-1;I++){
	  Xr[I]=Xr[I-1]*Cr+Xi[I-1]*Ci;
	  Xi[I]=Xr[I-1]*Ci+Xi[I-1]*Cr;
	  
  }
  
  for(L=1;L<=M;L++){
	  K1=1<<(M+1-L);
	  K2=K1/2;
	  K3=N/K1;
	  K4=K3*2;
	  for(IP=1;IP<=K2;IP++){
		  I1=(IP-1)*K4;
		  for(JP=1;JP<=K2;JP++){
			  J1=(JP-1)*K4;
			  for(I=1;I<=K3;I++){
				  IR=I+I1;
				  IM=IR+K3;
				  KI=(I-1)*K2+1;
				  for(J=1;J<=K3;J++){
					  JR=J+J1;
					  JM=JR+K3;
					  KJ=(J-1)*K2+1;
					  KK=KI+KJ-1;
					  S1r=Fr[IR-1][JR-1];
					  S1i=Fi[IR-1][JR-1];
					  S3r=Fr[IR-1][JM-1]*Xr[KJ]+Fi[IR-1][JM-1]*Xi[KJ];
					  S3i=Fi[IR-1][JM-1]*Xr[KJ]+Fr[IR-1][JM-1]*Xi[KJ];
					  S5r=Fr[IM-1][JR-1]*Xr[KI]+Fi[IM-1][JR-1]*Xi[KI];
					  S5i=Fi[IM-1][JR-1]*Xr[KI]+Fr[IM-1][JR-1]*Xi[KI];
					  S7r=Fr[IM-1][JM-1]*Xr[KK]+Fi[IM-1][JM-1]*Xi[KK];
					  S7i=Fi[IM-1][JM-1]*Xr[KK]+Fr[IM-1][JM-1]*Xi[KK];
					  Fr[IR-1][JR-1]=S1r+S3r+S5r+S7r;
					  Fi[IR-1][JR-1]=S1i+S3i+S5i+S7i;
					  Fr[IR-1][JM-1]=S1r-S3r+S5r-S7r;
					  Fi[IR-1][JM-1]=S1i-S3i+S5i-S7i;
					  Fr[IM-1][JR-1]=S1r+S3r-S5r-S7r;
					  Fi[IM-1][JR-1]=S1i+S3i-S5i-S7i;
					  Fr[IM-1][JM-1]=S1r-S3r-S5r+S7r;
					  Fi[IM-1][JM-1]=S1i-S3i-S5i+S7i;
				  }
			  }
		  }
	  }
  }
  
}

int main()
{

double Fr[128][128],Fi[128][128];
int T,S,N;
int M=2;
int i,j;

N=1 << M;

for (i=0;i<N;i++){
	for (j=0;j<N;j++){
	
		if(i<3 && j<2) Fr[i][j]=1;Fi[i][j]=0; 
		if(i>=3 || j>=2) Fr[i][j]=0;Fi[i][j]=0; 
		
	}
}

printf("FFT start\n");
printf("number of data N*N N=%d\n",N);



if(N>128){
  printf("M is too large\n");
  exit(1);
}

  S=0;//S=0 FFT S!=0 iFFT

  /*for(i=0;i<N;i++){
	printf("%lf   %lf   %lf   %lf\n",Fr[0][i],Fr[1][i],Fr[2][i],Fr[3][i]);
  }*/
  
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
		printf("%d  %d    %lf    %lf\n",i,j,Fr[i][j],Fi[i][j]);
	}
  }

  FFT2JT(S,M,Fr,Fi);
  
  printf("\n");
  
    for(i=0;i<N;i++){
    for(j=0;j<N;j++){
		printf("%d  %d    %lf    %lf\n",i,j,Fr[i][j],Fi[i][j]);
	}
  }


}
