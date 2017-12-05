#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sys/time.h>
#include <sys/resource.h>

#define N 1024

double Fr[N][N],Fi[N][N];


double get_time(void)
 {
 struct timeval tv;
 gettimeofday(&tv, NULL);
 //ミリ秒を計算
return ((double)(tv.tv_sec)*1000 + (double)(tv.tv_usec)*0.001);
 }

void FFT1JT(int S,int M,double Fr[],double Fi[])
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

void FFT2JS(int S,int M,double Fr[][N],double Fi[][N])
{
  int I,J,P,Q;
  double Xr[N],Xi[N];
  for (Q=0;Q<N;Q++){
    for(I=0;I<N;I++){
      Xr[I]=Fr[I][Q];
      Xi[I]=Fi[I][Q];
    }
    FFT1JT(S,M,Xr,Xi);
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
    FFT1JT(S,M,Xr,Xi);
    for(J=0;J<N;J++){
      Fr[P][J]=Xr[J];
      Fi[P][J]=Xi[J];
    }
  }
}

int main()
{

int T,S;
int M=log2(N);
int i,j;

char filename[20] = {};

FILE *fp;
FILE *fp1;
FILE *fp2;

double f;
double starttime1,endtime1,starttime2,endtime2;

printf("please input a filename  :  ");
scanf("%s",filename);

fp=fopen(filename,"r");
if(fp==NULL){
  printf("file open error¥n");
}

for (i=0;i<N;i++){
  for (j=0;j<N;j++){

    if (fscanf(fp, "%lf", &f) != 1){
     printf("file read error\n");
   }
    Fr[i][j] = f;
    Fi[i][j] = 0.0;
  }
}
fclose(fp);

printf("FFT start\n");
printf("number of bit M=%d\n",M);


  S=0;//S=0 FFT S!=0 IFFT

  starttime1=get_time();

  FFT2JS(S,M,Fr,Fi);

  endtime1=get_time();


  fp1=fopen("fft-1d.txt","w");
  if(fp1==NULL){
    printf("file open error\n");
  }

  for (i=0;i<N;i++){
    for (j=0;j<N;j++){
      fprintf(fp1,"%lf, %lf\n",Fr[i][j],Fi[i][j]);
    }
  }
  fclose(fp1);


  S=1;

  starttime2=get_time();

  FFT2JS(S,M,Fr,Fi);

  endtime2=get_time();

  fp2=fopen("ifft-1d.txt","w");
  if(fp2==NULL){
    printf("file open error\n");
  }

  for (i=0;i<N;i++){
    for (j=0;j<N;j++){
      fprintf(fp2,"%lf, %lf\n",Fr[i][j],Fi[i][j]);
    }
  }
  fclose(fp2);


  printf("FFT Calculation time is %lf\n",endtime1-starttime1);

  printf("IFFT Calculation time is %lf\n",endtime2-starttime2);

}
