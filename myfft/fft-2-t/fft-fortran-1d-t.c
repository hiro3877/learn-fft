#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sys/time.h>
#include <sys/resource.h>

#define N_DATA 128

double Fr[N_DATA*N_DATA],Fi[N_DATA*N_DATA];
//double Tr[N_DATA*N_DATA],Ti[N_DATA*N_DATA];

double getrusage_sec()
{
    struct rusage t;
    struct timeval tv;

    getrusage(RUSAGE_SELF,&t);
    tv = t.ru_utime;

    return tv.tv_sec + (double)tv.tv_usec*1e-6;
}

void FFT1JT(int i,int S,int M,double Tr[N_DATA],double Ti[N_DATA])
{

  double PI=3.14159265358979324;

  //printf("test\n");

  //bit reverse
  int K,L,N2,I,J;
  double Ar,Ai,Cr,Ci,Wr,Wi;
  double Q;
  int j;

  N2=N_DATA/2;

  L=N2;

  for(K=1;K<=N_DATA-2;K++){
    if(K<L){
      Wr=Tr[L];
      Tr[L]=Tr[K];
      Tr[K]=Wr;
      Wi=Ti[L];
      Ti[L]=Ti[K];
      Ti[K]=Wi;
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
    L1=1<<(L+1);
    L2=L1/2;
    Ar=1.0;
    Ai=0.0;
    SIG=1.0;
    if(S=0){
      SIG=-SIG;
    }
    Q=PI/L2;
    Cr=cos(Q);
    Ci=SIG*sin(Q);
    for(I=0;I<L2;I++){
      for(J=I;J<N_DATA-1;J=J+L1){
        K=(J+1)+L2;
        Wr=Tr[K]*Ar+Ti[K]*Ai;
        Wi=Ti[K]*Ar+Tr[K]*Ai;
        Tr[K]=Tr[J]-Wr;
        Ti[K]=Ti[J]-Wi;
        Tr[J]=Tr[J]+Wr;
        Ti[J]=Ti[J]+Wi;
      }
      Ar=Ar*Cr+Ai+Ci;
      Ai=Ai*Cr+Ar+Ci;
    }
  }

  if(S!=0){
    for(I=1;I<=N_DATA;I++){
      Tr[I-1]=Tr[I-1]/N_DATA;
      Ti[I-1]=Ti[I-1]/N_DATA;
    }
  }

  for(j=0;j<N_DATA;j++){
    Fr[i*N_DATA+j]=Tr[j];
    Fi[i*N_DATA+j]=Ti[j];
  }

  for(j=0;j<N_DATA;j++){
    Tr[j]=0;
    Ti[j]=0;
  }


}

int main()
{

int T,S;
int M=log2(N_DATA);
int i,j;

char filename[20] = {};

FILE *fp;
FILE *fp1;

double f;
double starttime1,endtime1;

printf("please input a filename  :  ");
scanf("%s",filename);

fp=fopen(filename,"r");
if(fp==NULL){
  printf("file open error¥n");
}

for (i=0;i<N_DATA*N_DATA;i++){

  if (fscanf(fp, "%lf", &f) != 1){
   printf("file read error\n");
 }
 Fr[i] = f;
 Fi[i] = 0.0;
}
fclose(fp);

/*for(i=0;i<N_DATA*N_DATA;i++){
  printf("%lf\n",Fr[i]);
}*/

printf("FFT start\n");
printf("number of bit M=%d\n",M);


  S=0;//S=0 FFT S!=0 iFFT

  /*for(i=0;i<N_DATA;i++){
    printf("%d   %lf   %lf\n",i,Fr[i],Fi[i]);
  }*/


  //FFT縦

  starttime1 = getrusage_sec();

  for(i=0;i<N_DATA;i++){
  //printf("Fr=%lf,Fi=%lf\n",Fr[i*N_DATA],Fi[i*N_DATA]);
  FFT1JT(i,S,M,Fr[i*N_DATA],&Fi[i*N_DATA]);
  printf("Fr=%lf,Fi=%lf\n",Fr[i*N_DATA],Fi[i*N_DATA]);

}

  //入れ替え

  double tmp;

  for(i=0;i<N_DATA;i++){
    for(j=i;j<N_DATA;j++){
      tmp=Fr[i*N_DATA+j];
      Fr[i*N_DATA+j]=Fr[j*N_DATA+i];
      Fr[j*N_DATA+i]=tmp;
      tmp=Fi[i*N_DATA+j];
      Fi[i*N_DATA+j]=Fi[j*N_DATA+i];
      Fi[j*N_DATA+i]=tmp;
    }
  }


//FFT横
  for(i=0;i<N_DATA;i++){

  FFT1JT(i,S,M,&Fr[i*N_DATA],&Fi[i*N_DATA]);

}


//入れ替え

for(i=0;i<N_DATA;i++){
  for(j=i;j<N_DATA;j++){
    tmp=Fr[i*N_DATA+j];
    Fr[i*N_DATA+j]=Fr[j*N_DATA+i];
    Fr[j*N_DATA+i]=tmp;
    tmp=Fi[i*N_DATA+j];
    Fi[i*N_DATA+j]=Fi[j*N_DATA+i];
    Fi[j*N_DATA+i]=tmp;
  }
}


endtime1 = getrusage_sec();



  fp1=fopen("fft-1d.txt","w");
  if(fp1==NULL){
    printf("file open error\n");
  }

  for (i=0;i<N_DATA*N_DATA;i++){
    fprintf(fp1,"%lf, %lf\n",Fr[i],Fi[i]);
  }
  fclose(fp1);

  printf("Calculation time is %lf\n",endtime1-starttime1);

  /*printf("\n");

  for(i=0;i<N;i++){
    printf("%d   %lf   %lf\n",i,Fr[i],Fi[i]);
  }*/


}
