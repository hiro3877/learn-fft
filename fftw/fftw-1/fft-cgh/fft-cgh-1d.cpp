#include <stdio.h>
#include<stdint.h>
#include<stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <fftw3.h>
#include <complex.h>

#define WID 4096
#define HEI 4096


double getrusage_sec()
{
    struct rusage t;
    struct timeval tv;

    getrusage(RUSAGE_SELF,&t);
    tv = t.ru_utime;

    return tv.tv_sec + (double)tv.tv_usec*1e-6;
}


int main(){

    double starttime1,endtime1;
    double starttime2,endtime2;
    double starttime3,endtime3;
    double starttime4,endtime4;



    FILE *fp;
    FILE *fp1;
    FILE *fp2;
    int i,j,k;



    fftw_plan plan;
    fftw_complex *ibuf, *obuf,*tbuf;

    double f;



// ファイル読み込み

    fp=fopen("data2000.txt","r");
    if(fp==NULL){
      printf("file open error¥n");
    }

    ibuf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*WID*HEI);
    tbuf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*WID*HEI);
    obuf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*WID*HEI);

    for (i = 0; i < WID*HEI; i++)
    {
        if (fscanf(fp, "%lf", &f) != 1){
         printf("file read error\n");
       }
        ibuf[i][0] = f;
        //printf("%lf\n",ibuf[i][0]);
        ibuf[i][1] = 0.0;
    }
    fclose(fp);




//FFT横

    starttime1 = getrusage_sec();

    for(i=0;i<WID;i++){

    plan = fftw_plan_dft_1d(WID, &ibuf[i*WID], &tbuf[i*WID], FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    }

    endtime1 = getrusage_sec();

    printf("Calculation time is %lf\n",endtime1-starttime1);



//入れ替え
    double tmp;

    for(k=0;k<2;k++){

    for(i=0;i<HEI;i++){
      for(j=i;j<WID;j++){
        tmp=tbuf[i*WID+j][k];
        tbuf[i*WID+j][k]=tbuf[j*WID+i][k];
        tbuf[j*WID+i][k]=tmp;
      }
    }

    }



//FFT縦


    starttime2 = getrusage_sec();

    for(i=0;i<WID;i++){

    plan = fftw_plan_dft_1d(WID, &tbuf[i*WID], &obuf[i*WID], FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    }

    endtime2 = getrusage_sec();

    printf("Calculation time is %lf\n",endtime2-starttime2);


//入れ替え

    for(k=0;k<2;k++){

    for(i=0;i<HEI;i++){
      for(j=i;j<WID;j++){
        tmp=obuf[i*WID+j][k];
        obuf[i*WID+j][k]=obuf[j*WID+i][k];
        obuf[j*WID+i][k]=tmp;
      }
    }

  }


    fp1=fopen("ffta-1d.txt","w");
    if(fp1==NULL){
      printf("file open error\n");
    }

    for (i = 0; i < WID*HEI; i++)
    {
      fprintf(fp1,"%lf, %lf\n",obuf[i][0],obuf[i][1]);
    }

    fclose(fp1);







//IFFT横

    starttime3 = getrusage_sec();


    for(i=0;i<WID;i++){

    plan = fftw_plan_dft_1d(WID, &obuf[i*WID], &tbuf[i*WID], FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    }

    endtime3 = getrusage_sec();

    printf("Calculation time is %lf\n",endtime3-starttime3);


//入れ替え

    for(k=0;k<2;k++){

    for(i=0;i<HEI;i++){
      for(j=i;j<WID;j++){
        tmp=tbuf[i*WID+j][k];
        tbuf[i*WID+j][k]=tbuf[j*WID+i][k];
        tbuf[j*WID+i][k]=tmp;
      }
    }

    }


//IFFT縦

    starttime4 = getrusage_sec();

    for(i=0;i<WID;i++){

    plan = fftw_plan_dft_1d(WID, &tbuf[i*WID], &ibuf[i*WID], FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    }

    endtime4 = getrusage_sec();

    printf("Calculation time is %lf\n",endtime4-starttime4);


    for(k=0;k<2;k++){

    for(i=0;i<HEI;i++){
      for(j=i;j<WID;j++){
        tmp=ibuf[i*WID+j][k];
        ibuf[i*WID+j][k]=ibuf[j*WID+i][k];
        ibuf[j*WID+i][k]=tmp;
      }
    }

  }


    fp2=fopen("iffta-1d.txt","w");
    if(fp2==NULL){
      printf("file open error\n");
    }

    for (i = 0; i < WID*HEI; i++)
    {
      fprintf(fp2,"%lf, %lf\n",ibuf[i][0],ibuf[i][1]);
    }

    fclose(fp2);




    fftw_free(ibuf);
    fftw_free(obuf);


    return 0;



}
