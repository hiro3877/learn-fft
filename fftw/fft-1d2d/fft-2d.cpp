#include <stdio.h>
#include<stdint.h>
#include<stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <fftw3.h>
#include <complex.h>

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


    FILE *fp;
    FILE *fp1;
    FILE *fp2;
    int i,j,k;
    int WID,HEI;
    char filename[20] = {};


    fftw_plan plan;
    fftw_complex *ibuf, *obuf;

    double f;

    printf("please input a number of width  :  ");
    scanf("%d",&WID);

    printf("please input a number of heith  :  ");
    scanf("%d",&HEI);

    printf("please input a filename  :  ");
    scanf("%s",filename);


    fp=fopen(filename,"r");
    if(fp==NULL){
      printf("file open error¥n");
    }

    ibuf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*WID*HEI);
    obuf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*WID*HEI);

    for (i = 0; i < WID*HEI; i++)
    {
        if (fscanf(fp, "%lf", &f) != 1){
         printf("file read error\n");
       }
        ibuf[i][0] = f;
        ibuf[i][1] = 0.0;
    }
    fclose(fp);



    starttime1 = getrusage_sec();

    plan = fftw_plan_dft_2d(HEI,WID, ibuf, obuf, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);	

    endtime1 = getrusage_sec();



    printf("Calculation time is %lf\n",endtime1-starttime1);


    fp1=fopen("fft-2d.txt","w");
    if(fp1==NULL){
      printf("file open error\n");
    }

    for (i = 0; i < WID*HEI; i++)
    {
      fprintf(fp1,"%lf, %lf\n",obuf[i][0],obuf[i][1]);
    }

    fclose(fp1);


    starttime2 = getrusage_sec();

    plan = fftw_plan_dft_2d(HEI,WID, obuf, ibuf, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    for(k=0;k<2;k++){
    	for(i=0;i<WID*HEI;i++){
    		ibuf[i][k]=ibuf[i][k]/(WID*HEI);
    	}
    }

    endtime2 = getrusage_sec();


    printf("Calculation time is %lf\n",endtime2-starttime2);


    fp2=fopen("ifft-2d.txt","w");
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
