#include <stdio.h>
#include<stdint.h>
#include<stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <fftw3.h>
#include <complex.h>

#define WID 2000
#define HEI 2000

/*#pragma pack(push,1)
typedef struct tagBITMAPFILEHEADER
{
  unsigned short bfType;
	uint32_t  bfSize;
	unsigned short bfReserved1;
	unsigned short bfReserved2;
	uint32_t  bf0ffBits;
}BITMAPFILEHEADER;
#pragma pack(pop)

typedef struct tagBITMAPINFOHEADER
{
	uint32_t  biSize;
	int32_t	biWidth;
	int32_t	biHeight;
	unsigned short  biPlanes;
	unsigned short  biBitCount;
	uint32_t   biCompression;
	uint32_t   biSizeImage;
	int32_t	 biXPelsPerMeter;
	int32_t	 biYPelsPerMeter;
	uint32_t   biCirUsed;
	uint32_t   biCirImportant;
}BITMAPINFOHEADER;

typedef struct tagRGBQUAD
{
	unsigned char  rgbBlue;
	unsigned char  rgbGreen;
	unsigned char  rgbRed;
	unsigned char  rgbReserved;
}RGBQUAD;

typedef struct tagBITMAPINFO
{
	BITMAPINFOHEADER bmiHeader;
	RGBQUAD          bmiColors[1];
}BITMAPINFO;
*/

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

    /*BITMAPFILEHEADER    BmpFileHeader;
    BITMAPINFOHEADER    BmpInfoHeader;
    RGBQUAD             RGBQuad[256];



    BmpFileHeader.bfType                =19778;
    BmpFileHeader.bfSize                =14+40+1024+(WID*HEI);
    BmpFileHeader.bfReserved1           =0;
    BmpFileHeader.bfReserved2           =0;
    BmpFileHeader.bf0ffBits             =14+40+1024;

    BmpInfoHeader.biSize                =40;
    BmpInfoHeader.biWidth               =WID;
    BmpInfoHeader.biHeight              =HEI;
    BmpInfoHeader.biPlanes              =1;
    BmpInfoHeader.biBitCount            =8;     //256髫手ｪｿ
    BmpInfoHeader.biCompression         =0L;
    BmpInfoHeader.biSizeImage           =0L;
    BmpInfoHeader.biXPelsPerMeter       =0L;
    BmpInfoHeader.biYPelsPerMeter       =0L;
    BmpInfoHeader.biCirUsed             =0L;
    BmpInfoHeader.biCirImportant        =0L;*/


    FILE *fp;
    FILE *fp1;
    FILE *fp2;
    int i,j,k;

    /*for(i=0;i<256;i++){
       RGBQuad[i].rgbBlue                =i;
       RGBQuad[i].rgbGreen               =i;
       RGBQuad[i].rgbRed                 =i;
       RGBQuad[i].rgbReserved            =0;
    }*/


    fftw_plan plan;
    fftw_complex *ibuf, *obuf;

    double f;


    fp=fopen("data2000.txt","r");
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
        //printf("%lf\n",ibuf[i][0]);
        ibuf[i][1] = 0.0;
    }
    fclose(fp);



    starttime1 = getrusage_sec();

    plan = fftw_plan_dft_1d(WID*HEI, ibuf, obuf, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    endtime1 = getrusage_sec();

    /*for (i = 0; i < WID*HEI; i++){
      printf("%lf, %lf\n",obuf[i][0],obuf[i][1]);
    }*/


    printf("Calculation time is %lf\n",endtime1-starttime1);


    fp1=fopen("fft-1d.txt","w");
    if(fp1==NULL){
      printf("file open error\n");
    }

    for (i = 0; i < WID*HEI; i++)
    {
      fprintf(fp1,"%lf, %lf\n",obuf[i][0],obuf[i][1]);
    }

    fclose(fp1);


    starttime2 = getrusage_sec();

    plan = fftw_plan_dft_1d(WID*HEI, obuf, ibuf, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    endtime2 = getrusage_sec();


  /*  for (int i = 0; i < WID*HEI; i++)
    {
        ibuf[i][0] /= (double)WID*HEI;
        ibuf[i][1] /= (double)WID*HEI;
    }*/


    /*for (i = 0; i < WID*HEI; i++){
      printf("%lf, %lf\n",ibuf[i][0],ibuf[i][1]);
    }*/

    printf("Calculation time is %lf\n",endtime2-starttime2);


    fp2=fopen("ifft-1d.txt","w");
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

/*
    plan = fftw_plan_dft_1d((int)WID*HEI, ibuf, obuf, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);


    FILE *fp1;
    fp1=fopen("cgh_root.bmp","wb");
    if(fp1==NULL){
      printf("file open error¥n");
    }


    FILE *fp2;
    fp2=fopen("cgh.txt","w");
    if(fp2==NULL){
      printf("file open error¥n");
    }

    for(i=0;i<WID*HEI;i++){
      fprintf(fp2,"%d¥n",img[i]);
    }


    fwrite(&BmpFileHeader, sizeof(BmpFileHeader) , 1 ,fp1);
    fwrite(&BmpInfoHeader, sizeof(BmpInfoHeader) , 1 ,fp1);
    fwrite(&RGBQuad[0], sizeof(RGBQuad[0]) , 256 ,fp1);
    fwrite(img,sizeof(unsigned char),WID*HEI,fp1);
    //fwrite(img,sizeof(unsigned char),1,fp2);


    printf("Calculation time is %lf¥n",endtime1-starttime1);

*/

    return 0;



}
