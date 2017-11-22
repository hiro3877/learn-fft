#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void fftt(int,int,double *,double *,int);
void dump(int,double *,double *);

#define NX 1024

void main()
{
  static double xr[NX],xi[NX];
  int i,m,n;

  printf("please input a number of data   N=2**M   M= : ");
  scanf("%d",&m);

  n= 1 << m;

  if(n>NX){
    printf("M is too large\n");
    exit(1);
  }






  FILE *fp;
  double f;

  fp=fopen("input.txt","r");
  if(fp==NULL){
    printf("file open error¥n");
  }

  for (i = 0; i < n; i++)
  {
      if (fscanf(fp, "%lf", &f) != 1){
       printf("file read error\n");
     }
      xr[i] = f;
      xi[i] = 0.0;
  }
  fclose(fp);









  printf("start fft\n");
  fftt(m,n,xr,xi,0);

  printf("\n\nresult\n");
  dump(n,xr,xi);

  printf("start ifft\n");
  fftt(m,n,xr,xi,1);

  printf("\n\nresult\n");
  dump(n,xr,xi);
}


void fftt(int m,int n,double xr[],double xi[],int iv)
{
  static double pi = 3.14159265358979324;
  double dp,px,temp,ti,tr,w,wi,wr;
  int i,ip,j,k,n2,p,pp;

  n2=n/2;

  //reverse bit

  j=n2;

  for(i=1;i<=n-3;i++){
    if(i<j){
      temp =xr[i];
      xr[i]=xr[j];
      xr[j]=temp;
      temp =xi[i];
      xi[i]=xi[j];
      xi[j]=temp;
    }
    k=n2;
    while(k<=j){
      j -= k;
      k /=2;
    }
    j+=k;
  }

  //switch check

  if(iv){
    temp = 1.0 / n;
    for(i=0;i<n;i++){
      xr[i] *= temp;
      xi[i] *= temp;
    }
    px = pi;
  }
  else{
    px = -pi;
  }

  //batterfly

  p=1;
  for(k=1;k<=m;k++){
    pp=p+p;
    dp=px/p;
    for(j=0;j<p;j++){
      w=j*dp;
      wr=cos(w);
      wi=sin(w);

      for(i=j;i<n;i+=pp){
        ip=i+p;
        tr=xr[ip]*wr-xi[ip]*wi;
        ti=xr[ip]*wi+xi[ip]*wr;
        xr[ip]=xr[i]-tr;
        xi[ip]=xi[i]-ti;
        xr[i]=xr[i]+tr;
        xi[i]=xi[i]+ti;
        printf("Fr[K]=Fr[%d]=%lf,Fi[K]=Fi[%d]=%lf\n",ip,xr[ip],ip,xi[ip]);
        printf("Fr[J]=Fr[%d]=%lf,Fi[J]=Fi[%d]=%lf\n",i,xr[i],i,xi[i]);
        printf("\n");
      }
    }
    p=pp;
  }

  return;
}

void dump(int n,double xr[],double xi[])
{
  int i;

  for(i=0;i<n;i++){
    printf("%d   %lf   %lf\n",i,xr[i],xi[i]);
  }

  return;

}
