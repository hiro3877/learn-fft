/*
  n x n 二次元 DFT :
  F[k1][k2]=Σ_j1=0^n-1 Σ_j2=0^n-1
                a[j1][a2] * exp(±2*pi*i*j1*k1/n)
                          * exp(±2*pi*i*j2*k2/n)
  を計算する.
  n はデータ数で 2 の整数乗, theta は ±2*pi/n .
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void FFT2JT(int n, double theta, double ar[][128], double ai[][128])
{
    int m, mh, i1, i2, j1, j2, k1, k2;
    double w1r, w1i, w2r, w2i, w12r, w12i;
    double xr, xi, xjjr, xjji, xkjr, xkji, xjkr, xjki, xkkr, xkki;

    for (m = n; (mh = m >> 1) >= 1; m = mh) {
        for (i1 = 0; i1 < mh; i1++) {
            w1r = cos(theta * i1);
            w1i = sin(theta * i1);
            for (i2 = 0; i2 < mh; i2++) {
                w2r = cos(theta * i2);
                w2i = sin(theta * i2);
                w12r = cos(theta * (i1 + i2));
                w12i = sin(theta * (i1 + i2));
                for (j1 = i1; j1 < n; j1 += m) {
                    k1 = j1 + mh;
                    for (j2 = i2; j2 < n; j2 += m) {
                        k2 = j2 + mh;
                        xjjr = ar[j1][j2] + ar[k1][j2];
                        xjji = ai[j1][j2] + ai[k1][j2];
                        xkjr = ar[j1][j2] - ar[k1][j2];
                        xkji = ai[j1][j2] - ai[k1][j2];
                        xjkr = ar[j1][k2] + ar[k1][k2];
                        xjki = ai[j1][k2] + ai[k1][k2];
                        xkkr = ar[j1][k2] - ar[k1][k2];
                        xkki = ai[j1][k2] - ai[k1][k2];
                        ar[j1][j2] = xjjr + xjkr;
                        ai[j1][j2] = xjji + xjki;
                        xr = xjjr - xjkr;
                        xi = xjji - xjki;
                        ar[j1][k2] = w2r * xr - w2i * xi;
                        ai[j1][k2] = w2r * xi + w2i * xr;
                        xr = xkjr + xkkr;
                        xi = xkji + xkki;
                        ar[k1][j2] = w1r * xr - w1i * xi;
                        ai[k1][j2] = w1r * xi + w1i * xr;
                        xr = xkjr - xkkr;
                        xi = xkji - xkki;
                        ar[k1][k2] = w12r * xr - w12i * xi;
                        ai[k1][k2] = w12r * xi + w12i * xr;
                    }
                }
            }
        }
        theta *= 2;
    }
    /* ---- unscrambler ---- */
    i1 = 0;
    for (j1 = 1; j1 < n - 1; j1++) {
        for (k1 = n >> 1; k1 > (i1 ^= k1); k1 >>= 1);
        if (j1 < i1) {
            for (j2 = 0; j2 < n; j2++) {
                xr = ar[j1][j2];
                xi = ai[j1][j2];
                ar[j1][j2] = ar[i1][j2];
                ai[j1][j2] = ai[i1][j2];
                ar[i1][j2] = xr;
                ai[i1][j2] = xi;
            }
        }
    }
    for (j1 = 0; j1 < n; j1++) {
        i2 = 0;
        for (j2 = 1; j2 < n - 1; j2++) {
            for (k2 = n >> 1; k2 > (i2 ^= k2); k2 >>= 1);
            if (j2 < i2) {
                xr = ar[j1][j2];
                xi = ai[j1][j2];
                ar[j1][j2] = ar[j1][i2];
                ai[j1][j2] = ai[j1][i2];
                ar[j1][i2] = xr;
                ai[j1][i2] = xi;
            }
        }
    }
}

void main()
{

  double ar[128][128],ai[128][128];
  int N;
  int M=2;
  int i,j;

  N=1 << M;

  int n = N*N;

  double theta = -6.283185307179584/N;

  for (i=0;i<N;i++){
  	for (j=0;j<N;j++){

  		if(i<3 && j<2) ar[i][j]=1;ai[i][j]=0;
  		if(i>=3 || j>=2) ar[i][j]=0;ai[i][j]=0;

  	}
  }

  printf("FFT start\n");
  printf("number of data N*N N=%d\n",N);



  if(N>128){
    printf("M is too large\n");
    exit(1);
  }

    /*for(i=0;i<N;i++){
  	printf("%lf   %lf   %lf   %lf\n",Fr[0][i],Fr[1][i],Fr[2][i],Fr[3][i]);
    }*/

    for(i=0;i<N;i++){
      for(j=0;j<N;j++){
  		printf("%d  %d    %lf    %lf\n",i,j,ar[i][j],ai[i][j]);
  	}
    }

    FFT2JT(n,theta,ar,ai);

    printf("\n");

      for(i=0;i<N;i++){
      for(j=0;j<N;j++){
  		printf("%d  %d    %lf    %lf\n",i,j,ar[i][j],ai[i][j]);
  	}
    }














}
