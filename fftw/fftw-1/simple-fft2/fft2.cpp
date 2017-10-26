//====================================================================
// fft 2D , パワースペクトル表示
//
// program名 <画像ファイル名>
//
// (c)Copyright Spacesoft corp., 2017 All rights reserved.
//                                    Kitayama, Hiroyuki
//====================================================================
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <string.h>
#include <iostream>
#include "Cbmp.h"
#include <fftw3.h>
#include <complex.h>
#include <math.h>


//--------------------------------------------------------------------
//  マクロの宣言
#define SP_DELETE(p)        if(p) { delete p; p=NULL;}

using namespace std;

//--------------------------------------------------------------------
// read image, RGB->Grayscale, then save it
unsigned char*
getImage(const char *fname, int *inWidth, int *inHeight)
{
    unsigned char* gs = NULL;
    Cbmp inBmp;
    inBmp.loadFromFile(fname);                      // load bitmap

    int width = inBmp.getWidth();
    int height = inBmp.getAbsHeight();

    gs = new unsigned char[width*height];
    inBmp.getGSData(gs);                            // rgb to gray

    Cbmp outBmp;
    outBmp.create24Dib(width, height);
    for (int y = 0; y < height; y++)
    {
        unsigned char* p = outBmp.getScanRow(y);
        for (int x = 0; x < width; x++)
        {
            p[0] = p[1] = p[2] = gs[y*width + x];
            p += 3;
        }
    }
    outBmp.saveToFile("in.bmp");                     // save gray bmp

    *inWidth = width;
    *inHeight = height;

    return gs;
}


//--------------------------------------------------------------------
// save IFFT image
void
saveIFFTimage(const fftw_complex *ifft, const int width, const int height)
{
    Cbmp outBmp;

    outBmp.create24Dib(width, height);
    for (int y = 0; y < height; y++)
    {
        unsigned char* p = outBmp.getScanRow(y);
        for (int x = 0; x < width; x++)
        {
            p[0] = p[1] = p[2] = (unsigned char)ifft[y*width + x][0];
            p += 3;
        }
    }
    outBmp.saveToFile("ifft.bmp");                      // save gray bmp
}


//--------------------------------------------------------------------
//
// swap: 1 <-> 4, 2 <-> 3
//
//              |                           |
//         1    |    2                 4    |    3
//              |                           |
//   -----------+------------     ----------+------------
//              |                           |
//         3    |    4                 2    |    1
//              |                           |
//
void swap(double *a, double *b)
{
    double temp = *a;
    *a = *b;
    *b = temp;
}
void swap1234(double *real, int width, int height)
{
    int hHeight = height / 2;
    int hWidth = width / 2;

    for (int y = 0; y < hHeight; y++)
    {
        for (int x = 0; x < hWidth; x++)
        {
            swap(&real[width*y + x], &real[width*(hHeight + y) + hWidth + x]);
            swap(&real[width*y + hWidth + x], &real[width*(hHeight + y) + x]);
        }
    }
}


//--------------------------------------------------------------------
// power spectol 2D
double *
powerSpectol(fftw_complex *fft, int width, int height)
{
    double max, min, scale; // max/min of powerspectol

    double *out = new double[height * width];

    // scaled powerspectol : log(sqrt(real^2 + image^2))
    for (int i = 0; i < height*width; i++)
    {
        out[i] = log10(1.f + sqrt(pow(fft[i][0], 2) + pow(fft[i][1], 2)));
    }

    // normalization,  search max, min
    max = min = out[0];
    for (int i = 0; i < height*width; i++)
    {
        if (out[i] > max) max = out[i];
        if (out[i] < min) min = out[i];
    }

    // normalize, 0.0 - 1.0
    scale = (double)(1. / (max - min));
    for (int i = 0; i < height*width; i++)
    {
        out[i] = (out[i] - min) * scale;
    }

    swap1234(out, width, height);

    return out;
}


//--------------------------------------------------------------------
// save PowerSpectol
void
savePowerSpectol(fftw_complex *fft, int width, int height)
{
    double *out = powerSpectol(fft, width, height);

    Cbmp powerBmp;
    powerBmp.create24Dib(width, height);
    for (int y = 0; y < height; y++)
    {
        unsigned char* p = powerBmp.getScanRow(y);
        for (int x = 0; x < width; x++)
        {
            // 0.0 - 1.0 -> 0 - 255
            p[0] = p[1] = p[2] = (unsigned char)(out[y*width + x]*255.0);
            p += 3;
        }
    }
    powerBmp.saveToFile("power.bmp");
    SP_DELETE(out);
}

//-------------------------------------------------------------------
// main
int
main(int argc, char* argv[])
{
    unsigned char* gs = NULL;
    int width, height;

    fftw_complex *in, *fft, *ifft;
    fftw_plan plan;

    try
    {
        if (argc < 2)
            throw "引数に<入力ファイル>を指定してください.¥n";

        //画像読み込み
        if((gs = getImage(argv[1], &width, &height))==NULL)
            throw "読み込み失敗.";

        //メモリ割り付け
        in  = (fftw_complex *)fftw_malloc(height * width * sizeof(fftw_complex));
        fft = (fftw_complex *)fftw_malloc(height * width * sizeof(fftw_complex));

        // inの実部に画像を詰め込む
        for (int i = 0; i < height*width; i++)
        {
            in[i][0] = (double)gs[i];
            in[i][1] = (double)0.0;
        }
        SAFE_DELETE(gs);

        // FFT
        plan = fftw_plan_dft_2d(height, width, in, fft, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        // パワースペクトルの保存
        savePowerSpectol(fft, width, height);

        // 正規化
        for (int i = 0; i < height * width; i++)
        {
            fft[i][0] /= height * width;
            fft[i][1] /= height * width;
        }

        // IFFT
        ifft = (fftw_complex *)fftw_malloc(height * width * sizeof(fftw_complex));
        plan = fftw_plan_dft_2d(height, width, fft, ifft, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        // IFFTの結果を書き込む
        saveIFFTimage(ifft, width, height);

        // 領域の開放
        fftw_free(fft);
        fftw_free(ifft);
    }
    catch (const char* str)
    {
        cerr << str << endl;
    }
    return 0;
}
