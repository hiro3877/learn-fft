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

#include <stdlib.h>
#include <assert.h>
#include <sys/stat.h>           // for SIDBA
#include "bitmapStruct.h"


//--------------------------------------------------------------------
//  マクロの宣言
#define SP_DELETE(p)        if(p) { delete p; p=NULL;}

using namespace std;


//-------------------------------------------------------------------
// コンストラクタ
Cbmp::Cbmp()
: mPdib(NULL), mPbitmap(NULL), mDibSize(0), mRowPitch(0),
                        mPixelPitch(0), mImageSize(0), mAbsHeight(0)
{
    assert(sizeof(char) ==1);
    assert(sizeof(short)==2);
    assert(sizeof(int)  ==4);
}


//-------------------------------------------------------------------
// デストラクタ
Cbmp::~Cbmp()
{
    SAFE_FREE(mPdib);                                // free bmp
}


//================ vvvvvv private vvvvvv ============================

//-------------------------------------------------------------------
// read bitmap file header
//
// return true :0
//        false:!0=error #
int
Cbmp::readHeader(FILE* fp)
{
    if(fread(&mBmpFileHdr, sizeof(bmpFileHdr), 1, fp)!=1)
        return -1;

    if(mBmpFileHdr.bfType!='B'+'M'*256)
        return -2;                                  // not bitmap file

    return 0;
}


//-------------------------------------------------------------------
// read bitmap body
int
Cbmp::readDib(FILE* fp)
{
    if(fread(mPdib , mDibSize, 1, fp)!=1)           // read body
        return -1;

	if(mPdib->biBitCount!=24)
		return -2;                                  // not 24bpp

    return 0;
}


//-------------------------------------------------------------------
// write bitmap file header
int
Cbmp::writeHeader(FILE* fp)
{
    if(fwrite(&mBmpFileHdr, sizeof(bmpFileHdr), 1, fp)!=1)
        return -1;

    return 0;
}


//-------------------------------------------------------------------
// write bitmap file body
int
Cbmp::writeDib(FILE* fp)
{
    if(fwrite(mPdib , mDibSize, 1, fp)!=1)          // write bitmap body
        return -1;

    return 0;
}


//-------------------------------------------------------------------
// set bitmap file header
void
Cbmp::setBmpInfoHdr(const int width, const int height)
{
    mPdib->biSize         =sizeof(bmpInfoHdr);
    mPdib->biWidth        =width;
    mPdib->biHeight       =height;
    mPdib->biPlanes       =1;
    mPdib->biBitCount     =24;                      // 24 bpp
    mPdib->biCompression  =0;
    mPdib->biSizeImage    =0;
    mPdib->biXPelsPerMeter=0;
    mPdib->biYPelsPerMeter=0;
    mPdib->biClrUsed      =0;
    mPdib->biClrImportant =0;
}


//-------------------------------------------------------------------
// set bitmap info header
//
// set bitmap file header
// set mAbsHeight
// set mPixelPitch
// set mRowPitch
//
void
Cbmp::setBmpFileHdr(const int width, const int height)
{
    mAbsHeight=height>0 ? height : -(height);       //abs

    mPixelPitch=3;                                  // 24 bpp

    mRowPitch=width*mPixelPitch;                    // to 4byte boundary
    if(mRowPitch%4)
        mRowPitch=mRowPitch+(4-(mRowPitch%4));

    mBmpFileHdr.bfType='B'+'M'*256;
    mBmpFileHdr.bfSize=(mRowPitch*mAbsHeight)+sizeof(bmpFileHdr)+sizeof(bmpInfoHdr);
    mBmpFileHdr.bfReserved1=0;
    mBmpFileHdr.bfReserved2=0;
    mBmpFileHdr.bfOffBits=sizeof(bmpFileHdr)+sizeof(bmpInfoHdr);
}


//================ ^^^^^^ private ^^^^^^ ====================================




//-------------------------------------------------------------------
// load bitmap image from file
void
Cbmp::loadFromFile(const char* bmpFName)
{
    FILE* fp;
    struct stat statbuf;                            // for SIDBA

    SAFE_FREE(mPdib);                               // delete image

    if ((fp = fopen(bmpFName, "rb")) == 0)          // open bitmap file
        throw "input file open failed.";

    if (stat(bmpFName, &statbuf) != 0)              // for SIDBA
        throw "function stat() failed.";            // for SIDBA

    if (readHeader(fp) != 0)                        // read file header
    {
        fclose(fp);
        throw "failed to read bitmap file header.";
    }

    //mDibSize=mBmpFileHdr.bfSize-sizeof(bmpFileHdr); // size of dib
    mDibSize = statbuf.st_size - sizeof(bmpFileHdr);    // for SIDBA
    mPdib = (bmpInfoHdr *)malloc(mDibSize);         // alloc dib memory

    if (readDib(fp) != 0)                           // read dib
    {
        SAFE_FREE(mPdib);
        fclose(fp);
        throw "failed to read bitmap file body.";
    }
    fclose(fp);                                     // close bitmap file

    mPbitmap = (unsigned char *)(mPdib)             // move pos. to body
        +mBmpFileHdr.bfOffBits - sizeof(bmpFileHdr);

    mPixelPitch = mPdib->biBitCount / 8;

    mRowPitch = (mPdib->biWidth*mPixelPitch);       // clac. row pitch by bytes
    if (mRowPitch % 4 != 0)
        mRowPitch += (4 - (mRowPitch % 4));

    mAbsHeight = mPdib->biHeight > 0 ? mPdib->biHeight : -(mPdib->biHeight);   //abs
    mImageSize = mRowPitch*mAbsHeight;
}


//-------------------------------------------------------------------
// get mem addr of specified scanrow#
unsigned char*
Cbmp::getScanRow(const int rowNo) const
{
    int absrowNo;

    if(mPdib==0)
        return 0;

    absrowNo=rowNo;
    if(mPdib->biHeight<0)
        absrowNo=mPdib->biHeight-rowNo-1;

    return (mPbitmap+(absrowNo*mRowPitch));
}


//-------------------------------------------------------------------
// save to bitmap file
void
Cbmp::saveToFile(const char* bmpFName)
{
    FILE* fp;

    if((fp=fopen(bmpFName, "wb"))!=0)               // open file
    {
        if(writeHeader(fp)==0)                      // write header
        {
            if(writeDib(fp)!=0)                     // write dib
                throw "failed to write dib.";
        }
        else
            throw "failed to write header.";
    }
    else
        throw "failed to open file.";

    fclose(fp);
}


//-------------------------------------------------------------------
// color to grayscale
void
Cbmp::getGSData(unsigned char* gs) const
{
    unsigned char* pRow = mPbitmap;
    unsigned char* pDest = gs;

    for (int y = 0; y < mAbsHeight; y++)
    {
        for (int x = 0; x < getWidth(); x++)
        {
            float m = (float)pRow[(x*mPixelPitch) + 0] * 0.114478f  // blue
                + (float)pRow[(x*mPixelPitch) + 1] * 0.586611f		// green
                + (float)pRow[(x*mPixelPitch) + 2] * 0.298912f;		// red

            *pDest = (unsigned char)m;								// gray scale
            pDest++;
        }
        pRow += mRowPitch;
    }
}


//-------------------------------------------------------------------
// create 24 bit DIB
int
Cbmp::create24Dib(const int width, const int height)
{
    setBmpFileHdr(width, height);

    SAFE_FREE(mPdib);                               // delete bmp
    mDibSize=mBmpFileHdr.bfSize-sizeof(bmpFileHdr); // size of dib
    mPdib=(bmpInfoHdr *)malloc(mDibSize);           // alloc dib memory

    setBmpInfoHdr(width, height);

    mPbitmap=(unsigned char *)(mPdib)               // move pos. to body
                            +mBmpFileHdr.bfOffBits
                                    -sizeof(bmpFileHdr);

    mImageSize=mRowPitch*mAbsHeight;

    memset(mPbitmap, 0xFF, mImageSize);             // init. image data

    return 0;
}






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
            throw "引数に<入力ファイル>を指定してください.\n";

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
