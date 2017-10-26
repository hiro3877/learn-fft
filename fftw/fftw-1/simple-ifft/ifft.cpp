#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <fftw3.h>
#include <complex.h>

//--------------------------------------------------------------------
//countLines
size_t
countComplexLines(const char* fname)
{
    FILE  *fp;
    double data[2];

    if ((fp = fopen(fname, "rt")) == NULL)
        throw "入力ファイルをオープンできません.";

    int count = 0;
    while (fscanf(fp, "%lf,%lf", &data[0], &data[1]) == 2)
        count++;

    fclose(fp);

    if (count <= 0)
        throw "入力ファイルの読み込み失敗.";

    return count;
}

//--------------------------------------------------------------------
//readComplexData
void
readComplexData(const char* fname, fftw_complex * buf, const size_t length)
{
    FILE *fp;
    double f[2];

    if ((fp = fopen(fname, "rt")) == NULL)
        throw "エラー：入力ファイルをオープンできません.";

    for (int i = 0; i < length; i++)
    {
        if (fscanf(fp, "%lf,%lf", &f[0], &f[1]) != 2)
            throw "エラー：入力ファイルの読み込み失敗.";

        buf[i][0] = f[0];
        buf[i][1] = f[1];
    }
    fclose(fp);
}

//--------------------------------------------------------------------
// main
int
main(int argc, char *argv[])
{
    fftw_plan plan;
    fftw_complex *ibuf, *obuf;

    try
    {
        if (argc != 2)
            throw "<入力ファイル名> を指定してください.";

        size_t length = countComplexLines(argv[1]);

        ibuf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*length);
        obuf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*length);

        readComplexData(argv[1], ibuf, length);

        //１次元のIDFTを実行,IFFT
        plan = fftw_plan_dft_1d((int)length, ibuf, obuf, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        //正規化
        for (int i = 0; i < length; i++)
        {
            obuf[i][0] /= (double)length;
            obuf[i][1] /= (double)length;
        }

        //出力
        for (size_t i = 0; i < length; i++)
        {
            fprintf(stdout, "%lf,%lf\n", obuf[i][0], obuf[i][1]);
        }

        fftw_free(ibuf);    // メモリ解放
        fftw_free(obuf);
    }
    catch (char *str)
    {
        fputs(str, stderr);
    }
    return 0;
}
