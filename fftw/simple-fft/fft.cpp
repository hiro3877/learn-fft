#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <fftw3.h>
#include <complex.h>

//--------------------------------------------------------------------
//countLines
size_t
countLines(const char* fname)
{
    FILE  *fp;
    double data;

    if ((fp = fopen(fname, "rt")) == NULL)
        throw "入力ファイルをオープンできません.";

    int count = 0;
    while (fscanf(fp, "%lf", &data) == 1)
        count++;

    fclose(fp);

    if (count <= 0)
        throw "入力ファイルの読み込み失敗.";

    return count;
}

//--------------------------------------------------------------------
//readData
void
readData(const char* fname, fftw_complex * buf, const size_t length)
{
    FILE *fp;
    double f;

    if ((fp = fopen(fname, "rt")) == NULL)
        throw "エラー：入力ファイルをオープンできません.";

    for (int i = 0; i < length; i++)
    {
        if (fscanf(fp, "%lf", &f) != 1)
            throw "エラー：入力ファイルの読み込み失敗.";

        buf[i][0] = f;
        buf[i][1] = 0.0;
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

        size_t length = countLines(argv[1]);

        ibuf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*length);
        obuf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*length);

        readData(argv[1], ibuf, length);

        //１次元のDFTを実行,FFT
        plan = fftw_plan_dft_1d((int)length, ibuf, obuf, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

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
