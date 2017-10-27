//====================================================================
// fft
//
// program�� ����.txt > foo.csv
//
// (c)Copyright Spacesoft corp., 2017 All rights reserved.
//                                    Kitayama, Hiroyuki
//====================================================================
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <fftw3.h>
#include <complex.h>

//--------------------------------------------------------------------
//countLines
size_t countLines(const char* fname)
{
    FILE  *fp;
    double data;

    if ((fp = fopen(fname, "rt")) == NULL)
        throw "cant open.";

    int count = 0;
    while (fscanf(fp, "%lf", &data) == 1)
        count++;

    fclose(fp);

    if (count <= 0)
        throw "cant read.";

    return count;
}

//--------------------------------------------------------------------
//readData
void readData(const char* fname, fftw_complex* buf, const size_t length)
{
    FILE *fp;
    double f;

    if ((fp = fopen(fname, "rt")) == NULL)
        throw "cant open.";

    for (int i = 0; i < length; i++)
    {
        if (fscanf(fp, "%lf", &f) != 1)
            throw "cant read.";

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
            throw "please input file name.";

        size_t length = countLines(argv[1]);


        readData(argv[1], ibuf, length);

        //1 dimension FFT
        plan = fftw_plan_dft_1d((int)length, ibuf, obuf, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        //output
        for (size_t i = 0; i < length; i++)
        {
            fprintf(stdout, "%lf,%lf\n", obuf[i][0], obuf[i][1]);
        }


        fftw_free(ibuf);    // memory free
        fftw_free(obuf);
    }

    catch (char *str)
    {
        fputs(str, stderr);
    }

    return 0;
}
