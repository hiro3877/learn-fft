#include <fftw3.h>
#include <complex.h>

int main(int argc, char *argv[])
{
  fftwf_complex *buf = (fftwf_complex *)fftw_malloc(10);

  buf[0][1] = 0.0f;
  buf[0][1] = 0.0f;

  fftw_free(buf);

  return 0;
}
