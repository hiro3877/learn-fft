1d : fft-fortran-1d-t1.c
	gcc fft-fortran-1d-t1.c -lm -O3 -o fft-1d.out
 
2d : fft-fortran-2d-t.c
	gcc fft-fortran-2d-t.c -lm -O3 -o fft-2d.out

2d-1 : fft-2d.c
	gcc fft-2d.c -lm -O3 -o fft-2d1.out

1d-icc : fft-fortran-1d-t1.c
	icc fft-fortran-1d-t1.c -lm -O3 -o fft-1d-icc.out

2d-icc : fft-fortran-2d-t.c
	icc fft-fortran-2d-t.c -lm -O3 -o fft-2d-icc.out

1d-iccp : fft-fortran-1d-t1.c
	icc fft-fortran-1d-t1.c -lm -parallel -O3 -o fft-1d-iccp.out

2d-iccp : fft-fortran-2d-t.c
	icc fft-fortran-2d-t.c -lm -parallel -O3 -o fft-2d-iccp.out

