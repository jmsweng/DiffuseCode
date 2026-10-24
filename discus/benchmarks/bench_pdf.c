// 3D-PDF FFT timing: what out_prep_3dpdf_3d does now (complex-to-complex,
// FFTW_ESTIMATE, 1 thread) vs. real-to-complex and multithreaded FFTW.
// gcc -O2 bench_pdf.c -lfftw3_omp -lfftw3 -fopenmp -lm
#include <fftw3.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static double now(void){ return omp_get_wtime(); }

int main(int argc, char **argv){
  int n = argc > 1 ? atoi(argv[1]) : 201;
  int nt = omp_get_max_threads();
  size_t N = (size_t)n*n*n, Nc = (size_t)n*n*(n/2+1);
  double *dsi = fftw_malloc(sizeof(double)*N);
  for(size_t i=0;i<N;i++) dsi[i] = (double)rand()/RAND_MAX;
  fftw_complex *in = fftw_malloc(sizeof(fftw_complex)*N), *out = fftw_malloc(sizeof(fftw_complex)*N);
  fftw_complex *rc = fftw_malloc(sizeof(fftw_complex)*Nc);
  double t0, t1, tp;
  printf("grid %d^3, %d threads available\n", n, nt);

  // current: copy real->complex, c2c, 1 thread, ESTIMATE (plan+exec timed like DISCUS)
  t0 = now();
  for(size_t i=0;i<N;i++){ in[i][0]=dsi[i]; in[i][1]=0; }
  fftw_plan p = fftw_plan_dft_3d(n,n,n,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_execute(p); fftw_destroy_plan(p);
  t1 = now(); printf("%-45s %8.3f s  mem %.2f GB\n","c2c, 1 thread, ESTIMATE (as in repo)", t1-t0, 2.0*16*N/1e9);
  double ref = out[1+n*(2+n*3)][0];

  // r2c, 1 thread
  t0 = now();
  p = fftw_plan_dft_r2c_3d(n,n,n,dsi,rc,FFTW_ESTIMATE);
  fftw_execute(p); fftw_destroy_plan(p);
  t1 = now(); printf("%-45s %8.3f s  mem %.2f GB\n","r2c, 1 thread, ESTIMATE", t1-t0, 16.0*Nc/1e9);
  // r2c output index (z fastest, half-length): element [3][2][1] in C order of (x,y,z) = out[1+n*(2+n*3)] in Fortran sense
  double chk = rc[3*(size_t)n*(n/2+1) + 2*(n/2+1) + 1][0];

  fftw_init_threads(); fftw_plan_with_nthreads(nt);
  t0 = now();
  for(size_t i=0;i<N;i++){ in[i][0]=dsi[i]; in[i][1]=0; }
  p = fftw_plan_dft_3d(n,n,n,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_execute(p); fftw_destroy_plan(p);
  t1 = now(); printf("%-45s %8.3f s\n","c2c, all threads", t1-t0);

  t0 = now();
  p = fftw_plan_dft_r2c_3d(n,n,n,dsi,rc,FFTW_ESTIMATE);
  fftw_execute(p); fftw_destroy_plan(p);
  t1 = now(); printf("%-45s %8.3f s\n","r2c, all threads", t1-t0);
  printf("check same coefficient: c2c %.10g  r2c %.10g\n", ref, chk);
  return 0;
}
