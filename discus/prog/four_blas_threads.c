/*
 * discus_blas_set_threads(n): set the number of threads the BLAS library uses
 * for its own calls, return the previous setting (or -1 if the BLAS library
 * offers no such control or it cannot be detected).
 *
 * Used by four_strucf_gemm: when the Fourier layers are distributed over
 * OpenMP threads, each thread calls ZGEMM and the BLAS must not start its own
 * thread pool as well (oversubscription; with a pthreads OpenBLAS this made
 * 16 threads ~16x slower than 1 thread).
 *
 * The OpenBLAS / MKL entry points are referenced as weak symbols, so this file
 * links against any BLAS. Weak references are only used on Linux/ELF; on other
 * platforms the function is a no-op and the user controls BLAS threads with
 * the environment (OPENBLAS_NUM_THREADS=1, MKL_NUM_THREADS=1, VECLIB_MAXIMUM_THREADS=1).
 */
#if defined(__linux__) && (defined(__GNUC__) || defined(__clang__))
#define DISCUS_WEAK_BLAS 1
extern int  openblas_get_num_threads(void)  __attribute__((weak));
extern void openblas_set_num_threads(int)   __attribute__((weak));
extern int  mkl_set_num_threads_local(int)  __attribute__((weak));
#endif

int discus_blas_set_threads(int n)
{
#ifdef DISCUS_WEAK_BLAS
    if (openblas_set_num_threads && openblas_get_num_threads) {
        int old = openblas_get_num_threads();
        openblas_set_num_threads(n);
        return old;
    }
    if (mkl_set_num_threads_local) {
        return mkl_set_num_threads_local(n);   /* previous thread-local value, 0 = global */
    }
#endif
    (void)n;
    return -1;
}
