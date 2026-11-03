module four_strucf_gemm_mod
!+
!  Drop-in alternative for the inner Fourier loop of four_strucf_omp /
!  four_strucf_serial (discus/prog/four_strucf_OMP_mod.f90).
!
!  Uses the separability of the phase on a regular grid
!     exp(2pi i [p0 + (j-1)u.x + (i-1)v.x + (h-1)w.x])
!   = U(k,j) * V(k,i) * W_h(k)
!  so that every layer h of the structure factor is one complex matrix product
!     tcsf(:,:,h) += (U o W_h)^T * V        (ZGEMM, atoms are the summed index)
!  The phases are evaluated exactly with sin/cos (no 2**16 lookup table), the
!  start phase p0 is built from lmn/off_shift exactly as iarg0 is in DISCUS.
!
!  Threading: for volumes (several layers) and more than one OpenMP thread the
!  layers are distributed over the OpenMP threads, each thread calls ZGEMM for
!  its own layers, and the BLAS library is switched to one thread meanwhile
!  (see four_blas_threads.c). Many small GEMMs threaded inside BLAS were much
!  slower than one thread. For a single layer the one GEMM keeps the BLAS
!  library's own threading.
!
!  Used by four_strucf (four_strucf_OMP_mod.f90) when DISCUS is configured
!  with -DDISCUS_BLAS=ON and the grid has more than one point along both
!  abscissa and ordinate. Set the environment variable
!  DISCUS_FOURIER_KERNEL=table to fall back to the lookup-table kernel.
!-
use iso_c_binding, only: c_int
implicit none
private
public :: four_strucf_gemm
integer, parameter :: dp = kind(0.d0)
integer, parameter :: NBLOCK = 2048            ! atoms per GEMM; ~ 2048*num*16 bytes per work array
real(dp), parameter :: ZPI = 8.0_dp*atan(1.0_dp)
!
interface
   integer(c_int) function discus_blas_set_threads(n) bind(C, name='discus_blas_set_threads')
      import :: c_int
      integer(c_int), value :: n
   end function discus_blas_set_threads
end interface
!
contains
!
subroutine four_strucf_gemm(nxat, xat, uin, vin, win, off_shift, lmn, fnum, num, tcsf)
!$ use omp_lib
integer                  , intent(in)    :: nxat           ! number of atoms in xat
real(dp), dimension(:,:) , intent(in)    :: xat            ! xat(k,1:3) fractional coordinates
real(dp), dimension(3)   , intent(in)    :: uin, vin, win  ! grid increments (abscissa, ordinate, top)
real(dp), dimension(3,3) , intent(in)    :: off_shift      ! as in diffuse_mod
integer , dimension(6)   , intent(in)    :: lmn            ! as in diffuse_mod
integer , dimension(3)   , intent(in)    :: fnum           ! points to compute (Friedel-reduced)
integer , dimension(3)   , intent(in)    :: num            ! leading dimensions of tcsf
complex(dp), dimension(num(1),num(2),num(3)), intent(inout) :: tcsf
!
complex(dp), parameter :: ONE = (1.0_dp,0.0_dp)
complex(dp), allocatable :: U(:,:), V(:,:), A(:,:), w(:)
real(dp), allocatable    :: xw(:)
real(dp) :: xu, xv, p0, x(3), o(3)
integer  :: kb, nb, kk, k, j, i, h, n1, n2, n3
integer  :: nthreads
logical  :: layer_par                       ! distribute layers over OpenMP threads
integer(c_int) :: blas_old
!
n1 = fnum(1); n2 = fnum(2); n3 = fnum(3)
tcsf = (0.0_dp,0.0_dp)     ! whole array: four_run adds tcsf(1:num) to csf, Friedel half stays zero
if(nxat < 1) return
!
nthreads = 1
!$ nthreads = omp_get_max_threads()
layer_par = nthreads > 1 .and. n3 > 1
blas_old  = -1
if(layer_par) blas_old = discus_blas_set_threads(1_c_int)
!
allocate(U(NBLOCK,n1), V(NBLOCK,n2), xw(NBLOCK))
!
do kb = 1, nxat, NBLOCK
   nb = min(NBLOCK, nxat-kb+1)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(k, x, o, p0, xu, xv, j, i)
   do kk = 1, nb
      k  = kb + kk - 1
      x  = [dot_product(uin, xat(k,1:3)), dot_product(vin, xat(k,1:3)), dot_product(win, xat(k,1:3))]
      o  = matmul(xat(k,1:3), off_shift)                  ! oincu, oincv, oincw
      x  = x - aint(x)
      o  = o - aint(o)
      p0 = sum(lmn(1:3)*x) + sum(lmn(4:6)*o)              ! same as iarg0 in DISCUS (in cycles)
      xu = x(1); xv = x(2); xw(kk) = x(3)
      do j = 1, n1
         U(kk,j) = cis(p0 + real(j-1,dp)*xu)
      enddo
      do i = 1, n2
         V(kk,i) = cis(real(i-1,dp)*xv)
      enddo
   enddo
!$OMP END PARALLEL DO
!
!  Layers are independent: each thread accumulates into its own tcsf(:,:,h)
!$OMP PARALLEL IF(layer_par) PRIVATE(A, w, h, kk, j)
   allocate(A(NBLOCK,n1), w(NBLOCK))
!$OMP DO SCHEDULE(STATIC)
   do h = 1, n3
      do kk = 1, nb
         w(kk) = cis(real(h-1,dp)*xw(kk))
      enddo
      do j = 1, n1
         A(1:nb,j) = U(1:nb,j)*w(1:nb)
      enddo
!     tcsf(1:n1,1:n2,h) += A^T * V       (leading dimension num(1) writes straight into tcsf)
      call zgemm('T', 'N', n1, n2, nb, ONE, A, NBLOCK, V, NBLOCK, ONE, tcsf(1,1,h), num(1))
   enddo
!$OMP END DO
   deallocate(A, w)
!$OMP END PARALLEL
enddo
deallocate(U, V, xw)
!
if(blas_old >= 0) blas_old = discus_blas_set_threads(blas_old)
!
contains
   pure complex(dp) function cis(cycles)
   real(dp), intent(in) :: cycles
   real(dp) :: t
   t   = ZPI*(cycles - anint(cycles))
   cis = cmplx(cos(t), sin(t), kind=dp)
   end function cis
end subroutine four_strucf_gemm
!
end module four_strucf_gemm_mod
