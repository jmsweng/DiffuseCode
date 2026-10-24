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
!  Threading comes from the BLAS library (OpenBLAS/MKL/BLIS).
!
!  Call it where four_strucf_omp now does its two loops, i.e.
!     call four_strucf_gemm(nxat, xat, uin, vin, win, off_shift, lmn, fnum, num, tcsf)
!     call tcsf_form_generic(lnufft, ldiscamb, lform, is_anis, iscat, isym, ientry, fnum, tcsf)
!  Needs a BLAS at link time (-lopenblas, or MKL).
!-
implicit none
private
public :: four_strucf_gemm
integer, parameter :: dp = kind(0.d0)
integer, parameter :: NBLOCK = 2048            ! atoms per GEMM; ~ 2048*num*16 bytes per work array
real(dp), parameter :: ZPI = 8.0_dp*atan(1.0_dp)
!
contains
!
subroutine four_strucf_gemm(nxat, xat, uin, vin, win, off_shift, lmn, fnum, num, tcsf)
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
!
n1 = fnum(1); n2 = fnum(2); n3 = fnum(3)
tcsf(1:n1,1:n2,1:n3) = (0.0_dp,0.0_dp)
if(nxat < 1) return
allocate(U(NBLOCK,n1), V(NBLOCK,n2), A(NBLOCK,n1), w(NBLOCK), xw(NBLOCK))
!
do kb = 1, nxat, NBLOCK
   nb = min(NBLOCK, nxat-kb+1)
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
enddo
deallocate(U, V, A, w, xw)
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
