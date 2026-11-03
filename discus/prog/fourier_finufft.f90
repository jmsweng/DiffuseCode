module four_finufft
!-
!  FOURIER transformatons with finufft library
!          1/2/3 dimensional non-uniform-FFT
!          At the moment axes must be parallel to a*; b*; c*
!          If less than 1/cr_icc(:) grid points per reciprocal unt cell are
!          used, the coordinates are multiplied with scales(:) 
!          Locally an allocatable array fk is used, data are copied into 
!          the "diffuse" variable tcsf in sequence with top the fastes and
!          abscissa the slowest index.
!
!+
use finufft_mod
!
contains
!
!*******************************************************************************
!
subroutine four_strucf_1d(MAXATOMS, nat, inc, MAXQXY, xpos, occ, iscales, lscales, tcsf)
!-
!  Calculate a 1D-NUFFT via finufft
!+
!
use iso_c_binding
!
use errlist_mod
use precision_mod
use wink_mod
!
implicit none
!
integer                                   , intent(in)  :: MAXATOMS ! Dimension of xpos
integer                                   , intent(in)  :: nat      ! Actual atom number
integer                                   , intent(in)  :: inc      ! Number data points along axes
integer              , dimension(3)       , intent(in)  :: MAXQXY   ! Actual Reciprocal space dimensions
real(kind=PREC_DP)   , dimension(MAXATOMS), intent(inout)  :: xpos     ! Actual atom coordinates
real(kind=PREC_DP)                        , intent(in)  :: occ      ! Atom occupancies
integer                                   , intent(in)  :: iscales  ! Scaling if DELTA (hkl) /= 1/cr_icc
integer                                   , intent(in)  :: lscales  ! Scaling if scale > MAXSCALE
complex(kind=PREC_DP), dimension(MAXQXY(1),1,1), intent(out) :: tcsf     ! Resulting structure factor
!
integer                     :: ier   = 0  ! Return error flag from finufft
integer                     :: iflag = 1  ! Flag to FINUFFT sign of i
integer(kind=PREC_INT_LONG) :: m          ! Number of non-uniform points == number atoms
integer(kind=PREC_INT_LONG) :: n          ! Number of uniform points in Fourier space
real(kind=PREC_DP)                               :: tol ! Required tolerance
complex(kind=PREC_DP), dimension(:), allocatable :: cj  ! Real space amplitudes = 1*occ
complex(kind=PREC_DP), dimension(:), allocatable :: fk  ! Complex structure factor
!
type(finufft_opts) :: opts
!
iflag = 1
M   = nat
N   = (inc-1) * lscales          + 1
allocate(cj(M))
allocate(fk(N))
!
cj = cmplx(occ, 0.0D0, kind=PREC_DP)
tol = 1d-9
!
call finufft_default_opts(opts)                            ! Use default options 
call finufft1d1(M, xpos(1:nat), cj, iflag, tol, N, fk, opts, ier)   ! Do 1D NUFFT
!
if(ier==0) then
   tcsf(1:MAXQXY(1), 1, 1) = fk( 1: N:lscales)
endif
!
deallocate(cj)
deallocate(fk)
!
end subroutine four_strucf_1d
!
!******************************************************************************
!
subroutine four_strucf_2d(MAXATOMS, nat, inc, MAXQXY, xpos, ypos, occ, iscales, lscales, tcsf)
!-
!  Calculate a 2D-NUFFT via finufft
!+
!
use errlist_mod
use precision_mod
use wink_mod
!
implicit none
!
integer                                   , intent(in)  :: MAXATOMS ! Dimension of xpos
integer                                   , intent(in)  :: nat      ! Actual atom number
integer              , dimension(2)       , intent(in)  :: inc      ! Number data points along axes
integer              , dimension(3)       , intent(in)  :: MAXQXY   ! Total number in reciprocal space
real(kind=PREC_DP)   , dimension(MAXATOMS), intent(inout)  :: xpos     ! Actual atom coordinates
real(kind=PREC_DP)   , dimension(MAXATOMS), intent(inout)  :: ypos     ! Actual atom coordinates
real(kind=PREC_DP)                        , intent(in)  :: occ      ! Atom occupancies
integer              , dimension(2)       , intent(in)  :: iscales  ! Scaling if DELTA (hkl) /= 1/cr_icc
integer              , dimension(2)       , intent(in)  :: lscales  ! Scaling if scales > MAXSCALE
complex(kind=PREC_DP), dimension(MAXQXY(1), MAXQXY(2), 1)  , intent(out) :: tcsf     ! Resulting structure factor
!
integer                     :: ier   = 0  ! Return error flag from finufft
integer                     :: iflag = 1  ! Flag to FINUFFT sign of i
integer(kind=PREC_INT_LONG) :: m          ! Number of non-uniform points == number atoms
integer(kind=PREC_INT_LONG) :: ms, mt     ! Number of uniform points in Fourier space
integer(kind=PREC_INT_LONG) :: h, k, i    ! Loop indices             in Fourier space
integer(kind=PREC_INT_LONG) :: hh, kk     ! Loop indices             in Fourier space
integer(kind=PREC_INT_LONG) :: nj         ! Number of uniform points in Fourier space
real(kind=PREC_DP)                               :: tol ! Required tolerance
complex(kind=PREC_DP), dimension(:), allocatable :: cj  ! Real space amplitudes = 1*occ
complex(kind=PREC_DP), dimension(:), allocatable :: fk  ! Complex structure factor
complex(kind=PREC_DP), dimension(:,:), allocatable :: tt  ! Complex structure factor
!
type(finufft_opts) :: opts                  ! Use default options
!
iflag = 1
M   = nat
ms  = (inc(1)-1) * lscales(1) + 1
mt  = (inc(2)-1) * lscales(2) + 1
nj = ms*mt
!
allocate(cj(M))
allocate(fk(nj))
!
cj = cmplx(occ, 0.0D0, kind=PREC_DP)
tol = 1d-9
!
call finufft_default_opts(opts)
call finufft2d1(M, xpos(1:nat), ypos(1:nat), cj, iflag, tol, ms, mt, fk, opts, ier)
deallocate(cj)
i = 0
!
! Copy fk onto a 2D grid
!
allocate(tt(ms, mt))
do h=1, ms
   do k=1, mt
      tt(h,k) = fk((k-1)*ms+h)
   enddo
enddo
!
! Cut the 2D grid with steps of lscales(*)
!
hh = 0
do h=1, ms, lscales(1)
   hh = hh + 1
   kk = 0
   do k=1, mt, lscales(2)
      kk = kk + 1
      tcsf(hh,kk,1) = tt(h,k)
   enddo
enddo
!
deallocate(fk)
deallocate(tt)
!
end subroutine four_strucf_2d
!
!******************************************************************************
!
subroutine four_strucf_3d(MAXATOMS, nat, inc, MAXQXY, xpos, ypos, zpos, occ, iscales, lscales, tcsf)
!-
!  Calculate a 3D-NUFFT via finufft
!+
!
use errlist_mod
use precision_mod
use wink_mod
use support_mod
!
implicit none
!
integer                                   , intent(in)  :: MAXATOMS ! Dimension of xpos
integer                                   , intent(in)  :: nat      ! Actual atom number
integer              , dimension(3)       , intent(in)  :: inc      ! Number data points along exes
integer              , dimension(3)       , intent(in)  :: MAXQXY   ! Total data points reciprocal space
real(kind=PREC_DP)   , dimension(MAXATOMS), intent(inout)  :: xpos     ! Actual atom coordinates
real(kind=PREC_DP)   , dimension(MAXATOMS), intent(inout)  :: ypos     ! Actual atom coordinates
real(kind=PREC_DP)   , dimension(MAXATOMS), intent(inout)  :: zpos     ! Actual atom coordinates
real(kind=PREC_DP)                        , intent(in)  :: occ      ! Atom occupancies
integer              , dimension(3)       , intent(in)  :: iscales  ! Scaling if DELTA (hkl) /= 1/cr_icc
integer              , dimension(3)       , intent(in)  :: lscales  ! Scaling if scale > MAXSCALE
complex(kind=PREC_DP), dimension(MAXQXY(1),MAXQXY(2),MAXQXY(3))  , intent(out) :: tcsf     ! Resulting structure factor
!
integer                     :: ier   = 0  ! Return error flag from finufft
integer                     :: iflag = 1  ! Flag to FINUFFT sign of i
integer(kind=PREC_INT_LONG) :: m          ! Number of non-uniform points == number atoms
integer(kind=PREC_INT_LONG) :: ms, mt, mu ! Number of uniform points in Fourier space
integer(kind=PREC_INT_LONG) :: h, k, l, i ! Loop indices             in Fourier space
integer(kind=PREC_INT_LONG) :: hh,kk,ll   ! Loop indices             in Fourier space
integer(kind=PREC_INT_LONG) :: nj         ! Number of uniform points in Fourier space
real(kind=PREC_DP)                               :: tol ! Required tolerance
complex(kind=PREC_DP), dimension(:), allocatable :: cj  ! Real space amplitudes = 1*occ
complex(kind=PREC_DP), dimension(:), allocatable :: fk  ! Complex structure factor
complex(kind=PREC_DP), dimension(:,:,:), allocatable :: ttcsf  ! Complex structure factor
!
type(finufft_opts) :: opts                  ! Use default options
!
iflag = 1
M   = nat
ms  = inc(1)
mt  = inc(2)
mu  = inc(3)
ms  = (inc(1)-1) * lscales(1) + 1
mt  = (inc(2)-1) * lscales(2) + 1
mu  = (inc(3)-1) * lscales(3) + 1
nj = ms*mt*mu
allocate(cj(M))
allocate(fk(nj))
!
cj = cmplx(occ, 0.0D0, kind=PREC_DP)
tol = 1.0D-7
!
call finufft_default_opts(opts)
call finufft3d1(M, xpos(1:nat), ypos(1:nat), zpos(1:nat), cj, iflag, tol, ms, mt, mu, fk, opts, ier)
deallocate(cj)
!
! copy fk ont 3D grid
!
allocate(ttcsf(ms,mt,mu))
i = 0
do l=1, mu
   do k=1, mt
      do h=1,ms
         i = i+1
         ttcsf(h ,k ,l ) = fk(i)
      enddo
   enddo
enddo
!
! Cut 3D grid at steps lscales(*)
l = 0
do ll=1,mu,lscales(3)
   l = l + 1
   k = 0
   do kk=1,mt,lscales(2)
      k = k + 1
      h = 0
      do hh=1,ms,lscales(1)
         h = h + 1
         tcsf(h,k,l) = ttcsf(hh,kk,ll)
      enddo
   enddo
enddo
!
deallocate(fk)
deallocate(ttcsf)
!
!
end subroutine four_strucf_3d
!
!******************************************************************************
!
subroutine four_strucf_nufft(nat, pos, occ, eck1, vi, inc, tol, tcsf, ier)
!-
!  Structure factor of nat atoms on an arbitrary regular grid
!     h(t1,t2,t3) = eck1 + (t1-1)*vi(:,1) + (t2-1)*vi(:,2) + (t3-1)*vi(:,3)
!  via a type-1 NUFFT of dimension = number of axes with inc>1.
!
!  With hc = eck1 + sum_a floor(inc(a)/2)*vi(:,a) (the grid point FINUFFT
!  maps to mode 0) the phase splits into
!     2 pi h.x = 2 pi hc.x + sum_a k_a * (2 pi vi(:,a).x)
!  so the NUFFT coordinates are 2 pi vi(:,a).x wrapped into [-pi,pi] and the
!  strengths are occ*exp(2 pi i hc.x). No condition on step size, crystal
!  size or grid centre, increment vectors need not be parallel to a*,b*,c*.
!+
use precision_mod
use wink_mod
!
implicit none
!
integer                                  , intent(in)  :: nat      ! Number of atoms
real(kind=PREC_DP), dimension(3, nat)    , intent(in)  :: pos      ! Fractional coordinates
real(kind=PREC_DP)                       , intent(in)  :: occ      ! Occupancy
real(kind=PREC_DP), dimension(3)         , intent(in)  :: eck1     ! Lower left bottom corner
real(kind=PREC_DP), dimension(3,3)       , intent(in)  :: vi       ! Increment vectors (columns)
integer           , dimension(3)         , intent(in)  :: inc      ! Grid points along abs, ord, top
real(kind=PREC_DP)                       , intent(in)  :: tol      ! FINUFFT tolerance
complex(kind=PREC_DP), dimension(inc(1)*inc(2)*inc(3)), intent(out) :: tcsf  ! abscissa fastest
integer                                  , intent(out) :: ier      ! FINUFFT error flag
!
integer                     :: iflag      ! Sign of i in FINUFFT
integer                     :: nd         ! Number of non-flat axes
integer, dimension(3)       :: axes       ! The non-flat axes
integer(kind=PREC_INT_LONG) :: m          ! Number of atoms
integer(kind=PREC_INT_LONG), dimension(3) :: nm   ! Modes along non-flat axes
real(kind=PREC_DP), dimension(3)          :: hc   ! Grid point at FINUFFT mode 0
real(kind=PREC_DP), dimension(:,:), allocatable :: xc  ! NUFFT coordinates
complex(kind=PREC_DP), dimension(:), allocatable :: cj ! Strengths
integer :: a, d, k
real(kind=PREC_DP) :: arg
type(finufft_opts) :: opts
!
ier   = 0
iflag = 1
nd    = 0
axes  = 0
nm    = 1
hc    = eck1
do a = 1, 3
   hc = hc + real(inc(a)/2, kind=PREC_DP)*vi(:,a)
   if(inc(a) > 1) then
      nd       = nd + 1
      axes(nd) = a
      nm(nd)   = inc(a)
   endif
enddo
!
allocate(cj(nat))
allocate(xc(nat, max(nd,1)))
do k = 1, nat
   arg   = zpi*dot_product(hc, pos(:,k))
   arg   = arg - zpi*anint(arg/zpi)
   cj(k) = occ*cmplx(cos(arg), sin(arg), kind=PREC_DP)
   do d = 1, nd
      arg     = zpi*dot_product(vi(:,axes(d)), pos(:,k))
      xc(k,d) = arg - zpi*anint(arg/zpi)            ! exp(i k x), k integer: 2 pi periodic
   enddo
enddo
!
m = nat
call finufft_default_opts(opts)
select case(nd)
case(0)
   tcsf(1) = sum(cj)
case(1)
   call finufft1d1(m, xc(:,1),                   cj, iflag, tol, nm(1),               tcsf, opts, ier)
case(2)
   call finufft2d1(m, xc(:,1), xc(:,2),          cj, iflag, tol, nm(1), nm(2),        tcsf, opts, ier)
case(3)
   call finufft3d1(m, xc(:,1), xc(:,2), xc(:,3), cj, iflag, tol, nm(1), nm(2), nm(3), tcsf, opts, ier)
end select
!  FINUFFT returns modes with the first axis fastest; axes of length 1 do not
!  change the linear order, so tcsf already has the layout of tcsf(inc(1),inc(2),inc(3)).
!
deallocate(cj, xc)
!
end subroutine four_strucf_nufft
!
!******************************************************************************
!
end module four_finufft
