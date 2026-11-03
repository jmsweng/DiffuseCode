! Validate four_strucf_gemm against the repo's table kernel and an exact sum,
! using DISCUS-style lmn / off_shift start-phase bookkeeping and a skewed grid.
program test_gemm
use omp_lib
use four_strucf_gemm_mod
implicit none
integer, parameter :: dp = kind(0.d0), il = selected_int_kind(18)
integer(il), parameter :: I2PI = 2_il**16, MASK = I2PI-1, SHIFT = -6
real(dp), parameter :: zpi = 8.d0*atan(1.d0)
integer :: nxat, num(3), lmn(6), k, s, j, i, h
real(dp) :: uin(3), vin(3), win(3), off(3,3), r(3), q(3), t0, t1, e1, e2, big, ph
real(dp), allocatable :: xat(:,:)
complex(dp), allocatable :: cex(:), ft(:,:,:), fg(:,:,:)
complex(dp) :: ref
nxat = 20000; num = [81, 71, 61]
allocate(xat(nxat,3)); call random_number(xat); xat = 40.d0*xat - 20.d0
uin = [0.05d0, 0.0d0, 0.0d0]; vin = [0.025d0, 0.05d0, 0.0d0]; win = [0.0d0, 0.0d0, 0.05d0]  ! skewed grid
lmn = [-40, -20, -30, 1, 1, 1]
off = reshape([0.013d0,0.0d0,0.0d0, 0.0d0,0.021d0,0.0d0, 0.0d0,0.0d0,0.007d0],[3,3])
allocate(cex(0:MASK)); do k=0,MASK; cex(k)=cmplx(cos(zpi*k/I2PI),sin(zpi*k/I2PI),dp); end do
allocate(ft(num(1),num(2),num(3)), fg(num(1),num(2),num(3)))
t0=omp_get_wtime(); call table(ft); t1=omp_get_wtime(); write(*,'(a,f8.3,a)') 'repo table kernel (omp): ', t1-t0,' s'
t0=omp_get_wtime(); call four_strucf_gemm(nxat, xat, uin, vin, win, off, lmn, num, num, fg); t1=omp_get_wtime()
write(*,'(a,f8.3,a)') 'four_strucf_gemm       : ', t1-t0,' s'
e1=0; e2=0; big=0
do s=1,200
  call random_number(r); j=1+int(r(1)*num(1)); i=1+int(r(2)*num(2)); h=1+int(r(3)*num(3))
  ref=0
  do k=1,nxat
    q = matmul(xat(k,:),off)
    ph = zpi*( lmn(1)*dot_product(uin,xat(k,:)) + lmn(2)*dot_product(vin,xat(k,:)) + lmn(3)*dot_product(win,xat(k,:)) &
             + lmn(4)*q(1) + lmn(5)*q(2) + lmn(6)*q(3) &
             + (j-1)*dot_product(uin,xat(k,:)) + (i-1)*dot_product(vin,xat(k,:)) + (h-1)*dot_product(win,xat(k,:)) )
    ref = ref + cmplx(cos(ph),sin(ph),dp)
  end do
  big=max(big,abs(ref)); e1=max(e1,abs(ft(j,i,h)-ref)); e2=max(e2,abs(fg(j,i,h)-ref))
end do
write(*,'(a,es10.2,a,es10.2,a,es10.2)') 'max|F| ', big, '   err table ', e1, '   err gemm ', e2
contains
subroutine table(f)           ! four_strucf_omp inner loops, verbatim logic
  complex(dp), intent(out) :: f(num(1),num(2),num(3))
  integer(il) :: iu(nxat),iv(nxat),iw(nxat),ju(nxat),jv(nxat),jw(nxat), iarg, iarg0
  real(dp) :: xu,xv,xw,o(3)
  complex(dp), allocatable :: p(:,:)
  integer :: j,i,h,k
  do k=1,nxat
    xu=dot_product(uin,xat(k,:)); xv=dot_product(vin,xat(k,:)); xw=dot_product(win,xat(k,:)); o=matmul(xat(k,:),off)
    iu(k)=nint(64*I2PI*(xu-int(xu)),il); iv(k)=nint(64*I2PI*(xv-int(xv)),il); iw(k)=nint(64*I2PI*(xw-int(xw)),il)
    ju(k)=nint(64*I2PI*(o(1)-int(o(1))),il); jv(k)=nint(64*I2PI*(o(2)-int(o(2))),il); jw(k)=nint(64*I2PI*(o(3)-int(o(3))),il)
  end do
  f=0
  !$omp parallel private(iarg,iarg0,p,i,h,k)
  allocate(p(num(2),num(3)))
  !$omp do schedule(static)
  do j=1,num(1)
    p=0
    do k=1,nxat
      iarg0 = lmn(1)*iu(k)+lmn(2)*iv(k)+lmn(3)*iw(k)+lmn(4)*ju(k)+lmn(5)*jv(k)+lmn(6)*jw(k)+(j-1)*iu(k)
      do i=1,num(2)
        iarg = iarg0 + iv(k)*(i-1)
        do h=1,num(3); p(i,h)=p(i,h)+cex(iand(ishft(iarg,SHIFT),MASK)); iarg=iarg+iw(k); end do
      end do
    end do
    f(j,:,:)=p
  end do
  !$omp end do
  !$omp end parallel
end subroutine
end program
