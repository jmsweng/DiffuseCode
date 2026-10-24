! Standalone benchmark replicating the DISCUS "turbo" Fourier kernel
! (four_strucf_serial / four_strucf_omp) and candidate replacements.
! Build:  gfortran -O3 -ffast-math -fopenmp [-march=native] bench_strucf.f90 -lopenblas
program bench_strucf
use omp_lib
implicit none
integer, parameter :: dp = kind(0.d0), il = selected_int_kind(18)
integer(il), parameter :: I2PI = 2_il**16, MASK = I2PI-1, SHIFT = -6
real(dp), parameter :: zpi = 8.d0*atan(1.d0)
integer :: nc, nxat, n1, n2, n3, k, ix, iy, iz, ii, it
real(dp), allocatable :: xat(:,:)
complex(dp), allocatable :: cex(:), f_ser(:,:,:), f_omp(:,:,:), f_fix(:,:,:), f_rec(:,:,:), f_gem(:,:,:)
real(dp) :: eck(3), uin(3), vin(3), win(3), r(3), t0, t1
character(len=32) :: arg
integer :: nthr

nc = 30; n1 = 61
if (command_argument_count() >= 1) then; call get_command_argument(1,arg); read(arg,*) nc; end if
if (command_argument_count() >= 2) then; call get_command_argument(2,arg); read(arg,*) n1; end if
n2 = n1; n3 = n1
nxat = nc**3
eck = -3.d0; uin = [6.d0/(n1-1),0.d0,0.d0]; vin = [0.d0,6.d0/(n2-1),0.d0]; win = [0.d0,0.d0,6.d0/(n3-1)]
allocate(xat(nxat,3))
call random_seed()
k = 0
do iz=0,nc-1; do iy=0,nc-1; do ix=0,nc-1
  k = k+1; call random_number(r)
  xat(k,:) = [real(ix,dp),real(iy,dp),real(iz,dp)] + 0.1d0*(r-0.5d0) - nc/2
end do; end do; end do

allocate(cex(0:MASK))
do ii=0,MASK; cex(ii) = cmplx(cos(zpi*ii/I2PI), sin(zpi*ii/I2PI), dp); end do
allocate(f_ser(n1,n2,n3), f_omp(n1,n2,n3), f_fix(n1,n2,n3), f_rec(n1,n2,n3), f_gem(n1,n2,n3))
nthr = omp_get_max_threads()
write(*,'(a,i9,a,3i5,a,i3,a,es10.3)') 'atoms', nxat, '  grid', n1,n2,n3, '  threads', nthr, &
      '  atom*points', real(nxat,dp)*n1*n2*n3

t0 = omp_get_wtime(); call k_serial_orig(f_ser); t1 = omp_get_wtime()
write(*,'(a40,f9.3,a)') 'A  four_strucf_serial (as in repo):', t1-t0, ' s'
t0 = omp_get_wtime(); call k_omp_orig(f_omp); t1 = omp_get_wtime()
write(*,'(a40,f9.3,a)') 'B  four_strucf_omp    (as in repo):', t1-t0, ' s'
t0 = omp_get_wtime(); call k_omp_fixed(f_fix); t1 = omp_get_wtime()
write(*,'(a40,f9.3,a)') 'C  omp, contiguous tcsfp(h,i):', t1-t0, ' s'
t0 = omp_get_wtime(); call k_recur(f_rec); t1 = omp_get_wtime()
write(*,'(a40,f9.3,a)') 'D  phase recurrence, SIMD over atoms:', t1-t0, ' s'
t0 = omp_get_wtime(); call k_gemm(f_gem); t1 = omp_get_wtime()
write(*,'(a40,f9.3,a)') 'E  separable phases + ZGEMM:', t1-t0, ' s'
call check()

contains

subroutine incs(iincu,iincv,iincw,iarg0)
  integer(il), intent(out) :: iincu(nxat), iincv(nxat), iincw(nxat), iarg0(nxat)
  real(dp) :: xu, xv, xw, x0
  integer :: k
  do k=1,nxat
    xu = dot_product(uin,xat(k,:)); xv = dot_product(vin,xat(k,:)); xw = dot_product(win,xat(k,:))
    x0 = dot_product(eck,xat(k,:))
    iincu(k) = nint(64*I2PI*(xu-int(xu)), il); iincv(k) = nint(64*I2PI*(xv-int(xv)), il)
    iincw(k) = nint(64*I2PI*(xw-int(xw)), il); iarg0(k) = nint(64*I2PI*(x0-int(x0)), il)
  end do
end subroutine

subroutine k_serial_orig(f)        ! = four_strucf_serial
  complex(dp), intent(out) :: f(n1,n2,n3)
  integer(il) :: iincu(nxat), iincv(nxat), iincw(nxat), ia0(nxat), iarg, iarg0
  integer :: j,i,h,k
  call incs(iincu,iincv,iincw,ia0)
  f = (0.d0,0.d0)
  do j=1,n1
    do k=1,nxat
      iarg0 = ia0(k) + (j-1)*iincu(k)
      do i=1,n2
        iarg = iarg0 + iincv(k)*(i-1)
        do h=1,n3
          f(j,i,h) = f(j,i,h) + cex(iand(ishft(iarg,SHIFT),MASK))
          iarg = iarg + iincw(k)
        end do
      end do
    end do
  end do
end subroutine

subroutine k_omp_orig(f)           ! = four_strucf_omp
  complex(dp), intent(out) :: f(n1,n2,n3)
  integer(il) :: iincu(nxat), iincv(nxat), iincw(nxat), ia0(nxat), iarg, iarg0
  complex(dp), allocatable :: p(:,:)
  integer :: j,i,h,k
  call incs(iincu,iincv,iincw,ia0)
  f = (0.d0,0.d0)
  !$omp parallel private(iarg,iarg0,p,i,h,k)
  allocate(p(n2,n3))
  !$omp do schedule(static)
  do j=1,n1
    p = (0.d0,0.d0)
    do k=1,nxat
      iarg0 = ia0(k) + (j-1)*iincu(k)
      do i=1,n2
        iarg = iarg0 + iincv(k)*(i-1)
        do h=1,n3
          p(i,h) = p(i,h) + cex(iand(ishft(iarg,SHIFT),MASK))
          iarg = iarg + iincw(k)
        end do
      end do
    end do
    f(j,:,:) = f(j,:,:) + p
  end do
  !$omp end do
  deallocate(p)
  !$omp end parallel
end subroutine

subroutine k_omp_fixed(f)          ! same, but scratch is (h,i): unit stride in the hot loop
  complex(dp), intent(out) :: f(n1,n2,n3)
  integer(il) :: iincu(nxat), iincv(nxat), iincw(nxat), ia0(nxat), iarg, iarg0
  complex(dp), allocatable :: p(:,:)
  integer :: j,i,h,k
  call incs(iincu,iincv,iincw,ia0)
  !$omp parallel private(iarg,iarg0,p,i,h,k)
  allocate(p(n3,n2))
  !$omp do schedule(static)
  do j=1,n1
    p = (0.d0,0.d0)
    do k=1,nxat
      iarg0 = ia0(k) + (j-1)*iincu(k)
      do i=1,n2
        iarg = iarg0 + iincv(k)*(i-1)
        do h=1,n3
          p(h,i) = p(h,i) + cex(iand(ishft(iarg,SHIFT),MASK))
          iarg = iarg + iincw(k)
        end do
      end do
    end do
    do h=1,n3; do i=1,n2; f(j,i,h) = p(h,i); end do; end do
  end do
  !$omp end do
  deallocate(p)
  !$omp end parallel
end subroutine

subroutine k_recur(f)
  ! No lookup table: per atom, exp(i*phase) is advanced by complex multiplication
  ! with the per-step rotation exp(2*pi*i*w.x). The innermost loop runs over a
  ! block of atoms (independent recurrences) -> vectorises; accumulators are
  ! reduced per grid point. Rows are restarted from an exact sincos every row,
  ! so error never accumulates over more than n3 steps.
  complex(dp), intent(out) :: f(n1,n2,n3)
  integer, parameter :: NB = 256
  real(dp), allocatable :: wr(:),wi(:), zr(:),zi(:), rr(:),ri(:), vr(:),vi(:)
  real(dp) :: sr, si, tr, ph
  integer :: j,i,h,k,kb,nb_,kk
  complex(dp), allocatable :: p(:,:)
  !$omp parallel private(wr,wi,zr,zi,rr,ri,vr,vi,sr,si,tr,ph,i,h,k,kb,nb_,kk,p)
  allocate(wr(NB),wi(NB),zr(NB),zi(NB),rr(NB),ri(NB),vr(NB),vi(NB),p(n3,n2))
  !$omp do schedule(dynamic)
  do j=1,n1
    p = (0.d0,0.d0)
    do kb=1,nxat,NB
      nb_ = min(NB,nxat-kb+1)
      do kk=1,nb_
        k = kb+kk-1
        ph = zpi*(dot_product(eck,xat(k,:)) + (j-1)*dot_product(uin,xat(k,:)))
        rr(kk) = cos(ph); ri(kk) = sin(ph)            ! row start
        ph = zpi*dot_product(vin,xat(k,:)); vr(kk) = cos(ph); vi(kk) = sin(ph)
        ph = zpi*dot_product(win,xat(k,:)); wr(kk) = cos(ph); wi(kk) = sin(ph)
      end do
      do i=1,n2
        zr(1:nb_) = rr(1:nb_); zi(1:nb_) = ri(1:nb_)
        do h=1,n3
          sr = 0.d0; si = 0.d0
          !$omp simd reduction(+:sr,si) private(tr)
          do kk=1,nb_
            sr = sr + zr(kk); si = si + zi(kk)
            tr = zr(kk)*wr(kk) - zi(kk)*wi(kk)
            zi(kk) = zr(kk)*wi(kk) + zi(kk)*wr(kk)
            zr(kk) = tr
          end do
          p(h,i) = p(h,i) + cmplx(sr,si,dp)
        end do
        !$omp simd private(tr)
        do kk=1,nb_                                   ! advance row start by v
          tr = rr(kk)*vr(kk) - ri(kk)*vi(kk)
          ri(kk) = rr(kk)*vi(kk) + ri(kk)*vr(kk)
          rr(kk) = tr
        end do
      end do
    end do
    do h=1,n3; do i=1,n2; f(j,i,h) = p(h,i); end do; end do
  end do
  !$omp end do
  deallocate(wr,wi,zr,zi,rr,ri,vr,vi,p)
  !$omp end parallel
end subroutine

subroutine k_gemm(f)
  ! exp(2pi i (e0 + j u + i v + h w).x) = A_j(k) * B_i(k) * C_h(k)
  ! => F(:,:,j) = (A_j o B)^T * C  : one complex GEMM per j-plane (atoms blocked).
  ! BLAS runs near peak; phases computed exactly (sincos), no table quantisation.
  complex(dp), intent(out) :: f(n1,n2,n3)
  integer, parameter :: NB = 4096
  complex(dp), allocatable :: B(:,:), C(:,:), Bj(:,:), P(:,:), a0(:), au(:)
  integer :: j, i, h, k, kb, nb_, kk
  real(dp) :: ph
  complex(dp), parameter :: one=(1.d0,0.d0), zero=(0.d0,0.d0)
  allocate(B(NB,n2), C(NB,n3), Bj(NB,n2), P(n2,n3), a0(NB), au(NB))
  f = zero
  do kb=1,nxat,NB
    nb_ = min(NB,nxat-kb+1)
    do kk=1,nb_
      k = kb+kk-1
      ph = zpi*dot_product(eck,xat(k,:)); a0(kk) = cmplx(cos(ph),sin(ph),dp)
      ph = zpi*dot_product(uin,xat(k,:)); au(kk) = cmplx(cos(ph),sin(ph),dp)
      do i=1,n2; ph = zpi*(i-1)*dot_product(vin,xat(k,:)); B(kk,i) = cmplx(cos(ph),sin(ph),dp); end do
      do h=1,n3; ph = zpi*(h-1)*dot_product(win,xat(k,:)); C(kk,h) = cmplx(cos(ph),sin(ph),dp); end do
    end do
    do j=1,n1
      do i=1,n2; Bj(1:nb_,i) = B(1:nb_,i)*a0(1:nb_); end do
      call zgemm('T','N', n2, n3, nb_, one, Bj, NB, C, NB, zero, P, n2)
      f(j,:,:) = f(j,:,:) + P
      a0(1:nb_) = a0(1:nb_)*au(1:nb_)
    end do
  end do
end subroutine

subroutine check()
  ! exact reference on random grid points
  integer :: s, j, i, h, k
  real(dp) :: q(3), r(3), ph, e(5), nrm
  complex(dp) :: ref
  e = 0; nrm = 0
  do s=1,300
    call random_number(r); j = 1+int(r(1)*n1); i = 1+int(r(2)*n2); h = 1+int(r(3)*n3)
    q = eck + (j-1)*uin + (i-1)*vin + (h-1)*win
    ref = 0
    do k=1,nxat; ph = zpi*dot_product(q,xat(k,:)); ref = ref + cmplx(cos(ph),sin(ph),dp); end do
    nrm = max(nrm, abs(ref))
    e(1) = max(e(1),abs(f_ser(j,i,h)-ref)); e(2) = max(e(2),abs(f_omp(j,i,h)-ref))
    e(3) = max(e(3),abs(f_fix(j,i,h)-ref)); e(4) = max(e(4),abs(f_rec(j,i,h)-ref))
    e(5) = max(e(5),abs(f_gem(j,i,h)-ref))
  end do
  write(*,'(a,es10.2,a)') 'max |error| vs exact sum (largest |F| sampled =', nrm, '):'
  write(*,'(a,5es10.2)') '   A..E: ', e
end subroutine
end program
