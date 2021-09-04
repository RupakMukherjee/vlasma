program vlasov
implicit none

include "fftw3.f"

integer*8, parameter :: N = 256, M = 1000, nt=200000
real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0

real*8 L, dt, x, dx, kx, CFL, dt_new,mass_ratio,temp_ratio,const
real*8  ve, Veth, Vemin, Vemax, dve, Edum, sum, kee
real*8  vi, Vith, Vimin, Vimax, dvi, kei
real*8, dimension (N):: fex,fex_m,fix,fix_m
real*8, dimension (M):: fev,fev_m,fiv,fiv_m
real*8, dimension (N,M):: fexv, fixv
integer*8 i, j, p, tmax, q, t, s

double precision :: rho, E
dimension rho(N), E(N)
double complex :: rho_k, E_k
dimension rho_k(N/2+1), E_k(N/2+1)
integer*8 :: plan
common/comm/ dx,dt,dt_new

! ====================== USER INPUTS ================================

L = 2.0d0 * pi/0.4d0
dx = L / dfloat(N)
temp_ratio=1.0d0
mass_ratio=500.0d0
const=mass_ratio/temp_ratio
Veth = 1.0d0
Vemin = - 6.0d0*Veth
Vemax = + 6.0d0*Veth
dve = ( Vemax - Vemin )/ dfloat(M)

Vith = veth*dsqrt(temp_ratio/mass_ratio)
Vimin = - 6.0d0*Vith
Vimax = + 6.0d0*Vith
dvi = ( Vimax - Vimin )/ dfloat(M)

dt = 0.10d0

! Step 1.
 CFL = Vemax*dt/dx

p = floor (CFL, 8)

dt_new = dt - dx * dfloat(p) / Vemax

write(*,*) dt_new,CFL

s = 100

! ====================== INITIAL CONDITION ==========================


do j = 1,M
ve = Vemin + dfloat(j-1) * dve
vi = Vimin + dfloat(j-1) * dvi
  do i = 1,N
  x = dfloat(i-1) * dx
  fex(i) = 1.0d0+0.001d0*dcos(2.0d0*pi*x/L)
  fix(i) = 1.0d0!+0.01d0*dcos(2.0d0*pi*x/L)
  fev(j) = (1.0d0 / (dsqrt(2.0d0*pi)*veth)) * dexp( -ve*ve / (2.0d0*veth*veth) )
  fiv(j) = (1.0d0 /( dsqrt(2.0d0*pi)*vith))*dexp( -vi*vi/ (2.0d0*vith*vith) )
  fexv(i,j) = fex(i)*fev(j)
  fixv(i,j) = fix(i)*fiv(j)
  write (100,*) x,ve,vi,fexv(i,j),fixv(i,j)
  enddo
write(100,*)
enddo
!stop
! ====================== TIME ITERATION =============================


do t = 1, nt
write(*,*) "time=", t
sum = 0.0d0
kee=0.0d0
kei=0.0d0

  do j = 1,M
  ve = Vemin + dfloat(j-1) * dve
  vi = Vimin + dfloat(j-1) * dvi

    do i=1,N
    fex_m(i)=fexv(i,j)
    fix_m(i)=fixv(i,j)
    enddo

  call ppm_x (N,ve,fex_m)
  call ppm_x (N,vi,fix_m)

    do i=1,N
    fexv(i,j)=fex_m(i)
    fixv(i,j)=fix_m(i)
    enddo
  enddo

  do i = 1,N
    rho(i) = 0.0d0
  enddo

  do i = 1,N
    do j = 2,M-1
    rho(i) = rho(i) + fixv(i,j) * dvi - fexv(i,j) * dve
    enddo
  rho(i)=rho(i)+0.5d0*dvi*(fixv(i,1)+fixv(i,M))-0.5d0*dve*(fexv(i,1)+fexv(i,M))

  enddo

  do i=1,N
    do j=1,M
      ve = Vemin + dfloat(j-1) * dve
      vi = Vimin + dfloat(j-1) * dvi
      kee=kee+0.5d0*dx*dve*ve*ve*fexv(i,j)
      kei=kei+0.5d0*dx*dvi*vi*vi*fixv(i,j)
    enddo
  enddo

  call dfftw_plan_dft_r2c_1d(plan,N,rho,rho_k,FFTW_ESTIMATE)
  call dfftw_execute_dft_r2c(plan,rho,rho_k)
  call dfftw_destroy_plan(plan)

  do i = 2,N/2+1
    kx = 2.0d0 * pi * dfloat(i-1) / L
    E_k(i) = - (0.0d0,1.0d0) * rho_k(i) / kx
  enddo
  E_k(1) = rho_k(1)

write(90,*) t*dt_new,(abs(E_k(i)),i=2,N/2+1)

  call dfftw_plan_dft_c2r_1d(plan,N,E_k,E,FFTW_ESTIMATE)
  call dfftw_execute_dft_c2r(plan,E_k,E)
  call dfftw_destroy_plan(plan)

  do i = 1,N
    E(i) = E(i)/dfloat(N)
    sum = sum + 0.50d0*E(i)*E(i)*dx
  enddo

  write(10,*) t*dt_new, sum, kee, kei
  call flush(10)

  do i = 1,N
  Edum = -E(i)

    do j=1,M
      fev_m(j)=fexv(i,j)
      fiv_m(j)=fixv(i,j)
    enddo

  call ppm_v (M,Edum,dve,fev_m)
  Edum=E(i)/mass_ratio
  call ppm_v (M,Edum,dvi,fiv_m)

    do j=1,M
    fexv(i,j)=fev_m(j)
    fixv(i,j)=fiv_m(j)
    enddo
  enddo

  do j = 1,M
  ve = Vemin + dfloat(j-1) * dve
  vi = Vimin + dfloat(j-1) * dvi
    do i=1,N
      fex_m(i)=fexv(i,j)
      fix_m(i)=fixv(i,j)
    enddo

  call ppm_x (N,ve,fex_m)
  call ppm_x (N,vi,fix_m)

    do i=1,N
      fexv(i,j)=fex_m(i)
      fixv(i,j)=fix_m(i)
    enddo
  enddo

enddo ! t

end program vlasov

!========================================================


!========================================================

subroutine ppm_v (M,Edum,dv,fv)
implicit none
integer*8, intent (in) :: M
real*8, intent (in):: Edum
real*8, intent (inout):: fv(M)
real*8  df, k,dt_new,dt,dx,dv
real*8, dimension (M):: diff, dmf, fL, fR, Delta_f, f_Bar, f6
integer*8 i
common/comm/ dx,dt,dt_new


! Step 2.
do i = 1,M
  if (i == 1) then
  diff(i) = fv(i) - fv(M)
  else
  diff(i) = fv(i) - fv(i-1)
  endif
enddo

! Step 3+4.
do i = 1,M
  if (i == 1) then
  df = 0.50d0 * (fv(i+1) - fv(M))
    if (diff(i) * diff(i+1) >= 0.0d0) then
    dmf(i) = min( abs(df), 2.0d0*abs(diff(i)), 2.0d0*abs(diff(i+1)) ) * dsign(1.0d0,df)
    else
    dmf(i) = 0.0d0
    endif
  elseif (i == M) then
  df = 0.50d0 * (fv(1) - fv(i-1))
    if (diff(i) * diff(1) >= 0.0d0) then
    dmf(i) = min( abs(df), 2.0d0*abs(diff(i)), 2.0d0*abs(diff(1)) ) * dsign(1.0d0,df)
    else
    dmf(i) = 0.0d0
    endif
  else
  df = 0.50d0 * (fv(i+1) - fv(i-1))
    if (diff(i) * diff(i+1) >= 0.0d0) then
    dmf(i) = min( abs(df), 2.0d0*abs(diff(i)), 2.0d0*abs(diff(i+1)) ) * dsign(1.0d0,df)
    else
    dmf(i) = 0.0d0
    endif
  endif
enddo

! Step 5.
do i = 1,M

  if (i == M) then
  fR(i) = fv(i) + 0.50d0 * diff(1) + (1.0d0/6.0d0) * (dmf(i) - dmf(1))
  fL(1) = fR(i)
  else
  fR(i) = fv(i) + 0.50d0 * diff(i+1) + (1.0d0/6.0d0) * (dmf(i) - dmf(i+1))
  fL(i+1) = fR(i)
  endif

enddo

! Step 6+7.
do i = 1,M
Delta_f(i) = fR(i) - fL(i)

  if ((fR(i)-fv(i)) * (fv(i)-fL(i)) <= 0.0d0) then
  fL(i) = fv(i)
  fR(i) = fv(i)
  elseif (Delta_f(i) * (fv(i) - 0.50d0*(fL(i)+fR(i))) > + Delta_f(i)*Delta_f(i)/6.0d0) then
  fL(i) = 3.0d0*fv(i) - 2.0d0*fR(i)
  elseif (Delta_f(i) * (fv(i) - 0.50d0*(fL(i)+fR(i))) < - Delta_f(i)*Delta_f(i)/6.0d0) then
  fR(i) = 3.0d0*fv(i) - 2.0d0*fL(i)
  endif

enddo

! Step 8+9.
do i = 1,M
Delta_f(i) = fR(i) - fL(i)
f6(i) = 6.0d0 * (fv(i) - 0.50d0 * (fL(i)+fR(i)))
enddo

! Step 10.

k = Edum * dt_new / dv

! Step 11.
do i = 1,M
  if (i == M) then
    if (k > 0) then
    f_Bar(i) = fR(i) - (k/2.0d0) * (Delta_f(i) - (1.0d0-2.0d0*k/3.0d0) * f6(i))
    else
    f_Bar(i) = fL(1) - (k/2.0d0) * (Delta_f(1) + (1.0d0+2.0d0*k/3.0d0) * f6(1))
    endif
  else
    if (k > 0) then
    f_Bar(i) = fR(i) - (k/2.0d0) * (Delta_f(i) - (1.0d0-2.0d0*k/3.0d0) * f6(i))
    else
    f_Bar(i) = fL(i+1) + (- (k/2.0d0)) * (Delta_f(i+1) + (1.0d0-(-2.0d0*k/3.0d0)) * f6(i+1))
    endif
  endif
enddo

! Step 12.
do i = 1,M
  if (i == 1) then
  fv(i) = ( fv(i) + k * ( f_Bar(M)-f_Bar(i)) )
  else
  fv(i) = ( fv(i) + k * ( f_Bar(i-1)-f_Bar(i)) )
  endif
enddo

end subroutine ppm_v

!========================================================
subroutine ppm_x (N,v,fx)
implicit none
integer*8, intent (in) :: N
real*8, intent (in):: v
real*8, intent (inout):: fx(N)
real*8 df, k,dx,dt,dt_new
real*8, dimension (N):: diff, dmf, fL, fR, Delta_f, f_Bar, f6
integer*8 i
common/comm/ dx,dt,dt_new

! Step 2.
do i = 1,N
  if (i == 1) then
  diff(i) = fx(i) - fx(N)
  else
  diff(i) = fx(i) - fx(i-1)
  endif
enddo

! Step 3+4.
do i = 1,N
  if (i == 1) then
  df = 0.50d0 * (fx(i+1) - fx(N))
    if (diff(i) * diff(i+1) >= 0.0d0) then
    dmf(i) = min( abs(df), 2.0d0*abs(diff(i)), 2.0d0*abs(diff(i+1)) ) * dsign(1.0d0,df)
    else
    dmf(i) = 0.0d0
    endif
  elseif (i == N) then
  df = 0.50d0 * (fx(1) - fx(i-1))
    if (diff(i) * diff(1) >= 0.0d0) then
    dmf(i) = min( abs(df), 2.0d0*abs(diff(i)), 2.0d0*abs(diff(1)) ) * dsign(1.0d0,df)
    else
    dmf(i) = 0.0d0
    endif
  else
  df = 0.50d0 * (fx(i+1) - fx(i-1))
    if (diff(i) * diff(i+1) >= 0.0d0) then
    dmf(i) = min( abs(df), 2.0d0*abs(diff(i)), 2.0d0*abs(diff(i+1)) ) * dsign(1.0d0,df)
    else
    dmf(i) = 0.0d0
    endif
  endif
enddo

! Step 5.
do i = 1,N

  if (i == N) then
  fR(i) = fx(i) + 0.50d0 * diff(1) + (1.0d0/6.0d0) * (dmf(i) - dmf(1))
  fL(1) = fR(i)
  else
  fR(i) = fx(i) + 0.50d0 * diff(i+1) + (1.0d0/6.0d0) * (dmf(i) - dmf(i+1))
  fL(i+1) = fR(i)
  endif

enddo

! Step 6+7.
do i = 1,N
Delta_f(i) = fR(i) - fL(i)

  if ((fR(i)-fx(i)) * (fx(i)-fL(i)) <= 0.0d0) then
  fL(i) = fx(i)
  fR(i) = fx(i)
  elseif (Delta_f(i) * (fx(i) - 0.50d0*(fL(i)+fR(i))) > + Delta_f(i)*Delta_f(i)/6.0d0) then
  fL(i) = 3.0d0*fx(i) - 2.0d0*fR(i)
  elseif (Delta_f(i) * (fx(i) - 0.50d0*(fL(i)+fR(i))) < - Delta_f(i)*Delta_f(i)/6.0d0) then
  fR(i) = 3.0d0*fx(i) - 2.0d0*fL(i)
  endif

enddo

! Step 8+9.
do i = 1,N
Delta_f(i) = fR(i) - fL(i)
f6(i) = 6.0d0 * (fx(i) - 0.50d0 * (fL(i)+fR(i)))
enddo

! Step 10.

k = v * dt_new/ (2.0d0 * dx)

! Step 11.
do i = 1,N
  if (i == N) then
    if (k > 0) then
    f_Bar(i) = fR(i) - (k/2.0d0) * (Delta_f(i) - (1.0d0-2.0d0*k/3.0d0) * f6(i))
    else
    f_Bar(i) = fL(1) - (k/2.0d0) * (Delta_f(1) + (1.0d0+2.0d0*k/3.0d0) * f6(1))
    endif
  else
    if (k > 0) then
    f_Bar(i) = fR(i) - (k/2.0d0) * (Delta_f(i) - (1.0d0-2.0d0*k/3.0d0) * f6(i))
    else
    f_Bar(i) = fL(i+1) + (- (k/2.0d0)) * (Delta_f(i+1) + (1.0d0-(-2.0d0*k/3.0d0)) * f6(i+1))
    endif
  endif
enddo

! Step 12.
do i = 1,N
  if (i == 1) then
  fx(i) = ( fx(i) + k * ( f_Bar(N)-f_Bar(i)) )
  else
  fx(i) = ( fx(i) + k * ( f_Bar(i-1)-f_Bar(i)) )
  endif
enddo

end subroutine ppm_x
