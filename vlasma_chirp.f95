program vlasov
implicit none

include "fftw3.f"
	
integer*8, parameter :: N = 256, M = 1024, nt=250000
real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0

real*8 L, dt, x, dx, kx, CFL, dt_new, v, Vth, Vmin, Vmax, dv, Edum, sum,ke
real*8, dimension (N):: fx,fx_m
real*8, dimension (M):: fv,fv_m
real*8, dimension (N,M):: fxv
integer*8 i, j, p, tmax, q, t
real*8 w0,E0,ks
double precision :: rho, E,w
dimension rho(N), E(N)
double complex :: rho_k, E_k
dimension rho_k(N/2+1), E_k(N/2+1)
integer*8 :: plan
common/comm/ dx,dt,dv,dt_new
! ====================== USER INPUTS ================================

L = 2.0d0 * pi/0.4d0
dx = L / dfloat(N)

Vth = 1.0d0
Vmin = - 6.0d0*Vth
Vmax = + 6.0d0*Vth
dv = ( Vmax - Vmin )/ dfloat(M)
dt = 0.10d0
ks=0.4d0
w0=0.6241
E0=0.025d0
! Step 1.
 CFL = Vmax*dt/dx 
 
p = floor (CFL, 8) 

dt_new = dt - dx * dfloat(p) / Vmax

write(*,*) dt_new,CFL

! ====================== INITIAL CONDITION ==========================


do j = 1,M
v = Vmin + dfloat(j-1) * dv
  do i = 1,N
  x = dfloat(i-1) * dx
  fx(i) = 1.0d0+0.0d0*dcos(2.0d0*pi*x/L)
  fv(j) = (1.0d0 / dsqrt(2.0d0*pi)) * dexp( -v*v / (2.0d0*Vth*Vth) )
  fxv(i,j) = fx(i)*fv(j)
!  write (100,*) fxv(i,1)
  enddo
!write(100,*)
enddo

! ====================== TIME ITERATION =============================


do t = 1, nt
write(*,*) "time=", t
sum = 0.0d0
ke=0.0d0
  do j = 1,M

  v = Vmin + dfloat(j-1) * dv

do i=1,N

fx_m(i)=fxv(i,j)

enddo

  call ppm_x (N,v,fx_m)

do i=1,N

fxv(i,j)=fx_m(i)

enddo
        
   enddo
 
  do i = 1,N
  rho(i) = 0.0d0
  enddo
  
  do i = 1,N
    do j = 2,M-1
    rho(i) = rho(i) + fxv(i,j) * dv
    enddo

  rho(i)=rho(i)+0.5d0*dv*(fxv(i,1)+fxv(i,M))
  rho(i) = rho(i) - 1.0d0 

  enddo  

do i=1,N

do j=1,M
v = Vmin + dfloat(j-1) * dv
ke=ke+0.5d0*dx*dv*v*v*fxv(i,j)
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
  
  call dfftw_plan_dft_c2r_1d(plan,N,E_k,E,FFTW_ESTIMATE)
  call dfftw_execute_dft_c2r(plan,E_k,E)
  call dfftw_destroy_plan(plan)

if (t.lt.10000) then  
  do i = 1,N
w=0.5d0*(-(t/10000)+2.0d0)
x = dfloat(i-1) * dx
  E(i) = E(i)/dfloat(N)
E(i)=E(i)+E0*(dsin(ks*x+w*t*dt_new)+dsin(ks*x-w*t*dt_new))  

sum = sum + 0.50d0*E(i)*E(i)*dx

  enddo
else
 do i = 1,N
  E(i) = E(i)/dfloat(N)

sum = sum + 0.50d0*E(i)*E(i)*dx

enddo

endif

  write(10,*) t*dt, sum,ke

  call flush(10)
  
  do i = 1,N

  Edum = E(i)

do j=1,M

fv_m(j)=fxv(i,j)

enddo


 call ppm_v (M,Edum,fv_m)

do j=1,M

fxv(i,j)=fv_m(J)

enddo


  enddo 

  do j = 1,M

  v = Vmin + dfloat(j-1) * dv

do i=1,N

fx_m(i)=fxv(i,j)

enddo


  call ppm_x (N,v,fx_m)

do i=1,N

fxv(i,j)=fx_m(i)

enddo
 

  enddo
  
enddo ! t

do j = 1,M
v = Vmin + dfloat(j-1) * dv
  do i = 1,N
  x = dfloat(i-1) * dx
  write (200,*) x,v,fxv(i,j)
  enddo
write(200,*)
enddo


end program vlasov

!========================================================


!========================================================

subroutine ppm_v (M,Edum,fv)
implicit none 
integer*8, intent (in) :: M
real*8, intent (in):: Edum
real*8, intent (inout):: fv(M)
real*8  df, k,dt_new,dt,dx,dv
real*8, dimension (M):: diff, dmf, fL, fR, Delta_f, f_Bar, f6
integer*8 i
common/comm/ dx,dt,dv,dt_new


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
real*8 df, k,dx,dt,dv,dt_new
real*8, dimension (N):: diff, dmf, fL, fR, Delta_f, f_Bar, f6
integer*8 i
common/comm/ dx,dt,dv,dt_new

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

