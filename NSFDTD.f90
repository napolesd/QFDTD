program FDTD2
implicit none
integer::n,nsteps,k,kc,ke,kstart
real*8::ra,ddx,dt,hbar,pi,melec,hbarra,lambda,sigma,vpot,ptot,kcenter,rb,E,uns
real*8,dimension(1:1000)::psi_rl,psi_im,vp
pi=4*atan(1.d0)
kstart=500
kcenter=int(kstart/2.0)
vpot=0.0!1.0  ! eV
nsteps=50000
ddx=0.1e-10
!ra=1.0/8.0
!rb=3.285e-4
ke=1000
melec=9.2e-31
hbarra=1.055e-34
dt=0.25*(melec/hbarra)*ddx**2
psi_rl=0.0
psi_im=0.0
vp=0.0
ptot=0.0

lambda=80.0
sigma=lambda*0.5

!----------NSFDTD-----------------------------
uns=sin((pi/2.0/lambda)**2)/sin(pi/lambda)
E=2.0*(pi*hbarra/lambda/ddx)**2/melec/1.602e-19
ra=0.5*uns/sin(pi/lambda)
rb=2.0*uns/E ! parece que aqui falta multiplicar por un sen(pi/lambda). Se agrega en la siguiente linea
rb=rb*sin(pi/lambda)
!----------NSFDTD-----------------------------
do k=1,kstart-1
psi_rl(k)=cos(2*pi*(k-kcenter)/lambda )*exp( -0.5*((k-kcenter)/sigma)**2 )
psi_im(k)=sin(2*pi*(k-kcenter)/lambda )*exp( -0.5*((k-kcenter)/sigma)**2 )
ptot=ptot+psi_rl(k)**2+psi_im(k)**2
enddo


do k=1,kstart
psi_rl(k)=psi_rl(k)/sqrt(ptot)
psi_im(k)=psi_im(k)/sqrt(ptot)
enddo

do k=kstart,kstart+20
vp(k)=vpot!/1.602e-19    ! vpot esta en eV,  convertir eV a Joules -----------------------------------------
enddo

write(1,*) hbarra**2*(2*pi/lambda/ddx)**2/(2*melec)/1.602e-19
kc=ke/2.0
!--------------------------------------------------------------------------
!----------------------ciclo temporal--------------------------------------
!--------------------------------------------------------------------------
do n=0,nsteps

do k=1,ke-1
psi_rl(k)=psi_rl(k)-ra*(psi_im(k+1)-2.0*psi_im(k)+psi_im(k-1))+rb*vp(k)*psi_im(k)
enddo


do k=1,ke-1
psi_im(k)=psi_im(k)+ra*( psi_rl(k+1)-2.d0*psi_rl(k) + psi_rl(k-1) )- rb*vp(k)*psi_rl(k) 
enddo

if (mod(n,100).eq.0) then
do k=1,ke-1
write(2,*)k,psi_rl(k)**2+psi_im(k)**2
enddo
write(2,*)""
write(2,*)""
endif

enddo
!--------------------------------------------------------------------------
!----------------------ciclo temporal--------------------------------------
!--------------------------------------------------------------------------
end PROGRAM FDTD2
