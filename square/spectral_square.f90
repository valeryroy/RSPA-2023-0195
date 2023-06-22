program spectral_square_reg
! usage: 
!        ./spectral_square_reg 0.5 [a] 6 [Ng] 40 [Ndpt] 10 [Nmom]
!
!  
! This code computes the iterate Gamma^n = Gamma o Gamma o ... o Gamma applied on u0= x or y
! and the moments mu[n+1]= < x, Gamma^n x>
!
! The unit cell is composed of:
!
! a matrix of conductivity sigma0
! inclusions (Npart) of conductivity sigma1 
!
! Let Gamma be the boundary between the 2 phases
!  
! with: 
!      s = sigma1/sigma0 -1
!      f = <theta>
!      n (n1,n2): unit outer normal to Gamma at P
!      e1=(1,0), e1= (0,1)      
!
!     G(P,Q) is the doubly periodic Green's function: Delta G(P,Q)= -delta(P-Q)+1
!
!
! The effective conductivity is then 
!
!       (sigma^*)_11 = 1 + f*s + <Gamma x, x> s^2 + <Gamma^2 x, x> s^3 + ...
!       (sigma^*)_22 = 1 + f*s + <Gamma y, y> s^2 + <Gamma^2 y, y> s^3 + ...
!  
  ! Inputs on command line:
  !          a: square side length
!            Ng: Number of series coefficients (integer)
!            Ndpt: Number of boundary points /unit length (integer)
!            Nmom: Number of moments (integer)
!
!
! Output: the effective conductivity sigma*_11 (or sigma*_22) as a series expansion
!
!          mu(1), ..., mu(Nmom)
!  
! Note #1:  code uses a series representation of G_n (P,Q) with fast converging properties. See subroutine GradGreen.
!           
! Note #2: quadrature is regular trapezoidal
!
! Note #3: the series must be summed by Pade approximants to accelerate the series.
! use: maple < pade_spect.mw  (the input data file is "moments.dat")
! Maple needs to be installed.
!
! Author: R.V. Roy, 07/21/2022
!
!....................................................................
!  
implicit double precision (a-h,o-z)
double precision, parameter::pi=3.141592653589793238d0
double precision, allocatable :: x(:),y(:),w(:),anormx(:),anormy(:),rho0(:),rho1(:),amu(:),amat(:,:)
CHARACTER(100) :: arg,num1char,num2char,num3char,num4char

double precision :: cc(10)
double complex :: tau,q,q2,qq,imagi,qc(10)

! output file for the moments

open(unit=7,file='moments.dat',status='unknown')

!
! Step 1: read parameters/data from terminal
!

narg = COMMAND_ARGUMENT_COUNT()

  if (narg < 4) then
     write(*,*) ""
     write(*,*) "Not a Multithread Version"
     write(*,*) ""
     write(*,*) "Usage:"
     write(*,*) "./spectral_square_reg  <a>  <Ng> <Ndpt>  <Nmom>"
     write(*,*) "   "
     write(*,*) " a: side length of square inclusion"
     write(*,*) " Ng : Number of Green function coefficients"
     write(*,*) " Ndpt : Number of boundary points per unit length"
     write(*,*) " Nmom: Number of moments (even integer)"
     write(*,*) " Note: quadrature is regular trapezoidal"
     write(*,*) " "
     stop
  end if
  
CALL GET_COMMAND_ARGUMENT(1,num1char)   !first, read in the 4 input values
CALL GET_COMMAND_ARGUMENT(2,num2char)
CALL GET_COMMAND_ARGUMENT(3,num3char)
CALL GET_COMMAND_ARGUMENT(4,num4char)



READ(num1char,*)a    ! square side length
READ(num2char,*)Ng   ! lattice size
READ(num3char,*)Npt ! number of points/side
READ(num4char,*)Nmom ! number of "moments"

wtime = omp_get_wtime ( )  ! wall time



! allocate arrays

allocate (amu(0:Nmom))
!
! compute series coefficients cc[i]:= q^(2*i)/(1-q^(2*i)), q=exp(-pi)

taux=0.d0
tauy=1.d0
tau=dcmplx(taux,tauy)
imagi=dcmplx(0.d0,1.d0)

q= exp(pi*imagi*tau)
print *,'q= ',q

area=tauy

do k=1,5
   q2=q**2
   qq= q2**k
   qc(k)= qq
end do

   qc(6) = dcmplx(0.42411511830160775422d-16,0.d0)
   qc(7) = dcmplx(0.79201069507981122597d-19,0.d0)
   qc(8) = dcmplx(0.14790346159617856517d-21,0.d0)
   qc(9) = dcmplx(0.27620124435223531531d-24,0.d0)
   qc(10)= dcmplx(0.51579000625428403468d-27,0.d0)

!
!     generate boundary points (x(k),y(k)) k=1..Npt1 (inner boundaries(
!
   
kk=0

dx= a/dfloat(2*NPT)
NP= 2*NPT-1
print *,'Npt = ',Np,' boundary points/side'

allocate(x(4*NP),y(4*NP),w(4*NP),anormx(4*NP),anormy(4*NP),rho0(4*NP),rho1(4*NP),amat(4*NP,4*NP))

do k=1,NP ! loop over # of boundary points side 1
   kk=kk+1
   anormx(kk)=1.d0 ! x-normal at (x(kk),y(kk))
   anormy(kk)=0.d0 ! y-normal at (x(kk),y(kk))
   x(kk) =   0.5d0*a  ! x-boundary point
   y(kk) =  -0.5d0*a + dfloat(k)*dx  ! y-boundary point
   w(kk) =   dx !  weight
   rho0(kk) =  x(kk) ! set rho0= x
end do

do k=1,NP ! loop over # of boundary points side 2
   kk=kk+1
   anormx(kk)=0.d0 ! x-normal at (x(kk),y(kk))
   anormy(kk)=1.d0 ! y-normal at (x(kk),y(kk))
   y(kk) =   0.5d0*a ! x-boundary point
   x(kk) =   y(NP-k+1) ! y-boundary point
   w(kk) =    dx !  weight
   rho0(kk) =  x(kk) ! set rho0= x
end do

do k=1,NP ! loop over # of boundary points side 3
   kk=kk+1
   anormx(kk)=-1.d0 ! x-normal at (x(kk),y(kk))
   anormy(kk)=0.d0 ! y-normal at (x(kk),y(kk))
   x(kk) =   -0.5d0*a ! x-boundary point
   y(kk) =   y(NP-k+1) ! y-boundary point
   w(kk) =   dx !  weight
   rho0(kk) =  x(kk) ! set rho0= x
end do

do k=1,NP ! loop over # of boundary points side 4
   kk=kk+1
   anormx(kk)=0.d0 ! x-normal at (x(kk),y(kk))
   anormy(kk)=-1.d0 ! y-normal at (x(kk),y(kk))
   y(kk) =   -0.5d0*a ! x-boundary point
   x(kk) =  y(k)  ! y-boundary point
   w(kk) =  dx  !  weight
   rho0(kk) =  x(kk) ! set rho0= x
end do

NPTS=kk

print *,'Npts = ',Npts,' total boundary points'
write(7,*)'Npts := ',Npts,'; # total boundary points'


! 
! find amu[1]= < theta> = f
!
amu(0)=1.d0 
amu(1)=0.d0 

do i=1,Npts
      amu(1)=amu(1)+ (x(i)*anormx(i)+y(i)*anormy(i))*w(i)
end do

amu(1)= 0.5d0*amu(1)/area
print *," "

! compare with area fraction

f=a**2

print *,'f:=',f,';'
print *,'mu[1]:=',amu(1),';'

!
! compute interaction matrix from gradient of Green's function
!

do i=1,Npts   
   do j=1,Npts
      if(j.eq.i)then
         amat(i,i)= 0.d0
      else
         xij=pi*(x(i)-x(j))
         yij=pi*(y(i)-y(j))

         call GradGreen(qc,Ng,xij,yij,dGdx,dGdy)
         dGdy=dGdy+(y(i)-y(j))/tauy
         amat(i,j)=dGdx*anormx(j)+dGdy*anormy(j)
      end if
   end do
end do

! main loop

do m=2,Nmom 

   do i=1,Npts
      rho1(i)=0.0
      do j=1,Npts
         rho1(i)=rho1(i)+amat(i,j)*w(j)*rho0(j)
      end do
      rho1(i)= rho1(i)-0.5d0*rho0(i)
   end do
   
!!!!!!!!!!!!!!!!!!!
! find mth moment
!!!!!!!!!!!!!!!!!!!
   amu(m)=0.0 
   do i=1,Npts
      amu(m)=amu(m)+ anormx(i)*rho1(i)*w(i)
   end do
   amu(m)=amu(m)/area
   write(*,*)'mu[',m,']:=',amu(m),';'
!
! set rho0 <- rho1
!     
   do i=1,Npts
      rho0(i)=rho1(i)
   end do
!         
end do                    ! end of loop
!
!     print moments
!
write(7,*)'M := ',Nmom,';'
write(7,*)'f := ',f,';'
write(7,*)'Npts := ',Npts,';'

do k=0,Nmom
   write(7,*)'mu[',k,'] := ',amu(k),';'
end do

!
wtime = omp_get_wtime ( ) - wtime
write(*,*)' '
write(*,*)'Exec time  = ',wtime,'seconds'   
end program spectral_square_reg

subroutine GradGreen(qc,N,dx,dy,dGdx,dGdy)
  implicit double precision (a-h,o-z)
  
  !compute dG=cot(z)+4*sum(qq^i*sin(2*z)/(1-2.0*qq^i*cos(2*z)+qq^(2*i)), i= 1..N)
  ! z=dx+ I*dy
  ! dGdx= Re(dG), dGdy= Im(dG)
  
  double complex :: qc(10),z,coz,snz,sin2z,cos2z,dG,dGb
  z=dcmplx(dx,dy)
  coz=cos(z)
  snz=sin(z)
  sin2z=2.d0*snz*coz
  cos2z=2.d0*coz**2-1.d0
  dG= coz/snz
  dGb=cmplx(0.d0,0.d0)
  do k=1,N
     dGb=dGb+ qc(k)*sin2z/(1.d0-2.d0*qc(k)*cos2z+qc(k)**2)
  end do
  dG=dG+4.d0*dGb
  dGdx= -0.5d0*real(dG)
  dGdy=  0.5d0*aimag(dG)
end subroutine GradGreen
