program spectral_disk
! usage: 
!        ./spectral_disk  0.4 [f] 6 [Ng] 40 [Ndpt] 10 [Nmom]
!
!  
! This code computes the iterate Gamma^n = Gamma o Gamma o ... o Gamma applied on u0= x 
! and the moments mu[n+1]= < x, Gamma^n x>
!
! The unit cell is composed of:
!
! a matrix of conductivity sigma0
! disk inclusion of conductivity sigma1 
!
! Let Gamma be the boundary between the 2 phases
!  
! with: 
!      s = sigma1/sigma0 -1
!      f = <theta>
!      n (n1,n2): unit outer normal to Gamma at P
!      e1=(1,0)     
!
!     G(P,Q) is the doubly periodic Green's function: Delta G(P,Q)= -delta(P-Q)+1
!
!
! The effective conductivity is then 
!
!       (sigma^*/sigma0) = 1 + f*s + <Gamma x, x> s^2 + <Gamma^2 x, x> s^3 + ...
!       
!  
! Inputs on command line:
!            geom: S (square) , T (hex)
!            Ng: Number of series coefficients (integer)
!            Ndpt: Number of boundary points /unit length (integer)
!            Nmom: Number of moments (integer)
!            filename: name of input file (character)
!
! Output: the effective conductivity sigma* as a series expansion
!
!          mu(1), ..., mu(Nmom)
!  
! Note #1:  code uses a series representation of G_n (P,Q) with fast converging properties. See subroutine GradGreen.
!           to compile, type: make spectral_disk
!
!
! Note #2: the series must be summed by Pade approximants to accelerate the series.
! use: maple < pade_spect.mw  (the input data file is "moments.dat")
! Maple needs to be installed.
!
! Author: R.V. Roy, 07/26/2022
!                   09/19/2022: added hexagonal array geometry.
!                             verified sigma* inthe limit sigma1/sigma0 ->infty
!....................................................................
!  
implicit double precision (a-h,o-z)
double precision, parameter::pi=3.141592653589793238d0
double precision, allocatable :: theta(:),x(:),y(:),w(:),anormx(:),anormy(:),rho0(:),rho1(:),amu(:),amat(:,:),ak(:)
CHARACTER(100) :: geom,arg,num1char,num2char,num3char,num4char,num5char

double precision :: cc(10)
double complex :: tau,q,q2,qq,imagi,qc(10)

! output file for the moments

open(unit=7,file='moments.dat',status='unknown')
open(unit=8,file='rho.dat',status='unknown')

!
! Step 1: read parameters/data from terminal
!

narg = COMMAND_ARGUMENT_COUNT()

  if (narg < 5) then
     write(*,*) ""
     write(*,*) "Not a Multithread Version"
     write(*,*) ""
     write(*,*) "Usage:"
     write(*,*) "./spectral_disk <geom> <f>  <Ng> <Ndpt>  <Nmom>"
     write(*,*) "   "
     write(*,*) " geom: S/T (square/hexagonal)  "
     write(*,*) " f: area fraction"
     write(*,*) " Ng : Number of Green function coefficients"
     write(*,*) " Ndpt : Number of boundary points per unit length"
     write(*,*) " Nmom: Number of moments (even integer)"
     write(*,*) " "
     stop
  end if
    
CALL GET_COMMAND_ARGUMENT(1,num1char)   !first, read in the 5 input values
CALL GET_COMMAND_ARGUMENT(2,num2char)
CALL GET_COMMAND_ARGUMENT(3,num3char)
CALL GET_COMMAND_ARGUMENT(4,num4char)
CALL GET_COMMAND_ARGUMENT(5,num5char)


READ(num1char,*)geom    ! geometry of array
READ(num2char,*)f    ! area fraction
READ(num3char,*)Ng   ! lattice size
READ(num4char,*)Ndpt ! number of points /unit length
READ(num5char,*)Nmom ! number of "moments"



wtime = omp_get_wtime ( )  ! wall time
!
! compute series coefficients cc[i]:= q^(2*i)/(1-q^(2*i)), q=exp(-pi)

if(geom .eq. 'S')then
   rad= dsqrt(f/pi)
   taux=0.d0
   tauy=1.d0
   tau=dcmplx(taux,tauy)
   imagi=dcmplx(0.d0,1.d0)

   q= exp(pi*imagi*tau)

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
end if

if(geom .eq. 'T')then
   sq3o2=dsqrt(3.d0)/2.d0
   rad=dsqrt(sq3o2*f/pi)
   taux=0.5d0
   tauy=sq3o2
   tau=dcmplx(taux,tauy)
   imagi=dcmplx(0.d0,1.d0)

   q= exp(pi*imagi*tau)

   area=tauy
   do k=1,5
      q2=q**2
      qq= q2**k
      qc(k)= qq
   end do
end if

perimeter=2.d0*pi*rad
Npt=Ndpt*perimeter

print *,'Npt = ',Npt,' boundary points'

! allocate arrays

allocate (theta(Npt),x(Npt),y(Npt),w(Npt),anormx(Npt),anormy(Npt),rho0(Npt),rho1(Npt),amat(Npt,Npt),ak(Npt))
allocate (amu(0:Nmom))

!
!     generate boundary points (x(k),y(k)) k=1..Npt1 (inner boundaries(
!
   
kk=0

do k=1,Npt ! loop over # of boundary points
   dte= 2.0*pi/Npt
   theta(k) = dte*float(k-1) ! angles from 0 to 2*pi
   csp=cos(theta(k))
   snp=sin(theta(k))
   anormx(k)=csp ! x-normal at (x(k),y(k))
   anormy(k)=snp ! y-normal at (x(k),y(k))
   x(k)=rad*csp
   y(k)=rad*snp
   w(k) =   rad*dte
   ak(k)=  -1.d0/rad ! curvature
      
   rho0(k) =  x(k) ! set rho0= x
   write(8,*)theta(k),rho0(k)
end do

write(8,*)

! 
! find amu[1]= < theta> = f
!
amu(0)=1.d0 
amu(1)=0.d0 

do i=1,Npt
      amu(1)=amu(1)+ (x(i)*anormx(i)+y(i)*anormy(i))*w(i)
end do

amu(1)= 0.5d0*amu(1)/area

print *,'series coefficients for sigma*'

! compare with area fraction


print *,'f:=',f,';'
print *,'mu[1]:=',amu(1),';'

!
! compute interaction matrix from gradient of Green's function
!

do i=1,Npt   
   do j=1,Npt
      if(j.eq.i)amat(i,i)= -ak(i)/(4.0*pi)
      if(j.ne.i)then
         xij=pi*(x(i)-x(j))
         yij=pi*(y(i)-y(j))

         call GradGreen(qc,Ng,xij,yij,dGdx,dGdy)
         dGdy=dGdy+(y(i)-y(j))/tauy
         amat(i,j)=dGdx*anormx(j)+dGdy*anormy(j)
      end if
   end do
end do

! main loop
id=1
do m=2,Nmom 

   do i=1,Npt
      rho1(i)=0.0
      do j=1,Npt
         rho1(i)=rho1(i)-amat(i,j)*w(j)*rho0(j)
      end do
      rho1(i)= rho1(i)+0.5d0*rho0(i)
      write(8,*)theta(i),rho1(i)
   end do
   write(8,*)
!!!!!!!!!!!!!!!!!!!
! find mth moment
!!!!!!!!!!!!!!!!!!!
   amu(m)=0.0 
   do i=1,Npt
      amu(m)=amu(m)+ anormx(i)*rho1(i)*w(i)
   end do
   id=-id
   amu(m)=dfloat(id)*amu(m)/area
   write(*,*)'mu[',m,']:=',amu(m),';'
!
! set rho0 <- rho1
!     
   do i=1,Npt
      rho0(i)=rho1(i)
   end do
!         
end do                    ! end of loop
!
!     print moments
!
write(7,*)'M := ',Nmom,';'
write(7,*)'f := ',f,';'
write(7,*)'Ndisk := ',Npart,';'
write(7,*)'Npt := ',Npt,';'

do k=0,Nmom
   write(7,*)'mu[',k,'] := ',amu(k),';'
end do

wtime = omp_get_wtime ( ) - wtime
write(*,*)' '
write(*,*)'Exec time  = ',wtime,'seconds'   
end program spectral_disk

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


  



