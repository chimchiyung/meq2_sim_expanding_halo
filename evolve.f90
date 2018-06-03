program fill
 implicit none
 integer::nparticle,navg=1*10**4,navgdqt=10**5,navgdf=5*10**3,irec
 integer,parameter::nintpl=51
 integer::printPosI=0
 real(8)::pi=acos(-1._8),q2,q2t=0._8,q2t0=0._8,q2tconst
 real(8)::omega,domega=0._8,domegat=0._8
 real(8)::printPosQ=.07_8
 real(8)::t=0._8,tend=1._8*10**6,dt=.1_8
 real(8),dimension(:),allocatable::pth,th
 real(8),dimension(:),allocatable::q2tRecord, domegaRecord
 real(8),dimension(nintpl)::pin,vqEr,efn,be,ce,de,bv,cv,dv
 real(8),dimension(7,7)::b
 real(8),dimension(7)::a,c
 real(8)::Pt,G,pmin,den,lambda,nh,Pw,D0

 call initran()
 call readparameter()
 call readpos()
 call loadrk6consts()
 call readefn()
 call readfreq()
 call readvqEr()
 call readq2tconst()
 call spline(pin,vqEr,bv,cv,dv,nintpl)
 call spline(pin,efn,be,ce,de,nintpl)
 write(*,*) size(q2tRecord)
 open(50,file="q2.dat")
 write(50,*) t,q2,q2t,q2t0,domega
 
 irec=1
! do while (t<3._8*10**5)
 do while ((.12_8>q2).and.(q2>.01_8).and.(t<tend))
  !q2t=((navg-1)*q2t+q2tint())/navg
  q2t0=((navgdqt-1)*q2t0+q2tint())/navgdqt
  q2tRecord(1:navg-1)=q2tRecord(2:navg)
  domegaRecord(1:navgdf-1)=domegaRecord(2:navgdf)
  q2tRecord(navg)=q2tint()
  domegaRecord(navgdf)=domegaint()
  q2t=sum(q2tRecord)/navg
  !domegat=(domegaint()-domega)/navgdf/dt
  domegat=(sum(domegaRecord)/navgdf-domega)/dt
  call pqnew(pth,th)
  call addpos()
  q2=q2+q2t*dt
  domega=domega+domegat*dt
  t=t+dt
  irec=irec+1
  if (mod(irec,10**5)==0) then
   write(50,*) t,q2,q2t,q2t0,domega
   irec=0
  end if
  if ((q2<printPosQ).and.(printPosI==0)) then
   call writepos(1)
   printPosI=1
  end if
 end do
 call writepos(0)
 close(50)
 deallocate(pth)
 deallocate(th)
 deallocate(q2tRecord)
 deallocate(domegaRecord)

contains

 subroutine addpos()
  integer::iadd,jadd,naddsum,ntemp,nadd
  real(8)::radd=0._8
  real(8),dimension(:),allocatable::padd,qadd
  real(8),dimension(:),allocatable::ptemp,qtemp

  nadd=floor(G*dt/2)
  radd=radd+G*dt/2-floor(G*dt/2)
  if (radd>1) then
   nadd=nadd+1
   radd=0._8
  end if
  allocate(padd(nadd),qadd(nadd))
  call random_number(qadd)
  qadd=pi*qadd
  padd=pmin

  ntemp=nparticle+nadd
  allocate(ptemp(ntemp),qtemp(ntemp))
  ptemp(1:nparticle)=pth
  qtemp(1:nparticle)=th
  ptemp(nparticle+1:ntemp)=padd
  qtemp(nparticle+1:ntemp)=qadd
  deallocate(pth)
  deallocate(th)
  nparticle=ntemp
  allocate(pth(nparticle),th(nparticle))
  pth=ptemp
  th=qtemp
  deallocate(ptemp)
  deallocate(qtemp)

  deallocate(padd)
  deallocate(qadd)

 end subroutine

 subroutine rmpos()
  integer::irm,jrm
  integer,dimension(:),allocatable::prm
  real(8),dimension(:),allocatable::ptemp,qtemp
  allocate(prm(size(pth)))
  where (pth>Pw)
   prm=1
  elsewhere
   prm=0
  end where
  nparticle=nparticle-sum(prm)
  allocate(ptemp(size(prm)),qtemp(size(prm)))
  ptemp=pth
  qtemp=th
  deallocate(pth)
  deallocate(th)
  allocate(pth(nparticle),th(nparticle))
  jrm=1
  do irm=1,size(prm)
   if (prm(irm)==0) then
    pth(jrm)=ptemp(irm)
    th(jrm)=qtemp(irm)
    jrm=jrm+1
   end if
  end do
  deallocate(prm)
  deallocate(ptemp)
  deallocate(qtemp)
 end subroutine

 subroutine pqnew(p0,q0)
  real(8),dimension(:)::p0,q0
  real(8),dimension(size(q0))::p1,q1
  real(8),dimension(size(q0))::dpRandom
  real(8),dimension(size(q0),7)::dp,dq
  integer::l1,l2
!  do l2=1,size(p0)
!   p0(l2)=p0(l2)+sqrt(2*difffn(p0(l2))*dt)*rannum()
!  end do
  do l2=1,size(p0)
   dpRandom(l2)=sqrt(2*difffn(p0(l2))*dt)*rannum()
  end do
!  pth=abs(pth)
  p1=p0
  q1=q0
  dp=0._8
  dq=0._8
  do l1=1,7
   do l2=1,size(q0)
    dp(l2,l1)=(vp(p1(l2),q1(l2),c(l1))&
               +difffndp(p1(l2))+Ptfn(sqrt(p1(l2)/Pw)))*dt
!    dp(l2,l1)=vp(p1(l2),q1(l2),c(l1))*dt
    dq(l2,l1)=vq(p1(l2),q1(l2),c(l1))*dt
   end do
   if (l1<7) then
    p1=p0
    q1=q0
    do l2=1,l1
     p1=p1+b(l2,l1+1)*dp(:,l2)
     q1=q1+b(l2,l1+1)*dq(:,l2)
    end do
   end if
  end do
  do l1=1,7
   p0=p0+a(l1)*dp(:,l1)
   q0=q0+a(l1)*dq(:,l1)
  end do
  q0=mod(q0+pi,pi)
  p0=abs(p0+dpRandom)
 end subroutine

 function q2tint()
  real(8)::q2tint
  integer::iq2t
  q2tint=0._8
  do iq2t=1,nparticle
   q2tint=q2tint+vp(pth(iq2t),th(iq2t),0._8)
  end do
  q2tint=q2tint*lambda/q2/4*q2tconst
 end function

 function domegaint()
  real(8)::domegaint
  integer::ido
  domegaint=0._8
  do ido=1,nparticle
   domegaint=domegaint-vp(pth(ido),th(ido),0._8)/tan(2*th(ido))
  end do
  domegaint=domegaint*lambda/q2**2/4*q2tconst
 end function

 function vp(p0,q0,cp)
  real(8)::p0,q0,vp,cp
  vp=ispline(p0,pin,efn,be,ce,de,nintpl,0)*2*sin(2*q0)*(q2+q2t*cp*dt)
 end function

 function vq(p0,q0,cq)
  real(8)::p0,q0,vq,cq
  vq=ispline(p0,pin,vqEr,bv,cv,dv,nintpl,0)-(omega+domega+domegat*cq*dt)/2&
     +ispline(p0,pin,efn,be,ce,de,nintpl,1)*cos(2*q0)*(q2+q2t*cq*dt)
 end function

 subroutine loadrk6consts()
  integer::iload
  open(30,file="rk6consts.dat",status="old")
  read(30,*) a(:)
  do iload=1,7
   read(30,*) b(iload,:)
  end do
  close(30)
  do iload=1,7
   c(iload)=sum(b(:,iload))
  end do
 end subroutine

 subroutine writepos(iIntermediate)
  integer::iwrite,iIntermediate
  if (iIntermediate==1) then
   open(30,file="posIntermediate.dat",status="replace")
  else
   open(30,file="pos.dat",status="replace")
  end if
  write(30,*) nparticle
  write(30,*) q2,domega
  do iwrite=1,nparticle
   write(30,*) pth(iwrite), th(iwrite)
  end do
  close(30)
 end subroutine

 subroutine readpos()
  integer::iread
  allocate(q2tRecord(navg))
  q2tRecord = 0._8
  allocate(domegaRecord(navgdf))
  open(30,file="posstart.dat",status="old")
  read(30,*) nparticle
  read(30,*) q2,domega
  domegaRecord = domega
  allocate(pth(nparticle),th(nparticle))
  do iread=1,nparticle
   read(30,*) pth(iread), th(iread)
  end do
  close(30)
 end subroutine

 subroutine readefn()
  real(8)::pread
  integer::iread
  open(31,file="dphitab.dat")
  read(31,*) pread,efn(1)
  do iread=2,nintpl
   read(31,*)
   read(31,*) pread,efn(iread)
  end do
  close(31)  
 end subroutine

 subroutine readvqEr()
  real(8)::pread
  integer::iread
  open(31,file="dthdtEr.dat")
  read(31,*) pin(1),vqEr(1)
  do iread=2,nintpl
   read(31,*)
   read(31,*) pin(iread),vqEr(iread)
  end do
  close(31)  
 end subroutine

 subroutine readfreq()
  open(31,file="omega.dat")
  read(31,*) omega
  close(31)
  
 end subroutine

 subroutine readq2tconst()
  open(31,file="dampconst.dat")
  read(31,*) q2tconst
  close(31)
 end subroutine

 function difffn(p1)
  real(8)::difffn,p1,p0=.9_8,dp0=.1_8
  difffn=D0*stepfn(p1-p0,dp0)
 end function

 function difffndp(p1)
  real(8)::difffndp,p1,p0=.9_8,dp0=.1_8
  difffndp=D0/2/dp0/cosh((p1-p0)/dp0)
 end function

 function Ptfn(x)
  real(8)::Ptfn,x
  Ptfn=Pt*nh/nhgen(x)
 end function

 function nhgen(x)
  real(8)::nhgen,x
  real(8)::nc0=1.5_8,rnh=.017_8,Rp=.33_8
  nhgen=nc0*rnh*(1-.75_8*x*stepfn(x-Rp,.05_8))&
               *(.8_8*stepfn(x-Rp,.03_8)+.2_8)
 end function

 function stepfn(x,del)
  real(8)::stepfn,x,del
  stepfn=(1+tanh(x/del))/2
 end function

 subroutine readparameter()
  open(31,file="parameter.dat",status="old")
  read(31,*)
  read(31,*) Pt
  read(31,*)
  read(31,*) G
  read(31,*)
  read(31,*) pmin
  read(31,*)
  read(31,*)
  read(31,*)
  read(31,*) nh
  read(31,*)
  read(31,*) Pw
  close(31)
  den=G/(2*pi*Pt)
  lambda=nh*Pt/G
  Pw=1/Pw
  !D0=Pt*3*.1_8**3
  D0=0._8
 end subroutine

 function rannum()
  real(8)::rannum,s,r
  real(8),save::x1,x2
  integer::ipair=0
  if (ipair==1) then
   rannum=x2
   ipair=0
  else
   s=2
   do while (s>1)
    call random_number(x1)
    call random_number(x2)
    x1=2*x1-1
    x2=2*x2-1
    s=x1**2+x2**2
   end do
   r=sqrt(-2*log(s)/s)
   x1=x1*r
   x2=x2*r
   rannum=x1
   ipair=1
  end if
 end function

 subroutine initran()
  integer::i,seedsize
  integer,allocatable::seedval(:)
  call random_seed(size=seedsize)
  allocate(seedval(seedsize))
  do i=1,seedsize
   call system_clock(count=seedval(i))
  end do
  call random_seed(put=seedval)
  deallocate(seedval)
 end subroutine

   subroutine spline (x, y, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
implicit none
integer n
double precision x(n), y(n), b(n), c(n), d(n)
integer i, j, gap
double precision h

gap = n-1
! check input
if ( n < 2 ) return
if ( n < 3 ) then
  b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
  c(1) = 0.
  d(1) = 0.
  b(2) = b(1)
  c(2) = 0.
  d(2) = 0.
  return
end if
!
! step 1: preparation
!
d(1) = x(2) - x(1)
c(2) = (y(2) - y(1))/d(1)
do i = 2, gap
  d(i) = x(i+1) - x(i)
  b(i) = 2.0*(d(i-1) + d(i))
  c(i+1) = (y(i+1) - y(i))/d(i)
  c(i) = c(i+1) - c(i)
end do
!
! step 2: end conditions 
!
b(1) = -d(1)
b(n) = -d(n-1)
c(1) = 0.0
c(n) = 0.0
if(n /= 3) then
  c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
  c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
  c(1) = c(1)*d(1)**2/(x(4)-x(1))
  c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
end if
!
! step 3: forward elimination 
!
do i = 2, n
  h = d(i-1)/b(i-1)
  b(i) = b(i) - h*d(i-1)
  c(i) = c(i) - h*c(i-1)
end do
!
! step 4: back substitution
!
c(n) = c(n)/b(n)
do j = 1, gap
  i = n-j
  c(i) = (c(i) - d(i)*c(i+1))/b(i)
end do
!
! step 5: compute spline coefficients
!
b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
do i = 1, gap
  b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
  d(i) = (c(i+1) - c(i))/d(i)
  c(i) = 3.*c(i)
end do
c(n) = 3.0*c(n)
d(n) = d(n-1)
end subroutine spline

  function ispline(u, x, y, b, c, d, n, idiff)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!=======================================================================
implicit none
double precision ispline
integer n, idiff
double precision  u, x(n), y(n), b(n), c(n), d(n)
integer i, j, k
double precision dx

! if u is ouside the x() interval take a boundary value (left or right)
!if(u <= x(1)) then
!  ispline = y(1)
!  return
!end if
!if(u >= x(n)) then
!  ispline = y(n)
!  return
!end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
!i = 1
!j = n+1
!do while (j > i+1)
!  k = (i+j)/2
!  if(u < x(k)) then
!    j=k
!    else
!    i=k
!   end if
!end do
i=min(nintpl,max(1,ceiling((nintpl-1)*u/Pw)))
!*
!  evaluate spline interpolation
!*
dx = u - x(i)

if (idiff==0) then
 ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
else
 ispline = b(i) + dx*(c(i)*2 + dx*d(i)*3)
end if
end function ispline

end program
