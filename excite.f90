program excite
 implicit none
 integer::nparticle,imain
 real(8)::pi=acos(-1._8),d2=.0_8,d2t,d2end=.05_8
 real(8)::domega=0._8,domegat=0._8
 real(8)::t=0._8,tend=200._8,dt=.1_8
 real(8),dimension(:),allocatable::pth,th
 real(8),dimension(7,7)::b
 real(8),dimension(7)::a,c
 real(8)::Pt,G,pmin,den,pr,lambda,nh,Pw
 integer::nexcite=10**4,navg=10**3

 call initran()
 call readparameter()
 call readpos()
 call loadrk6consts()
 do imain=1,nexcite*10
  d2t=(d2end-d2)/nexcite/dt
  if (imain>100) then
   domegat=(dOmegaInt()-domega)/dt/navg
  end if
  call pqnew(pth,th)
  d2=d2+d2t*dt
  domega=domega+domegat*dt
  t=t+dt
 end do
 call writepos()
 deallocate(pth)
 deallocate(th)

contains

 subroutine pqnew(p0,q0)
  real(8),dimension(:)::p0,q0
  real(8),dimension(size(q0))::p1,q1
  real(8),dimension(size(q0),7)::dp,dq
  integer::l1,l2
  p1=p0
  q1=q0
  dp=0._8
  dq=0._8
  do l1=1,6
   dp(:,l1)=vp(p1,q1,c(l1))*dt
   dq(:,l1)=vq(p1,q1,c(l1))*dt
   p1=p0
   q1=q0
   do l2=1,l1
    p1=p1+b(l2,l1+1)*dp(:,l2)
    q1=q1+b(l2,l1+1)*dq(:,l2)
   end do
  end do
  dp(:,7)=vp(p1,q1,c(l1))*dt
  dq(:,7)=vq(p1,q1,c(l1))*dt
  do l1=1,7
   p0=p0+a(l1)*dp(:,l1)
   q0=q0+a(l1)*dq(:,l1)
  end do
  q0=mod(q0+pi,pi)
 end subroutine

 function vp(p0,q0,cp)
  real(8),dimension(:)::p0,q0
  real(8),dimension(size(p0))::vp
  real(8)::cp
  integer::ip
  do ip=1,size(p0)
   vp(ip)=(d2+d2t*cp*dt)*(min(p0(ip),1._8)/max(p0(ip),1._8)-p0(ip)/Pw**2)*2*sin(2*q0(ip))
  end do
 end function

 function vq(p0,q0,cq)
  real(8),dimension(:)::p0,q0
  real(8),dimension(size(p0))::vq
  real(8)::cq
  integer::iq
  do iq=1,size(p0)
   if (p0(iq)<1) then
    vq(iq)=(-1+1/pr)&
           +(d2+d2t*cq*dt)&
            *(1-1/Pw**2)*cos(2*q0(iq))
   else
    vq(iq)=(-1/p0(iq)+1/pr)&
           +(d2+d2t*cq*dt)&
            *(-1/p0(iq)**2-1/Pw**2)*cos(2*q0(iq))
   end if
  end do
  vq=vq-(domega+cq*dt*domegat)/2
 end function

 function dOmegaInt()
  real(8)::dOmegaInt
  dOmegaInt=sum(vp(pth,th,0._8)/tan(2*th))*lambda/d2**2
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

 subroutine writepos()
  integer::iwrite
  open(30,file="posstart.dat",status="replace")
  write(30,*) nparticle
  write(30,*) d2,domega
  do iwrite=1,nparticle
   write(30,*) pth(iwrite), th(iwrite)
  end do
  close(30)
 end subroutine

 subroutine readpos()
  integer::iread
  open(30,file="posread.dat",status="old")
  read(30,*) nparticle
  allocate(pth(nparticle),th(nparticle))
  do iread=1,nparticle
   read(30,*) pth(iread), th(iread)
  end do
  close(30)
 end subroutine

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
  pr=2._8/(1+1/Pw**2)
 end subroutine

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

end program
