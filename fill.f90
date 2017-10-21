program fill
 implicit none
 integer::nparticle,nbd=10**2,irun,navg=10**5,irec
 real(8)::pi=acos(-1._8),d2,d2t=0._8
 real(8)::t=0._8,tend=1._8*10**6,dt=.1_8
 real(8),dimension(:),allocatable::pth,th
 real(8),dimension(:),allocatable::pbd,thbd
 real(8),dimension(7,7)::b
 real(8),dimension(7)::a,c
 real(8)::Pt,G,pmin,den,pr,lambda,nh,Pw

 call initran()
 call readparameter()
 call readpos()
 call loadrk6consts()
 open(50,file="d2.dat")
 write(50,*) t,d2
 irec=1
! do while (t<tend)
 do while ((.06_8>d2).and.(d2>.01_8).and.(t<tend))
  d2t=((navg-1)*d2t-sum(vp(pth,th,0._8))*lambda/(1-nh)/d2)/navg
!  d2t=((navg-1)*d2t-sum(sin(2*th))*lambda)/navg
  call pqnew(pth,th)
!  call renewbd()
!  call pqnew(pbd,thbd)
!  call rearrangebd()
!  call rmpos()
  call addpos()
  pth=abs(pth)
  d2=d2+d2t*dt
  t=t+dt
  irec=irec+1
  if (mod(irec,10**3)==0) then
   write(50,*) t,d2
   irec=0
  end if
 end do
 call writepos()
 close(50)
 deallocate(pth)
 deallocate(th)
 deallocate(pbd)
 deallocate(thbd)

contains

 subroutine addpos()
  integer::iadd,jadd,naddsum,ntemp,nadd
  real(8),dimension(size(pbd))::qlow,qhigh,plow,phigh
  real(8)::qnow,addprob,flux,radd=0._8
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
  where (pth<pmin)
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
  real(8),dimension(size(q0),7)::dp,dq
  integer::l1,l2
  p1=p0
  q1=q0
  dp=0._8
  dq=0._8
  do l1=1,6
   dp(:,l1)=(vp(p1,q1,c(l1))+Pt)*dt
   dq(:,l1)=vq(p1,q1,c(l1))*dt
   p1=p0
   q1=q0
   do l2=1,l1
    p1=p1+b(l2,l1+1)*dp(:,l2)
    q1=q1+b(l2,l1+1)*dq(:,l2)
   end do
  end do
  dp(:,7)=(vp(p1,q1,c(l1))+Pt)*dt
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
   vp(ip)=(d2+d2t*cp*dt)*(1-nh)&
           *(min(p0(ip),1._8)/max(p0(ip),1._8)-p0(ip)/Pw**2)*2*sin(2*q0(ip))
  end do
 end function

 function vq(p0,q0,cq)
  real(8),dimension(:)::p0,q0
  real(8),dimension(size(p0))::vq
  real(8)::cq
  integer::iq
  do iq=1,size(p0)
   if (p0(iq)<1) then
    vq(iq)=(1-nh)*(-1+1/pr)&
           +(d2+d2t*cq*dt)&
            *(1-nh)&
            *(1-1/Pw**2)*cos(2*q0(iq))
   else
    vq(iq)=(1-nh)*(-1/p0(iq)+1/pr)&
           +(d2+d2t*cq*dt)&
            *(1-nh)&
            *(-1/p0(iq)-1/Pw**2)*cos(2*q0(iq))
   end if
  end do
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

 subroutine rearrangebd()
  integer::istart,iend
  real(8),dimension(size(pbd))::ptmp,thtmp
  istart=1
  iend=size(pbd)
  do while ((istart.le.size(pbd)).and.(thbd(istart)<0))
   istart=istart+1
  end do
  do while ((iend.ge.0).and.(thbd(iend)>pi))
   iend=iend-1
  end do
  ptmp(1:size(pbd)-iend)=pbd(iend+1:size(pbd))
  ptmp(size(pbd)-iend+1:size(pbd)-istart+1)=pbd(istart:iend)
  ptmp(size(pbd)-istart+2:size(pbd))=pbd(1:istart-1)
  thtmp(1:size(pbd)-iend)=thbd(iend+1:size(pbd))-pi
  thtmp(size(pbd)-iend+1:size(pbd)-istart+1)=thbd(istart:iend)
  thtmp(size(pbd)-istart+2:size(pbd))=thbd(1:istart-1)+pi
  pbd=ptmp
  thbd=thtmp
 end subroutine

 subroutine renewbd()
  integer::irenew
  pbd=pmin
  do irenew=1,nbd
   thbd(irenew)=pi*irenew/nbd
  end do
 end subroutine

 subroutine writepos()
  integer::iwrite
  open(30,file="pos.dat",status="replace")
  write(30,*) nparticle
  write(30,*) d2
  do iwrite=1,nparticle
   write(30,*) pth(iwrite), th(iwrite)
  end do
  close(30)
 end subroutine

 subroutine readpos()
  integer::iread
  open(30,file="posstart.dat",status="old")
  read(30,*) nparticle
  read(30,*) d2
  allocate(pth(nparticle),th(nparticle))
  do iread=1,nparticle
   read(30,*) pth(iread), th(iread)
  end do
  close(30)
  allocate(pbd(nbd),thbd(nbd))
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
  pr=2._8/(1-1/Pw**2)
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
