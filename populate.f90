 program populate
 implicit none

 real(8)::pi=acos(-1._8),den
 real(8)::Pt,G,Pw,nh
 real(8)::pran,rparticle
 real(8),dimension(:),allocatable::ptot,qtot
 real(8)::pmax=1.0_8,pmin
 integer::ntot,nparticle=0,i

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
 Pw=1/Pw
 close(31)

 den=G/(2*pi*Pt)
! den=den*nhgen(0._8)/nh
 den=den*.02_8/nh
 ntot=idnint(den*(pmax-pmin)*pi)
 call initran()

 allocate(ptot(nparticle))
 do i=1,ntot
  call random_number(pran)
  pran=pmin+(pmax-pmin)*pran
  call random_number(rparticle)
!  if (rparticle<nhgen(sqrt(pran/Pw))/nhgen(0._8)) then
  if (rparticle<nhgen(sqrt(pran/Pw))/.02_8) then
   nparticle=nparticle+1
   call append()
  end if
 end do
 allocate(qtot(nparticle))
 call random_number(qtot)
 qtot=pi*qtot

 open(30,file="posread.dat",status="replace")
 write(30,*) nparticle
 do i=1,nparticle
  write(30,*) ptot(i),qtot(i)
 end do
 close(30)
 deallocate(ptot)
 deallocate(qtot)

 contains

 subroutine append()
  real(8),dimension(:),allocatable::ptemp
  allocate(ptemp(nparticle))
  ptemp(1:nparticle-1)=ptot
  ptemp(nparticle)=pran
  deallocate(ptot)
  allocate(ptot(nparticle))
  ptot=ptemp
  deallocate(ptemp)
 end subroutine 

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

 subroutine initran()
  integer::seedsize
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
