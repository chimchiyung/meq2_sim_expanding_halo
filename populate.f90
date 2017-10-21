 program populate
 implicit none

 real(8)::pi=acos(-1._8),den
 real(8)::Pt,G
 real(8)::pran,qran
 real(8)::pmax=1._8,pmin
 integer::nparticle,i

 open(31,file="parameter.dat",status="old")
 read(31,*)
 read(31,*) Pt
 read(31,*)
 read(31,*) G
 read(31,*)
 read(31,*) pmin
 close(31)

 den=G/(2*pi*Pt)
 nparticle=idnint(den*(pmax-pmin)*pi)
 call initran()

 open(30,file="posread.dat",status="replace")
 write(30,*) nparticle
 do i=1,nparticle
  call random_number(pran)
  call random_number(qran)
  pran=pmin+(pmax-pmin)*pran
  qran=pi*qran
  write(30,*) pran,qran
 end do
 close(30)

 contains

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
