program parmtrconv
implicit none

 real(8)::rc2,rp,vr,n0,B,rw,nh
 real(8)::Pt,G,pmin
 real(8)::omegac,omegap,omegaE,pi=acos(-1._8)

 rw=3.5_8
 rc2=1.9_8
 rp=rc2/sqrt(2._8)
 vr=.8_8*.1
 nh=10._8**5
 n0=10._8**7
 B=1.2*10._8**4

 omegac=1.76*10.**7*B
 omegap=5.64*10.**4*sqrt(n0)
 omegaE=omegap**2/(2*omegac)

 Pt=(2*rc2*vr/rp**2)/omegaE
 G=.1_8**2
! G=5*.1_8**2
 pmin=0._8

 open(30,file="parameter.dat",status="replace")
 write(30,*) "Pt"
 write(30,*) Pt
 write(30,*) "G"
 write(30,*) G
 write(30,*) "pmin"
 write(30,*) pmin
 write(30,*) "delta=nhRw^2/n0Rp^2"
 write(30,*) nh/n0*(rw/rp)**2
 write(30,*) "nh/n0"
 write(30,*) nh/n0
 write(30,*) "Rp^2/Rw^2"
 write(30,*) (rp/rw)**2
 close(30)

end program
