program rk6const
 implicit none
 real(8),dimension(7,7)::b
 real(8),dimension(7)::a
 integer::i,j

 a=0._8
 b=0._8
 b(1,2)=1._8/3
 b(1,3)=0._8
 b(2,3)=2._8/3
 b(1,4)=1._8/12
 b(2,4)=1._8/3
 b(3,4)=-1._8/12
 b(1,5)=-1._8/16
 b(2,5)=9._8/8
 b(3,5)=-3._8/16
 b(4,5)=-3._8/8
 b(1,6)=0._8
 b(2,6)=9._8/8
 b(3,6)=-3._8/8
 b(4,6)=-3._8/4
 b(5,6)=1._8/2
 b(1,7)=9._8/44
 b(2,7)=-9._8/11
 b(3,7)=63._8/44
 b(4,7)=18._8/11
 b(5,7)=0._8
 b(6,7)=-16._8/11
 a(1)=11._8/120
 a(2)=0._8
 a(3)=27._8/40
 a(4)=27._8/40
 a(5)=-4._8/15
 a(6)=-4._8/15
 a(7)=11._8/120

 open(30,file="rk6consts.dat")
 write(30,*) a(:)
 do i=1,7
  write(30,*) b(i,:)
 end do
 close(30)

end program