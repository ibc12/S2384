      character*60 nom,nom2
40    format(a60)
      read 40, nom
      read 40, nom2
   
      read *, fac
      open(7,file=nom,status='old')
      open(11,file=nom2,status='unknown')
      do i=1,100000
        read(7,*,err=99,end=99)ang,sig
c        print*,ang,sig, sig*fac
        write(11,*),ang,sig*fac
      enddo
99    continue
      print 40,nom2
      end
      
