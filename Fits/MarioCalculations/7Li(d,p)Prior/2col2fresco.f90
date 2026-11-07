      program transpot
      
      implicit none
      
      
      integer i,j,nr
      real*8 impot(100000),rpot(100000),r(100000),dr,rfirst
      character*40 header
      logical ifcmplx
      
      write(0,*) 'Put T on top of file if potential is complex or F if it is real'
      ifcmplx=.false.
      read(5,*) ifcmplx
      
      read(5,*) header
      write(0,*) 'Header: ',header
      do i=1,100000
      if (ifcmplx) then
      read(5,*,end=999) r(i),rpot(i),impot(i)
      else
      read(5,*,end=999) r(i),rpot(i)
      endif
      enddo
999   continue

      nr=i-1
      dr=r(3)-r(2)
      rfirst=r(1)
      write(6,*) '#Expressed as fresco ',header
      write(6,*) nr,dr,rfirst
      if (ifcmplx) then
      write(6,'(1P,6E12.4)') (rpot(i),impot(i),i=1,nr)   
      else
      write(6,'(1P,6E12.4)') (rpot(i),i=1,nr)   
      endif
      
      end program 
