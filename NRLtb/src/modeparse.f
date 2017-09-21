      subroutine modeparse(string,number)
c
c     Converts the two-digit string "string" to the integer "number",
c      treating all non-numeric symbols, including "+" and "-", as
c      nulls.
c
      implicit none
c
      character*2 string
      integer number,i(2),j,k
c
      character*1 char(0:9)
c
c-----------------------------------------------------------------------
c
      data char /'0','1','2','3','4','5','6','7','8','9'/
c
c-----------------------------------------------------------------------
c
      do j = 1,2
         i(j) = -1
         do k = 0,9
            if(string(j:j).eq.char(k)) i(j) = k
         end do
      end do
c
      if(i(1).lt.0) i(1) = 0
      if(i(2).lt.0) then
         number = i(1)
      else
         number = 10*i(1)+i(2)
      end if
      return
      end
