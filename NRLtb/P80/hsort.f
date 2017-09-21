      subroutine hsort(n,a,ind)
c
c @(#) hsort.f 1.1@(#)
c
c     Heap Sort (Faster than Metzner shell sort)
c
c     D. Singh (1986)
c
      implicit real*8(a-h,o-z)
      dimension a(*),ind(*)
      do j=1,n
         ind(j)=j
      end do
      if(n.le.1)return
      l=n/2+1
      ir=n
   10 continue
      if(l.gt.1)then
         l=l-1
         indxt=ind(l)
         q=a(indxt)
      else
         indxt=ind(ir)
         q=a(indxt)
         ind(ir)=ind(1)
         ir=ir-1
         if(ir.eq.1)then
            ind(1)=indxt
            return
         endif
      endif
      i=l
      j=l+l
 20   if(j.le.ir) then
         if(j.lt.ir) then
            if(a(ind(j)).lt.a(ind(j+1)))j=j+1
         endif
         if(q.lt.a(ind(j))) then
            ind(i)=ind(j)
            i=j
            j=j+j
         else
            j=ir+1
         endif
         go to 20
      endif
      ind(i)=indxt
      go to 10
      end
