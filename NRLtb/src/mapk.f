      subroutine mapk(rk,uk,inv)
      implicit real*8 (a-h,o-z)
c
C     ******************************************************************
C     *****
C     *****  MAPS CARTESIAN COORDINATES TO LATTICE COORDINATES AND
C     *****  THE REVERSE.
C     *****      RK(3)        K-POINT IN CARTESIAN COOR.
C     *****      UK(3)        K-POINT IN RECIP. LAT. COOR.
C     *****      INV= 1       CALCULATES UK FROM RK
C     *****      INV=-1       CALCULATES RK FROM UK
C     *****
C     *****              H. KRAKAUER  SEPT. 1984
C     ******************************************************************
c
      common /lattyp/ avec(3,3),bvec(3,3),bij(3,3),wsvol,lattic
c$$$      common /fundco/ pi,tpi,fpi
      parameter (pi = 3.14159265358979323846d0)
      parameter (tpi = 2d0*pi)
      parameter (tpiinv = 1d0/tpi)
      dimension rk(3),uk(3)
c
      if (inv .eq. 1) then
         do j=1,3
            uk(j)=(rk(1)*avec(1,j)+rk(2)*avec(2,j)+
     $           rk(3)*avec(3,j))*tpiinv
         end do
c
      else if (inv .eq.-1) then
         do j=1,3
            rk(j)= uk(1)*bvec(j,1)+uk(2)*bvec(j,2)+uk(3)*bvec(j,3)
         end do
c
      else
         stop 'MAPK:  INV .NE. TO 1 OR -1'
      end if
      return
      end
