      real*8 FUNCTION GDOT(A,B,GIJ)
      IMPLICIT REAL*8 (A-H,O-Z)
C     ***********************************************************
C     *****
C     *****   DOT PRODUCT OF TWO VECTORS WITH METRIC GIJ
C     *****
C     ***********************************************************
c
      real*8 A(3),B(3),GIJ(3,3)
      GDOT=0.0D0
c$$$      DO 1 J=1,3
c$$$      DO 1 I=1,3
c$$$    1 GDOT=GDOT+A(I)*GIJ(I,J)*B(J)
      do j = 1,3
         sum = 0d0
         do i = 1,3
            sum = sum + a(i)*gij(i,j)
         end do
         gdot = gdot + sum*b(j)
      end do
      RETURN
      END
