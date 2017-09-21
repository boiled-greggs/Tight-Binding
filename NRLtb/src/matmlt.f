*DECK MATMLT
      SUBROUTINE MATMLT(A,B,C,N)
      IMPLICIT REAL*8 (A-H,O-Z)
C @(#) matmlt.f 1.1@(#)
C     **************************************************************
C     *****  MULTIPLIES REAL SQUARE  N X N  MATRICES  A X B AND
C     *****  PUTS THE RESULT IN MATRIX C
C     *****
C     *****             H. KRAKAUER
C     **************************************************************
      dimension a(n,n),b(n,n),c(n,n)
      parameter (zero = 0d0)
      do i=1,n
         do j=1,n
            c(j,i)=zero
         end do
      end do
      do k=1,n
         do j=1,n
            s=b(j,k)
            do i=1,n
               c(i,k)=a(i,j)*s+c(i,k)
            end do
         end do
      end do
      return
      end
