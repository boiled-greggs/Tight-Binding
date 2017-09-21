      subroutine addspins(nv,spndiag,spnmix,hmat,smat,lprint)
c
c     Mixes spin up and spin down vectors by doubling the size of the
c      secular equation from nv to 2*nv, adding spndiag to the spin-up
c      diagonal term of hmat, subtracting it from the spin down part,
c      and using spnmix to mix the two parts of H (but not S).
c
c=======================================================================
c
c     REVISION HISTORY
c
c     Created as part of static version 1.30
c                                                    -- mjm  10 Aug 2000
c=======================================================================
c
      implicit none
c
c=======================================================================
c
      include 'P1'
c
c=======================================================================
c
c     External variables:
c
c-----------------------------------------------------------------------
c
c     nv is the size of the basis set, without spin
c
      integer nv
c
c     spndiag(i) is the diagonal correction to the ith orbital due
c      to spin.  It is added to spin up, and subtracted from spin down
c
      real*8 spndiag(mh/mcol)
c
c     spnmix(i) mixes orbital i and orbital nv+i
c
      complex*16 spnmix(mh/mcol)
c
c     hmat and smat are the Hamiltonian and Overlap matrices for the
c      current k-point, stored in triangular form.  On input they hold
c      the average-spin matrix elements, from element (1,1) to (nv,nv).
c      On output they contain the full spin-corrected matrix, from
c      (1,1) to (2*nv,2*nv).
c
      complex*16 hmat(mh*(mh+1)/2),smat(mh*(mh+1)/2)
c
c     If lprint is true, non-zero matrix elements will be printed
c
      logical lprint
c=======================================================================
c
c     Integers
c
c     Do loop variables:
c
      integer i,j,k
c
c=======================================================================
c
c     Constants
c
      real*8 zero,tol
c
      parameter (zero = 0d0)
      parameter (tol = 1d-10)
c
c=======================================================================
c
c     In triangular storage format, if j > i, then matrix element
c      H(j,i) is stored in hmat(j*(j+1)/2+i)
c
      do j = 1,nv
c
c        Off-diagonal corrections:
c
         do i = 1,j-1
c
c           Duplicate the spin-up block in the spin-down block:
c
            hmat((j+nv)*(j+nv-1)/2+i+nv) = hmat(j*(j-1)/2+i)
            smat((j+nv)*(j+nv-1)/2+i+nv) = smat(j*(j-1)/2+i)
c
c           Zero this element in the spin-up/spin-down mixing block:
c
            hmat((j+nv)*(j+nv-1)/2+i) = dcmplx(zero,zero)
            smat((j+nv)*(j+nv-1)/2+i) = dcmplx(zero,zero)
         end do
c
c        Diagonal corrections:
c
c        Spin-down:
c
         hmat((j+nv)*(j+nv+1)/2) = hmat(j*(j+1)/2) -
     $        dcmplx(spndiag(j),zero)
         smat((j+nv)*(j+nv+1)/2) = smat(j*(j+1)/2)
c
c        Spin-up:
c
         hmat(j*(j+1)/2) = hmat(j*(j+1)/2) + dcmplx(spndiag(j),zero)
c
c        Cross terms:
c
         hmat((j+nv)*(j+nv-1)/2+j) = spnmix(j)
         smat((j+nv)*(j+nv-1)/2+j) = dcmplx(zero,zero)
c
      end do
c
c     Printing options:
c
      if(lprint) then
         k = 0
         do j = 1,2*nv
            do i = 1,j
               k = k + 1
               if(abs(dreal(hmat(k)))+abs(dimag(hmat(k)))+
     $            abs(dreal(smat(k)))+abs(dimag(smat(k))).gt.tol)
     $              write(15,'(i5,i5,5x,2f15.10,5x,2f15.10)')
     $              j,i,hmat(k),smat(k)
            end do
         end do
      end if
c
      return
      end
