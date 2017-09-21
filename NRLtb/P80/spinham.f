      subroutine spinham (natoms,kkind,kbas,spinat,spinos,
     $     spndiag,spnmix,ltorque,dspth,dspphi,spndth,spnmxth,spndph)
c
c     Calculates the non-collinear spin corrections to the unpolarized
c      Hamiltonian matrix.  If ltorque is true, calculates derivatives
c      as well.
c
c=======================================================================
c
c     REVISION HISTORY
c
c     Created as part of static version 1.30
c                                                 -- mjm  10-15 Aug 2000
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
c     Input:
c
c     natoms is the number of atoms in the system
c
      integer natoms
c
c     Atom i is of type kkind(i)
c
      integer kkind(matom)
c
c     Atoms of type j have a basis set of size kbas(j)
c
      integer kbas(mkind)
c
c     Atom i is oriented in the direction spinat(l,i), l = 1,2,3 (x,y,z)
c      this is normalized to unity
c
      real*8 spinat(3,matom)
c
c     spinos(k,i) is the spin-orbit splitting for the ith atom, with
c      k = 1 -> s orbitals
c      k = 2 -> p orbitals
c      k = 3 -> t2g orbitals
c      k = 4 -> eg orbitals
c
      real*8 spinos(4,matom)
c
c     ltorque says that we should determine the derivatives of spndiag
c      and spnmix with respect to each angle.  By convention, we will
c      do this in degrees.
c
      logical ltorque
c
c     These are only used if ltorque is true:
c
c     dspth (i,iatom) = d spinat(i,iatom) / d theta(i), units energy/degree
c     dspphi(i,iatom) = d spinat(i,iatom) / d phi  (i), units energy/degree
c
      real*8 dspth(3,matom),dspphi(3,matom)
c
c-----------------------------------------------------------------------
c
c     Output:
c
c     spndiag(i) is the diagonal correction to the ith orbital of the
c      non-polarized Hamiltonian.  When constructing the final
c      Hamiltonian,
c        hmat(   i,   i) -> hmat(i,i) + spndiag(i) (spin up terms)
c        hmat(i+nv,i+nv) -> hmat(i,i) - spndiag(i) (spin down terms)
c
c     Note that mh will specify a doubled secular equation, so we
c      can dimension this to mh/mcol
c
      real*8 spndiag(mh/mcol)
c
c     spnmix(i) mixes the spin up and spin down calculation.  It is
c      complex, and we will have
c
c        hmat(i,i+nv) = spndiag(i)
c
c     with all other matrix elements (i,j+nv), 0<i,j<nv+1 set to zero.
c
      complex*16 spnmix(mh/mcol)
c
c     If ltorque is true, we also calculate
c
c     spndth (i) = d spndiag(i) / d theta(i)
c     spnmxth(i) = d spnmix (i) / d theta(i)
c     spndph (i) = d spnmix (i) / d phi  (i)
c
c     Note that the derivatives are calculated for ANGLES in DEGREES.
c
      real*8 spndth(mh/mcol)
      complex*16 spnmxth(mh/mcol), spndph(mh/mcol)
c
c=======================================================================
c
c     Integers
c
c     Do loop variables and other counters
c
      integer i,j,k
c
c     The "current" atom is of type kind
c
      integer kind
c
c     q is the current type of orbital (s,p,t2g,eg)
c
      integer q
c
c=======================================================================
c
c     i will represent the current global orbital, j the current atom, k
c      the current orbital for the current atom
c
      i = 0
      do j = 1,natoms
c
c        What kind of atom is this, and how big is its basis?
c
         kind = kkind(j)
c
ctemp
c$$$         write(15,*) 'Atom ',j,' kind ',kind,' basis',kbas(kind)
cend temp
         do k = 1,kbas(kind)
c
c           Update the global basis number:
c
            i = i + 1
c
c           Note that our selection of the orbital type depends on
c            the basis function
c
            if(k.eq.1) then
               q = 1
            else if (k.lt.5) then
               q = 2
            else if (k.lt.8) then
               q = 3
            else
               q  =4
            end if
            spndiag(i) = spinos(q,j)*spinat(3,j)
            spnmix(i)  = dcmplx(spinos(q,j)*spinat(1,j),
     $                          spinos(q,j)*spinat(2,j))
            if (ltorque) then
c
c              Update the angular derivatives if needed.  Note
c               that spndiag(i) is independent of phi(j).
c
               spndth(i)  = spinos(q,j)*dspth(3,j)
               spnmxth(i) = dcmplx(spinos(q,j)*dspth (1,j),
     $                             spinos(q,j)*dspth (2,j))
               spndph(i)  = dcmplx(spinos(q,j)*dspphi(1,j),
     $                             spinos(q,j)*dspphi(2,j))
            end if
         end do
      end do
ctemp
c$$$      write(15,*) 'Diagonal               Off-diagonal'
c$$$      do j = 1,i
c$$$         write(15,'(i5,f20.10,5x,2f20.10)') j,spndiag(j),spnmix(j)
c$$$      end do
cend temp
      return
      end
