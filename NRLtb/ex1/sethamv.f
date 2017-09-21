      subroutine sethamv(jsover,ak,dparos,dhop,dover,dhmat,dsmat,lprint)
c
c     This code has been yanked from the original program setvol.f.
c     sethamv.f takes the perturbed onsite, Slater-Koster, and
c      (jsover being true) overlap matrix elements and the k-point
c      information from ak(i) and turns this into perturbed
c      Hamiltonina (dhmat) and overlap (dsmat) matrices by repeated
c      calls to rotate.f.
c
c=======================================================================
c
c     REVISION HISTORY for sethamv.f
c
c-----------------------------------------------------------------------
c
c     version 1.11:
c
c     Split off from setvol.f.  Put on-site, Slater-Koster, and
c      optional overlap parameters into the call statement.  Rewrite
c      the pars, parp, ... into an array dparos(type,atom).  Set setpar.f
c      for more details.
c
c     Rename ham, over, paros to dham, dover, and dparos, so that they
c      will not conflict with setpar's variables in bande.f.  This also
c      makes it a bit easier to distinguish the regular Hamiltonian
c      from the perturbation.
c
c     Get all these values from setparv.f as passed through bands.f.
c
c                                                     mjm -- 23 Jul 1999
c=======================================================================
c
c     REVISION HISTORY for setvol.f
c
c-----------------------------------------------------------------------
c
c     Changed to accomodate multi-atom systems.   -- mjm  29 Sept 1994
c
c-----------------------------------------------------------------------
c
c     Allows "asymptotic" form of the onsite parameters  -- mjm  6 Dec 1994
c
c-----------------------------------------------------------------------
c
c     Cubic polynomials for onsite terms (4 parameters each) -- mjm  14 Dec 1994
c
c-----------------------------------------------------------------------
c
c     On-site d terms split, quadratic Slater-Koster prefactors -- mjm 31 Oct 95
c
c-----------------------------------------------------------------------
c
c     Replace (dhr,dhi), (dsr,dsi) by dhmat and dsmat -- mjm  5 Feb 1997
c
c-----------------------------------------------------------------------
c
c     Changed on-site form to match that of setup.f, where different
c      atom types contribute differently to the on-site terms through
c      their densities.                               -- mjm 15 Apr 1997
c
c-----------------------------------------------------------------------
c
c     search2 has been changed to eliminate the double counting of pairs.
c     This means that we'll have to make sure that the density
c      terms are incremented properly.
c                                                mjm --  2 May  1997
c
c-----------------------------------------------------------------------
c
c     Eliminated references to "asymptotic" on-site terms
c                                                mjm -- 23 Sept 1997
c
c-----------------------------------------------------------------------
c
c     If the parameter realov is true, then like atom overlap matricies
c      are constrained to have the correct behavior as R -> 0, i.e.
c      S_{l,l',m}(R) -> delta_{l,l'}
c
c     For these terms (and only these terms), the parametrization
c      is of the form
c
c     S(R) = [delta_{l,l'} + R( A + B R + C R^2)] Exp[-D^2 R]
c
c                                                mjm -- 25 Sept. 1997
c
c-----------------------------------------------------------------------
c
c     Version 1.05:
c
c     As part of the general upgrade, jsover is moved into the calling
c      parameters.  Remember that jsover = 1 for non-orthogonal
c      calculations, and 0 for orthogonal calculations.
c                                                mjm -- 20 April 1998
c
c-----------------------------------------------------------------------
c
c     Removed the "labels" common block, and pass the volume ("thisvol")
c      via the calling parameters.
c                                                mjm -- 22 April 1998
c
c-----------------------------------------------------------------------
c
c     Version 1.06:
c
c     jsover has been changed to a logical variable:
c        jsover = .true.  -> Non-orthogonal Hamiltonian (S <> identity)
c        jsover = .false. -> Orthogonal Hamiltonian (S = identity)
c
c                                                mjm --  6 July  1998
c
c     Change the name of the common block /overlap/ to /parcom/.  It now
c      includes the parametrization type variable nltype.
c     When nltype = 90000 use a special extended polynomial
c      representation for H_{ss sigma} and S_{ss sigma}.  This is
c      restricted to the H atom, so set all other Slater-Koster
c      parameters to 1.
c                                                mjm --  4 Aug   1998
c=======================================================================
c
      implicit none
      include 'P1'
c
c=======================================================================
c
c     Arguments
c
c     jsover is true if this is a non-orthogonal calculation
c
      logical jsover
c
c     ak(i) is the i-th k-point coordinate
c
      real*8 ak(3)
c
c     dparos(k,iat) contains the perturbed value of the onsite term
c       for the k-th orbital type (see below) for atom iat
c
      real*8 dparos(4,matom)
c
c     dhop(i,j) is the ith Slater-Koster perturbed matrix element
c      of the jth pair.  dover(i,j) is the corresponding overlap
c      element
c
      real*8 dhop(mpkind,mpair),dover(mpkind,mpair)
c
c     dhmat and dsmat are the packed complex Hamiltonian and overlap
c      matrices
c
      complex*16 dhmat(mh*(mh+1)/2),dsmat(mh*(mh+1)/2)
c
c     If lprint is true, print diagnostic information
c
      logical lprint
c
c=======================================================================
c
c     Real Arrays
c
c     vpar and vovl hold the Slater-Koster parameters in a form
c      understood by rotate.f
c
      real*8 vpar(mpkind,mkind,mkind),vovl(mpkind,mkind,mkind)
c
c     tt(i) is the local value of tt_list(i,ipair)
c
      real*8 tt(3)
c
c-----------------------------------------------------------------------
c
c     Real scalars
c
c-----------------------------------------------------------------------
c
c     Misc. integers.
c
      integer jk1,jk2,i,iat,j,kb,kk,kkj,l1,l2,mij,npk,nv
c
c-----------------------------------------------------------------------
c
c     Counter names:
c
      integer ipair
c
c-----------------------------------------------------------------------
c
c     Common blocks:
c
      real*8 tt_list,dlv_list,dist_list,screen_list
      integer jkind_list,jatm_list,npairs
c
      common /pairs/tt_list(3,mpair),
     $     dlv_list(3,mpair),dist_list(mpair),
     $     screen_list(mpair),
     $     jkind_list(mpair,2),
     $     jatm_list(mpair,2),
     $     npairs
c
      real*8 posn(matom,3)
      integer kkind(matom),natoms
      common /codes/posn,kkind,natoms
c
      real*8 valence(mkind)
      integer kinds,kbas(mkind)
      common /parinfo/ valence,kinds,kbas
c
      integer npkind(mptype)
      common /types/ npkind
c
      real*8 dlv
      integer jkind,jatm,kneigh
      common /relat/dlv(3),jkind(2),jatm(2),kneigh
c
      integer ksk(matom)
      common /kstuff/ ksk
c
c     Only needed when we are printing the Hamiltonian.  Otherwise
c      it can be commented out.
c
      common/struc/ nv
c
c=======================================================================
c
c     Parameters
c
      real*8 zero
      parameter (zero  = 0d0)
c
c     Tolerance for printing Hamiltonian matrix elements
c
      real*8 tol
      parameter (tol=1d-10)
c
c=======================================================================
c
c
c     Initialize the Perturbed Hamiltonian and Overlap matrices
c
      do mij=1,mh*(mh+1)/2
         dhmat(mij) = dcmplx(zero,zero)
         dsmat(mij) = dcmplx(zero,zero)
      end do
      kneigh=1
      do ipair=1,npairs
         jkind(1)=jkind_list(ipair,1)
         jkind(2)=jkind_list(ipair,2)
         jatm(1) =jatm_list(ipair,1)
         jatm(2) =jatm_list(ipair,2)
c
c        Properly order the atom types:
c
         jk1=min(jkind(1),jkind(2))
         jk2=max(jkind(1),jkind(2))
c
         if(jk2.eq.jk1) then
            npk = npkind(2)
         else
            npk = npkind(3)
         end if
         do j = 1,npk
            vpar(j,jk2,jk1) = dhop(j,ipair)
         end do
         if(jk1.ne.jk2)then
            do j=1,npkind(2)
               vpar(j,jk1,jk2)=vpar(j,jk2,jk1)
            end do
            vpar(2,jk1,jk2)=-vpar(11,jk2,jk1)
            vpar(5,jk1,jk2)= vpar(12,jk2,jk1)
            vpar(6,jk1,jk2)=-vpar(13,jk2,jk1)
            vpar(7,jk1,jk2)=-vpar(14,jk2,jk1)
         endif
         vpar(11,jk1,jk2)=-vpar(2,jk2,jk1)
         vpar(12,jk1,jk2)= vpar(5,jk2,jk1)
         vpar(13,jk1,jk2)=-vpar(6,jk2,jk1)
         vpar(14,jk1,jk2)=-vpar(7,jk2,jk1)
c
c        Overlap matrices
c
         if(jsover)then
c
            if(jk2.eq.jk1) then
               npk = npkind(2)
            else
               npk = npkind(3)
            end if
            do j = 1,npk
               vovl(j,jk2,jk1) = dover(j,ipair)
            end do
            if(jk1.ne.jk2)then
               do j=1,npkind(2)
                  vovl(j,jk1,jk2)=vovl(j,jk2,jk1)
               end do
               vovl(2,jk1,jk2)=-vovl(11,jk2,jk1)
               vovl(5,jk1,jk2)= vovl(12,jk2,jk1)
               vovl(6,jk1,jk2)=-vovl(13,jk2,jk1)
               vovl(7,jk1,jk2)=-vovl(14,jk2,jk1)
            endif
            vovl(11,jk1,jk2)=-vovl(2,jk2,jk1)
            vovl(12,jk1,jk2)= vovl(5,jk2,jk1)
            vovl(13,jk1,jk2)=-vovl(6,jk2,jk1)
            vovl(14,jk1,jk2)=-vovl(7,jk2,jk1)
         endif
c
c        Get the vector representing the pair separation
c         and orientation
c
         do j=1,3
            tt(j)=tt_list(j,ipair)
            dlv(j)=dlv_list(j,ipair)
         end do
c
c        Set up perturbation matricies
c
         call rotate(tt(1),tt(2),tt(3),ak,vpar,vovl,
     $        dhmat,jsover,dsmat,ksk,kbas)
      end do
c
c     On-site terms
c
      do iat=1,natoms
         jkind(1)=kkind(iat)
         jkind(2)=kkind(iat)
         jatm(1)=iat
         jatm(2)=iat
c
c        Directly include the onsite part of rotate here,
c         rather than calling it.  Note that only the real
c         part of dH is affected, since the onsite terms are real.
c         dS is not affected because the onsite overlap matrix is
c         the Identity and its volume derivative is zero.
c
         kk=ksk(jatm(1))
         kb = kbas(jkind(1))
         if(kb.gt.0) then
c
c           S states
c
            kkj = kk + 1
            dhmat(kkj*(kkj+1)/2) = dhmat(kkj*(kkj+1)/2) +
     $           dcmplx(dparos(1,iat),zero)
c
c           P states
c
            do j = 2,min(4,kb)
               kkj = kk + j
               dhmat(kkj*(kkj+1)/2) = dhmat(kkj*(kkj+1)/2) +
     $              dcmplx(dparos(2,iat),zero)
            end do
c
c           D-t2g states
c
            do j = 5,min(7,kb)
               kkj = kk + j
               dhmat(kkj*(kkj+1)/2) = dhmat(kkj*(kkj+1)/2) +
     $              dcmplx(dparos(3,iat),zero)
            end do
c
c           D-eg states
c
            do j = 8,min(9,kb)
               kkj = kk + j
               dhmat(kkj*(kkj+1)/2) = dhmat(kkj*(kkj+1)/2) +
     $              dcmplx(dparos(4,iat),zero)
            end do
         end if
      end do
C.......................................................
C--->  PRINTOUT OF HAMILTONIAN AND OVERLAP MATRICES
C.......................................................
      if(lprint) then
         write(15,*) 'Secular equation dimension = ',nv
         write(15,'(/''Hamiltonian'')')
         i = 0
         do l1 = 1,nv
            do l2 = 1,l1
               i = i + 1
               if(cdabs(dhmat(i)).gt.tol)
     $              write(15,51) i,l1,l2,dhmat(i)
 51            format(3i6,2f20.10)
            end do
         end do
         if(jsover) then
            write(15,'(/''Overlap Matrix'')')
            i = 0
            do l1 = 1,nv
               do l2 = 1,l1
                  i = i + 1
                  if(cdabs(dsmat(i)).gt.tol)
     $                 write(15,51) i,l1,l2,dsmat(i)
               end do
            end do
         end if
      end if
      return
      end
