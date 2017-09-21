      subroutine diag(mode,lnoncol,hmat,smat,en,
     $     lspdir,spev,n,zeig,jsover,
     $     ak,eigcut,weight,ioverr,natoms,kkind,kbas,ksk)
c
c     Arguments after ioverr are needed for a proper angular momentum
c      decomposition of an atom.
c
c=======================================================================
c
c     REVISION HISTORY:
c
c-----------------------------------------------------------------------
c
c     New version.  Set up to call LAPACK, which we'll let
c      solve for eigenvalues and/or eigenvectors
c
c-----------------------------------------------------------------------
c
c     If lqlmt = .true., calculate partial occupancies and
c      add to the QLMT file (for eigenvalues < eigcut)
c
c                                        --mjm  14 July 1994
c
c-----------------------------------------------------------------------
c
c     If lqapw = .true., print out the eigenvalues in
c      APW format.                       --mjm  12 Sept 1996
c
c-----------------------------------------------------------------------
c
c     Finally rewrite the code so that (hr,hi) and (sr,si) are
c      stored as complex numbers.  Also, note that hmat is used once
c      and not needed again, so we don't need to copy from (hr,hi) to
c      hmat.                             --mjm   5 Feb  1997
c
c-----------------------------------------------------------------------
c
c     Added "lqlmt0" option, which allows printing of a QLMT-like
c      file, with complete k-point information, with no "atomic"
c      information, i.e., angular momentum decomposition.  Note that
c      unlike lqlmt, this option allows mpress=0 and meigen=0 in the
c      P1 file.
c                                        --mjm  28 Oct  1997
c
c-----------------------------------------------------------------------
c
c     Output for the 'lq*' modes is now written to a binary file
c      on unit 55 instead of directly into ASCII form on unit 45.
c      This makes it easier to run the code in parallel.
c      The subroutine bande.f will concatinate the files produced
c      by different processors in multi-processor mode.
c
c     We've also added an 'ioverr' variable to the output.  If it is
c      non-zero it alerts the program to the fact that the overlap
c      matrix is bad.
c
c                                        --mjm  27 Feb  1998
c
c-----------------------------------------------------------------------
c
c     Fixed a bug in the angular momentum decomposition of the density
c      of states:  smat, the overlap matrix, is destroyed by zhpgv in
c      the determination of the eigenstates of the matrix.  On the
c      other hand, we need smat to do the Lowdin decomposition to get
c      the angular momentum decomposition.  To do this, we make a copy
c      smat1 == smat.  The bug was that we then used smat1 for the
c      original call to zhpgv, and then used the same copy for the
c      Lowdin decomposition.  This lead to unfortunate results, to say
c      the least.  EIGENVALUES AND TOTAL ENERGIES WERE NOT AFFECTED.
c      The bug has been fixed.  smat is destroyed where it should be
c      and smat1 is used only where it should be.
c                                        --mjm  20 May  1998
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
c-----------------------------------------------------------------------
c
c     Separated the eigenvector calculation parameters into a parameter
c      file "P2".
c
c     Moved k-point generation parameters into 'P3'
c
c                                                mjm --  6 July  1998
c-----------------------------------------------------------------------
c     
c     Speed up calculations for 1X1 matrix       mjm --  4 Aug   1998
c-----------------------------------------------------------------------
c
c     Code modified to print a nonsense angular momentum decomposition
c      if the Cholesky decomposition of the overlap matrix fails.  This
c      allows the complete printing of the QLMT file, highlighting where
c      the overlap matrix goes wrong.            mjm -- 27 Aug   1998
c-----------------------------------------------------------------------
c
c     Oops!  The nonsense decomposition was only done for the first
c      atom in the unit cell.  Fixed that.       mjm -- 17 Sept  1998
c-----------------------------------------------------------------------
c
c     Replaced the ever-increasing number of lqapw, lqlmt, etc. with
c      "mode".  Modes recognized by diag.f are
c
c     1,3 -- Only calculates the eigenvalues
c     2,4 -- Calculates eigenvalues and perturbation
c     5,7 -- Unused
c     6   -- Calculates eigenvalues and angular momentum decomposition
c     6,8,9 -- Dumps appropriate eigenvalue/decomposition information
c              into file 55:
c                                                    -- mjm  27 Jan 1999
c-----------------------------------------------------------------------
c
c     For the purposes of this subroutine, mode 5 is now identical to
c      mode 6.  To make it easier to use this fact, we've defined a
c      logical variable "mode56", true if mode=5 or mode=6.
c                                                    -- mjm  11 Feb 1999
c-----------------------------------------------------------------------
c
c     Fixed angular momentum decomposition of atoms in the case when
c      the number of basis points per atom was smaller than the mbas,
c      as well as changing the basis size for different kinds of atoms,
c      as in NbC.  This requires passing extra information into the
c      program.
c                                                    -- mjm  22 Jul 1999
c-----------------------------------------------------------------------
c
c     Version 1.30:
c
c     Fix the QLMT file angular momentum decomposition section to
c      handle the case of non-collinear spins,  properly decompose the
c      wave function along the global spin axis.
c
c     To do this we need to know that we are doing a non-collinear
c      calculation, hence "lnoncol" has been added to the calling
c      statement.
c                                                    -- mjm  14 Aug 2000
c-----------------------------------------------------------------------
c
c     Version 1.31:
c
c     If lspdir is true, we calculate the expectation value of the
c      spin operators (sigma_{x,y,z}) for each atom and eigenstate
c      these are stored in the array spev, which has dimensions
c
c     spev({x,y,z},atom,eigenstate,k-point)
c
c     Note that the k-point variable is not used here, so we've
c      used the old Fortran trick and eliminated it from the
c      dimensions.  As such, calls to diag should reference
c      spev as spev(1,1,1,kpt), where kpt is the current k-point
c      under discussion.
c
c     ioverr now does double duty.  If ioverr = -1, then we have
c      a logical inconsistency.
c                                                    -- mjm  29 Nov 2000
c-----------------------------------------------------------------------
c
c     Version 1.32:
c
c     Mode 10 prints out the occupation of every orbital in the system.
c      The variable mode56 will now also be true if mode=10.
c                                                    -- mjm   1 Nov 2001
c=======================================================================
c
      implicit real*8 (a-h,o-z)
c
      include 'P1'
      include 'P2'
c
c     lnoncol is true if the calculation is for a system with
c      possible non-collinear magnetization.  Note that this is only
c      important when doing a angular-momentum-and-spin decomposition
c      of the eigenvalues for a QLMT file.
c
      logical lnoncol
c
      logical jsover
c
c     leigen is true if we need the eigenvectors:
c
      logical leigen
c
c     lspdir is true if we are calculating spin directions:
c
      logical lspdir
c
c     spev gives the direction of each spin on each atom for
c      each eigenstate:
c
      real*8 spev(3,msdir*(matom-1)+1,msdir*(mh-1)+1)
c
c     orbocc(i) is the occupation of the ith orbital on the current
c      atom for the current eigenstate.
c
      real*8 orbocc(mbas)
c
c-----------------------------------------------------------------------
c
      complex*16 hmat(mh*(mh+1)/2),smat(mh*(mh+1)/2)
      real*8 en(mh)
      real*8 ak(3)
c
c     We need this to print the k-points in the QLMT and QAPW case:
c
c     Arrays needed for the LApack routines
c
c     smat1 is (maybe) a copy of smat used when we need to decompose
c      the eigenvalues.  Note that it is only properly defined
c      when meigen=1
c
      complex*16 smat1(mhe*(mhe+1)/2)
      complex*16 zwork(2*mh-1)
c$$$      complex*16 dcmplx,dconjg
      real*8 rwork(3*mh-2)
c
c     Dimension zeig to (mh,mh) if you want the eigenvectors.
c      Note that we must dimension it so in this version,
c      since we want to be able to find the eigenvalues
c      If you're doing a lot of runs in a large system, but
c      don't want the QLMT file, you might neglect this.  Note
c      that the P2 file has a special dimension for these vectors.
c
      complex*16 zeig(mhe,mhe)
c
c     Unfortunately, to do the partial occupancy writing we'll also
c      need the eigenvectors of s
c
      complex*16 svect(mhe,mhe),u(mhe,mhe)
c
c     We use seig elsewhere, so we'll leave this dimension alone
c
      real*8 seig(mh)
      complex*16 csum
c
c     mode56 is true if the mode is 5, 6, or 10
c
      logical mode56
c
c-----------------------------------------------------------------------
c
c     Arguments needed for angular momentum decomposition:
c
c     kkind(i) gives the type of the ith atom in the basis
c     kbas(j) gives the number of basis functions for atoms of type j,
c      i.e., atom i has kbas(kkind(i)) basis functions.
c
c     Just to make counting easier, ksk(i) tells the system
c      where the basis functions for atom i start (in the top half
c      of the Hamiltonian, if we have a non-collinear or spin-orbit
c      system).
c
      integer kkind(matom),kbas(mkind),ksk(matom)
c
c=======================================================================
c
      parameter (zero = 0d0,one = 1d0)
c
c=======================================================================
ctemp
c$$$      write(*,*) 'In diag'
cend temp
c
c     Consistency check.  If lspdir is true, lnoncol must also be true.
c
      if(lspdir) then
         if(.not.lnoncol) then
            write(*,*) 'lspdir cannot be true if lnoncol is false'
            ioverr = -1
            return
         end if
c
c        In addition, the secular equation dimension should be twice
c         the number of available real-space basis functions:
c
         if(n.ne.(2*(ksk(natoms)+kbas(kkind(natoms))))) then
            write(*,*) 'Improper basis size for lspdir = true'
            ioverr = -2
            return
         end if
c
c        We'll need this later:
c
         nhalf = n/2
      end if
c
c     mode 5 and 6 are identical in this program.
c     mode 10 is close
c
      mode56=(mode.eq.5).or.(mode.eq.6).or.(mode.eq.10)
c
c     If we are decomposing a non-collinear magnetic system
c      then there are twice as many basis functions as in
c      the collinear case (where spin down is handled
c      separately from spin up).  The spin-down orbitals are below
c      all of the spin-up orbitals.  We'll handle this by doubling
c      what the code thinks of as the number of atoms in the system.
c
      if(lnoncol) then
         natx = 2*natoms
      else
         natx = natoms
      end if
c     
c     Do we need the eigenvectors?
c
      leigen=(mode.eq.2).or.(mode.eq.4).or.mode56.or.lspdir
ctemp
c$$$      write(*,*) 'leigen = ',leigen
cend temp
c
      ioverr = 0
c
c     Special case for a 1X1 matrix:
c
      if(n.eq.1) then
         if(jsover) then
            h1 = dble(hmat(1))
            s1 = dble(smat(1))
c
c           Non-orthogonal.  Need on
c
            if(s1.gt.zero) then
               en(1) = h1/s1
            else
c
c              Negative Overlap Matrix Element.  This is
c               forbidden.
c
c              Create an eigenvalue anyway.  Make it rather large:
c
               en(1) = abs(h1/s1)
c
c              Send error message
c
               ioverr = 1
               write(*,*) 'Found s(1) = ',s1,' use ',en(1)
               write(0,*) 'Found s(1) = ',s1,' use ',en(1)
            end if
         else
c
c           Orthogonal:
c
            en(1) = dble(hmat(1))
         end if
c
c        Special cases:
c
         if(mode56) then
            ntot = 1
            write(55) ak(1),ak(2),ak(3),natx,ntot,weight
            write(55) en(1)
            write(55) one,zero,zero
         else if(mode.eq.8) then
            write(55) n,ak(1),ak(2),ak(3),en(1)
         else if(mode.eq.9) then
            write(55) ak(1),ak(2),ak(3),n,weight,en(1)
         end if
         return
      end if
c
c     Pack the arrays:
c
      if(.not.jsover) then
         do m = 1,n*(n+1)/2
            smat(m) = dcmplx(zero,zero)
         end do
      else if(meigen.eq.1) then
c
c        We do want to keep a copy of smat around for use in the
c         angular momentum decomposition phase
c
         do m = 1,n*(n+1)/2
            smat1(m) = smat(m)
         end do
      end if
c
c     To calculate eigenvectors, redimension zeig and change the 'N'
c      to 'V'
c
      if(leigen) then
c
c        We need the eigenvectors for QLMT or pressure calculations
c
         if(jsover) then
c
c           Non-orthogonal calculation:
c
ctemp
c$$$            write(*,*) 'Getting eigenvectors (leigen,jsover)'
cend temp
            call zhpgv(1,'V','U',n,hmat,smat,en,zeig,mhe,
     $           zwork,rwork,info)
         else
c
c           Orthogonal calculation:  don't need the overlap matrix
c
            call zhpev('V','U',n,hmat,en,zeig,mhe,zwork,rwork,info)
         end if
      else
         if(jsover) then
ctemp
c$$$            write(*,*) 'Getting eigenvectors (!leigen,jsover)'
cend temp
            call zhpgv(1,'N','U',n,hmat,smat,en,zeig,mhe,
     $           zwork,rwork,info)
         else
            call zhpev('N','U',n,hmat,en,zeig,mhe,zwork,rwork,info)
         end if
      end if
ctemp
c$$$      write(*,*) 'Got eigenvectors'
cend temp
ctemp
c$$$      write(*,'(''K = '',3f10.4)') (ak(idir),idir=1,3)
c$$$      write(*,*) 'info = ',info
c$$$      write(*,'(8f10.4)') (en(ieig),ieig=1,n)
c$$$      write(15,'(''K = '',3f10.4)') (ak(idir),idir=1,3)
c$$$      write(15,'(8f10.4)') (en(ieig),ieig=1,n)
cend temp
c
c     A non-zero info means an error, probably smat is not positive
c      definite.
c
      if(info.ne.0) then
         write(*,*) 'Error return from zhpgv:  info = ',info
         write(*,*) 'Probably a negative eigenvalue in S'
c
         if(meigen.eq.1) then
c
c           Since we have a copy of smat, might as well see how
c            bad things are going to get:
c
c$$$            write(0,*) 'Overlap matrix at k = ',ak
c$$$            i = 0
c$$$            do l1 = 1,n
c$$$               do l2 = 1,l1
c$$$                  i = i + 1
c$$$                  write(0,'(3i5,2f15.8)') l1,l2,i,smat(i)
c$$$               end do
c$$$            end do
            call zhpev('N','U',n,smat1,en,zeig,mhe,zwork,rwork,info)
            write(*,*) 'The eigenvalues of S are:  '
            write(*,'(i5,f15.8)') (j,en(j),j=1,n)
         end if
         ioverr = 1
c
c        Don't return yet.  We'll just note that ioverr is set.
c
c$$$         return
      end if
c
c
      if(lspdir) then
c
c        Calculate the expectation value of the spin operator for
c         each eigenstate for the current (unspecified) k-point
c
c        I {\em think} that this can be done without doing a
c         decomposition of the eigenvector zeig.  However, time
c         will tell.
c
c        Note that this only works when lnoncol is true, i.e.,
c         both sets of spins are in the Hamiltonian, and the
c         eigenvectors must be calculated (leigen = .true.)
c         These conditions have been checked/set above.
c
c        Scan over eigenvalues
c
         do ieig = 1,n
c
c           Scan over atoms
c
            do iat = 1,natoms
c
c              Zero the directions:
c
               spinx = zero
               spiny = zero
               spinz = zero
c
c              Sum over the (up,down) basis functions
c
               do ibas = ksk(iat)+1,ksk(iat)+kbas(kkind(matom))
c
c                 Note that spinx is real:
c
                  spinx = spinx
     $                 - dconjg(zeig(ibas,ieig))*
     $                 zeig(ibas+nhalf,ieig)
     $                 - dconjg(zeig(ibas+nhalf,ieig))
     $                 *zeig(ibas,ieig)
c
c                 Check the sign of spiny:
c
                  spiny = spiny + dimag(
     $                 dconjg(zeig(ibas,ieig))*
     $                 zeig(ibas+nhalf,ieig)
     $                 - dconjg(zeig(ibas+nhalf,ieig))*
     $                 zeig(ibas,ieig)   )
c
c                 spinz is again real:
c
                  spinz = spinz +
     $                 dconjg(zeig(ibas+nhalf,ieig))
     $                 * zeig(ibas+nhalf,ieig)
     $                 - dconjg(zeig(ibas,ieig))*zeig(ibas,ieig)
               end do
c
c              Assign the spins to the proper matrix element:
c
               spev(1,iat,ieig) = spinx
               spev(2,iat,ieig) = spiny
               spev(3,iat,ieig) = spinz
ctemp
c$$$               write(*,*) ieig,iat,spinx,spiny,spinz
cend temp
            end do
         end do
      end if
c
      if(mode56) then
c
c        If we have a non-orthogonal basis we've got some work to do.
c         I think this is equivalent to Lowdin's method.
c
         if(jsover.and.(ioverr.eq.0)) then
c
c           In this case we have a successful Cholesky decomposition
c            of smat, i.e., all of its eigenvalues are positive
c            definite.  If the decomposition fails, we'll print out
c            a set of defining junk in the QLMT file, but we don't
c            need to do any decomposition.
c
c           Used smat1, where we stored smat.  We'll need
c            to diagonalize it to get the S^(1/2) matrix and
c            get back to the correct partial occupancies
c
c           Diagonalize S.  It's similar to the diagonalization of
c            (H,S), but there is no overlap matrix
c
            call zhpev('V','U',n,smat1,seig,svect,mhe,zwork,rwork,info)
            if(info.ne.0) then
               write(0,*) 'INFO = ',info,' when diagonalizing S'
               ioverr = 1
               return
            end if
ctemp
c$$$            write(*,*) 'Eigenvalues of S = '
c$$$            write(*,*) (seig(i),i=1,n)
cend temp
c
c           We need to extract the s, p, and d "charges" on each
c            atom from each state.  If this is a non-collinear
c            calculation, we also need to do a spin decomposition
c            along the global spin axis.
c
c           The logic goes like this:
c
c           We're solving  H|n> = e(n) S|n>, normalized so that
c            <n|S|m> = delta(m,n) .
c
c           Since S is positive definite, its eigenvalues, o(i),
c            are all positive, and
c
c           S|a> = o(a) |a>, <a|b> = delta(a,b)
c
c           This makes it easy to define a "Square Root Matrix",
c
c           U = Sum[ o(a)^(1/2) |a><a| ,{a,1,n}]
c
c           Note that U is Hermitian
c
c           Such that S = UU.  Note that U^(-1) is also easy to define
c
c           Then our problem becomes
c
c           H |n> = e(n) U U |n> ,
c
c           which is easily written as
c
c           [U^(-1) H U^(-1)] [ U |n> ] = e(n) [ U |n> ]
c
c           Let   Hs = U^(-1) H U^(-1) and |N> = U |n>
c
c           Then
c
c           Hs |N>  = e(n) |N>, and <M|N> = delta(m,n) ,
c
c           so the |N> are our physical eigenstates, and the charge
c            on the ith "orbital" is just |<i|N>|^2
c
c           Note that we don't want charge by orbital, rather we want
c            charge by angular momentum state in each orbital.
c            I hope I've got all of this in the right order.
c
c           There are faster ways of doing this, but let's play safe for
c            now.  Form the U matrix
c
            do i = 1,n
               seig(i) = sqrt(seig(i))
            end do
            do i = 1,n
               do j = 1,i
                  csum = dcmplx(zero,zero)
                  do k = 1,n
                     csum = csum + seig(k)*svect(i,k)*dconjg(svect(j,k))
                  end do
                  u(i,j) = csum
                  u(j,i) = dconjg(u(i,j))
               end do
            end do
         else
c
c           The U matrix is just the identity.  We really shouldn't
c            have to worry about this, but for now, it's relatively
c            fast:
c
            do i = 1,n
               do j=1,n
                  u(j,i) = dcmplx(zero,zero)
               end do
               u(i,i) = dcmplx(one,zero)
            end do
         end if
c
c        Now that U is formed, find the charge for each orbital
c         with eigenvalue < eigcut.  First find out how many
c         states there are
c
         do k = 1,n
            if(en(k).gt.eigcut) then
               ntot = k-1
               go to 40
            end if
         end do
         ntot = n
c
c        Here's where we change the write statements
c
c 40      write(45,43) ak(1),ak(2),ak(3),weight,ntot,weight
c 43      format(4f10.5,2X,i10,d20.10,' K')
 40      write(55) ak(1),ak(2),ak(3),natx,ntot,weight
         do k = 1,ntot
c$$$            write(45,45) en(k),zero
c$$$ 45         format(f10.7,e15.6)
            write(55) en(k)
c
c           Only calculate the actual decomposition if it is
c            meaningful, i.e., when ioverr = 0.  If this isn't the
c            case, we'll set s,p, and d to -1, a flag to the user
c            to indicate that things aren't really working well.
c
c
c           i is the current orbital.  I'm assuming the orbitals are
c            stacked by atom, (s,px,py,pz,d1,d2,d3,d4,d5), with
c            kbas(kkind(j)) basis functions per atom.  If not,
c            your mileage may vary.
c
            i = 0
c
c           Note that in the non-collinear case this is a sum
c            over twice the number of atoms in the system.
c
            do j = 1,natx
c
c              Make sure during each part of the decomposition
c               that there are enough basis functions to
c               support it.  We'll take advantage of the fact
c               that a Fortran 77 or better do loop doesn't
c               execute if its final argument is smaller
c               than its original argument.
c
c              How many basis functions are there for this atom?
c
c              Note that in the non-collinear case there are the
c               same number of spin down as spin up basis functions.
c               Spin down basis functions are below the spin up
c               basis functions.
c
               if(j.gt.natoms) then
                  nbas = kbas(kkind(j-natoms))
               else
                  nbas = kbas(kkind(j))
               end if
c
c              Sum up each orbital for the current atom:
c
               do k1 = 1,nbas
                  i = i + 1
                  if(ioverr.eq.0) then
                     csum = dcmplx(zero,zero)
                     do m = 1,n
                        csum = csum + u(i,m)*zeig(m,k)
                     end do
                     orbocc(k1) = dble(csum*dconjg(csum))
                  else
                     orbocc(k1) = -one
                  end if
               end do
c
c              If mode=5 or 6, sum up orbitals of the same
c               angular momentum
c
               if(mode.ne.10) then
                  s = orbocc(1)
c
                  p = zero
                  do k1 = 2,min(4,nbas)
                     p = p + orbocc(k1)
                  end do
c
                  d = zero
                  do k1 = 5,nbas
                     d = d + orbocc(k1)
                  end do
c
c                 Write to disk, including the number of items
c                  written:
c
                  write(55) 3,s,p,d
c
               else
c
c                 write all of the occupied orbitals to disk:
c
                  write(55) nbas,(orbocc(k1),k1=1,nbas)
               end if
            end do
         end do
c
c     For the other modes, note that "n" is the size of the secular
c      equation, and is properly doubled before diag is called in
c      the non-collinear case.
c
      else if (mode.eq.8) then
c
c        Here we just want to print out the eigenvalues.  Note
c         that I'm using Lattice coordinates, but that probably
c         doesn't matter so long as we're consistent.
c
c$$$         write(45,'(1x,3f9.6,10x,4f9.5/(38x,4f9.5))')
c$$$     $        ak(1),ak(2),ak(3),(en(j),j=1,n)
         write(55) n,ak(1),ak(2),ak(3),(en(j),j=1,n)
      else if (mode.eq.9) then
c
c        Just a standard QLMT file, but no "atomic" information
c
c$$$         write(45,43) ak(1),ak(2),ak(3),weight,n,weight
c$$$         do j = 1,n
c$$$            write(45,45) en(j),zero
c$$$         end do
         write(55) ak(1),ak(2),ak(3),n,weight,(en(j),j=1,n)
      end if
      return
      end
