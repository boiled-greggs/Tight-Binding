      subroutine bande(mode,jspins,flipos,lnoncol,spinat,
     $     ltorque,dspth,dspphi,lspdir,
     $     thisvol,zfill,tkb,eigcut,weight,efermi,earg1,earg2,
     $     smag,lp_par,lp_parv,lp_klst,kptpr,label,iberr)
c
c=======================================================================
c
c     Calculates band structure and sum of eigenvalues, assuming
c      zfill electrons and a Fermi temperature of tkb.  If
c      lqlmt is true, calculates partial occupancies and writes
c      to the QLMT file.
c
c     thisvol is the volume of the unit cell.  It is only needed
c      during pressure calculations (lpress = .true.), when it is
c      passed to setvol.
c
c=======================================================================
c
c     Note:  The QLMT file and perturbations (e.g., the pressure) are
c      calculated only for states with eigenvalues less then eigcut
c
c=======================================================================
c
c     Fixed (hr,hi) -> hmat, (sr,si) -> smat, both complex*16
c
c=======================================================================
c
c     If lpress is true, calculate the pressure using first-order
c      perturbation theory.  Otherwise ignore it.
c
c     lqlmt is true when we want to produce a qlmt file
c     lqlmt0 produces a "qlmt-like" file (one with on 'atomic'
c      information)
c
c     lqapw is true when we want to produce an eigenvalue file
c      in APW format
c
c=======================================================================
c
c     This is a port of the serial code for an f77/MPI based
c      multiprocessor system.  Note that the number of processors is
c      arbitrary, and that the true parallelization takes place only
c      in the bande.f subroutine, where the procedure is to assign a
c      set of k-points to each processor.  Given the setup times
c      involved, this is probably highly scalable until the number
c      of processors gets to about 1/2 the number of k-points.
c
c     When possible, MPI commands will be set apart from the rest of
c      the code by "cMPI/cend MPI" comments.
c
c                                                -- mjm  24 Feb 1998
c=======================================================================
c
c     Version 1.05:
c
c     As part of the general upgrade, jsover is moved into the calling
c      parameters.  Remember that jsover = 1 for non-orthogonal
c      calculations, and 0 for orthogonal calculations.
c
c                                                mjm -- 20 April 1998
c=======================================================================
c
c     Version 1.06:
c
c     jsover has been changed to a logical variable:
c        jsover = .true.  -> Non-orthogonal Hamiltonian (S <> identity)
c        jsover = .false. -> Orthogonal Hamiltonian (S = identity)
c
c     Separated the eigenvector calculation parameters into a parameter
c      file 'P2'.
c
c     Moved k-point generation parameters into 'P3'
c
c                                                mjm --  6 July  1998
c
c     The variable iwrite is no longer with us.  RIP.
c                                                mjm --  4 Aug   1998
c=======================================================================
c
c     Version 1.10:
c
c     Uses Gillan's method (J Phys Cond. Mat. 1, 689-711 (1989),
c      see also Grotheer and Fahnle, PRB 58, 13459-64 (1998)) to
c      "extrapolate" the dependence of the energy on the Fermi
c      broadening temperature to T=0.  Grotheer&Fahnle show that
c      the correction term goes as T**4.
c
c     Pursuant to this, changed the calling sequence to bande,
c      passing "mode" instead of "lpress,lqlmt,lqlmt0,lqapw.  Note
c      that this means we have to change the calling sequence to
c      diag.f as well.
c
c     Also now the Fermi energy, energy, and pressure are passed
c      through the calling statement rather than through the common
c      block /fermcom/.  Note that we call "energy and pressure"
c      "earg1" and "earg2", since in different modes we output
c      different quantities.
c
c     Here are the current modes relevant to bande, with the
c      corresponding output of earg1 and earg2:
c
c     1 -- Prints Gillan's extrapolated E(0) and E(T) in the SKENG file
c          earg1 = extrapolated energy
c          earg2 = E(T)
c     2 -- Prints extrapolated E(0) and P(0) in the SKENG file
c          earg1 = extrapolated energy
c          earg2 = extrapolated pressure
c     3 -- Total Energy, no Pressure
c          earg1 = E(T)
c          earg2 = undefined
c     4 -- Total Energy and Pressure
c          earg1 = E(T)
c          earg2 = P(T)
c     6 -- Total Energy, no Pressure, and produces an LAPW like QLMT file
c          earg1 = E(T)
c          earg2 = undefined
c     8 -- Total Energy, and no Pressure, and produces an APW like
c             eigenvalue file (called QAPW).  (added 9/12/96)
c          earg1 = E(T)
c          earg2 = undefined
c     9 -- QLMT like file with "0" atoms
c          earg1 = E(T)
c          earg2 = undefined
c
c                                                    -- mjm  27 Jan 1999
c-----------------------------------------------------------------------
c
c     Added mode 5, which is identical to mode 6 except that the
c      QLMT file is unformatted.
c                                                    -- mjm  11 Feb 1999
c=======================================================================
c
c     Version 1.11:
c
c     Split the work of the old setup.f into two parts:
c      setpar.f determines the k-point independent quantities
c      setham.f determines the k-point dependency of the Hamiltonian
c
c     Now we can drop the "first, last, ..." variables.
c      The variable lprint is now used to determine the amount of
c      output from setpar and setham.
c
c     This version has to store the onsite, Hamiltonian, and overlap
c      parameters in bande, rather than in setup.
c                                                    -- mjm  30 Apr 1999
c
c     Fixed diag so that it can handle atoms with basis size less than
c      mbas and when the basis size changes between atoms, as in
c      something like NbC or PdH.  This requires kkind and kbas to
c      be passed as arguments to diag.
c                                                    -- mjm  22 Jul 1999
c
c     Split setvol.f into setparv.f and sethamv.f.  These are the
c      volume perturbation analogs to setpar.f and setham.f.
c                                                    -- mjm  23 Jul 1999
c
c     In "Magnetic" mode (jspins=2), calculate the spin up and spin
c      down populations, and so, by extension, the magnetization of
c      the system.  This is printed to SKOUT.        -- mjm  30 Jul 1999
c
c      Print the diagnostic messages from setpar[v]/setham[v] only
c      if we are on the "master" processor.  Note that this can be
c      done without adding MPI calls to those subroutines.
c                                                    -- mjm   4 Aug 1999
c=======================================================================
c
c     Version 1.20:
c
c     If the tight-binding parameters are magnetic (see input.f), then
c      allow for "spin flips."  That is, and the input for the atomic
c      positions is scanned to see if this particular atom is pointing
c      "up" or "down."  By default, the onsite part of spin 1 is
c      associated with up, and spin 2's onsites are associated with spin
c      down.  If an atom is pointed "up", or if "spinflip" is set to
c      .false., then the default situation holds.  If, however, an atom
c      is supposed to be "flipped," then the flipped atom uses the spin 2
c      onsites for "up" and spin 1 for "down."  In this way we can
c      generate an anti-ferromagnetic system.  See input1.f to see how
c      the input changes.
c
c     To help in this effort, the variable "jspins" has been moved from
c      to common block /lia/ and moved into the calling sequence of
c      input.f, input1.f, and bande.f.
c
c     Also note the new logical array, "flipos".  "flipos(i)" is false
c      in the default case, i.e., the spin of atom i points in the spin
c      direction we agree to be called "up."  If flipos(i) is false,
c      atom i is "pointing down," and so we flip its onsite spins.  See
c      "input1.f" and "setpar.f" to see how this happens.
c                                                    -- mjm  21 Dec 1999
c-----------------------------------------------------------------------
c
c     bande.f now reports back the magnetic moment, if any, in the
c      real*8 variable "smag".
c                                                     -- mjm  8 Feb 2000
c=======================================================================
c
c     Version 1.21:
c
c     The diagnostic keys lp_par, lp_parv, and kptpr are now read
c      in from the TWEAKS file:
c                                                     -- mjm  7 Aug 2000
c=======================================================================
c
c     Version 1.22:
c
c     While working on version 1.30, I realized that I was going to
c      have to schlep a bunch of eigenvalues and various eigenvalue
c      derivatives from the current processor back to the master
c      processor.  I then realized that this wasn't necessary.  Thus we
c      come back to this version, which will be folded into version
c      1.31.
c
c     Thus
c
c     When calculating perturbations on the eigenvalues, either for
c      pressure or when calculating torques in 1.3+ (or forces, in 1.?),
c      derivatives of the eigenvalues have to be stored.  In older
c      versions of the code these, and the eigenvalues, were then
c      transferred to the master processor in order to find the Fermi
c      level, the total energy, and whatever derivatives we wanted.
c      While these calculations were taking place the other processors
c      were idle.  While this isn't bad for small jobs, it is a major
c      drawback for large systems.  Thus here we rewrite bande.f and
c      fermi.f to keep eigenvalues and derivatives on the originating
c      processors and to only transfer the minimum data necessary back
c      to the original processor.
c
c                                                  -- mjm 19 October 2000
c=======================================================================
c
c     Version 1.30:
c
c     Adding non-collinear magnetization, following Pickett (1996).
c
c     Specifically,
c
c     If lnoncol is true and jspins = 2, then non-collinear
c      magnetization is enabled.  (lnoncol = true and jspins = 1
c      aborts the program.)  In this case, the array spinat(j,iatom)
c      holds the Cartesian coordinates of the spin direction unit
c      vector of atom i.  Note that bande does not check to see that
c      spinat is actually normalized.
c
c     Note that lnoncol overrides the flipos array.
c
c     Changes required to implement this.  Unless otherwise specified,
c      all of the following assumes that lnoncol is true:
c
c     (a) The array paros will now hold the <<average>> of majority
c         and minority onsite terms.
c
c     (b) The new array spinos = paros(majority) - paros(minority)
c         will hold the spin information.  This is the part that
c         will be rotated to tilt the spins.
c
c     (c) use spinos and spinat (the direction of each atomic spin)
c         to compute the diagonal correction to the spin up
c         on-site terms, spndiag, and the off-diagonal correction
c         spndiag.  These are calculated in the subroutine spinham,
c         added to the Hamiltonian (at each k-point) in addspins.
c
c     (d) The number of electrons/band, sfac, must be adjusted for
c         the various combinations of jspins and lnoncol.  sfac is
c         now given as an input to fermi.f.
c
c                                                  -- mjm 10-11 Aug 2000
c
c     (e) Note that lpress = .true. will calculate the pressure of
c         the system ASSUMING FIXED SPIN ORIENTATION.
c                                                     -- mjm 14 Aug 2000
c
c     (f) If ltorque is true, use perturbation theory to determine the
c         d E/d theta(i) and d E/d phi(i) for each atom, storing the
c         results in the file SKDERIV (unit 41).  Note that you will
c         get the wrong answers unless you use a symmetry file which
c         assumes that no two atoms are identical.  To make this work
c         we input two new arrays, dspth = d spinat/d theta, and
c         dspphi = d spinat/d phi, both in units of energy/degree.
c                                                  -- mjm 14-15 Aug 2000
c-----------------------------------------------------------------------
c
c     iberr returns a non-zero value if something has gone wrong in the
c      bande calculation.
c     iesp checks the error conditions assigned by setpar.
c                                                     -- mjm 10 Aug 2000
c-----------------------------------------------------------------------
c
c     Changed the formatting of the QLMT file when the eigenvalue is
c      <= -10 or >= 100 so that we don't get overflows
c                                                     -- mjm  3 Jul 2003
c=======================================================================
c
c     Version 1.31:
c
c     Merger of the 1.21->1.22 and 1.21->1.30 changes.
c                                                     -- mjm 23 Oct 2000
c
c     If lspdir is true:
c
c     We specify the direction of the spin of each atom in the input.
c      But this does not guarantee that we know the direction of each
c      spin in the output.  Here we turn on an option to determine the
c      "output spin direction" of each atom.  Note that at
c      self-consistency the input and output spin directions should be
c      aligned, so this is one way of determining the directions of
c      spins at a domain wall, or in a frustrated anti-ferromagnet.
c      See input1.f for directions on how to turn on this calculation.
c
c      Results are in the file SPNDIR (unit 42).  Note that this is
c       open only on the master node.
c                                                -- mjm  26 October 2000
c
c     Updated "tweaks" to provide more information on how the program
c      is progressing through the k-point list.  This adds kp_klst
c      to this routine, which forces printing of the current k-point
c      and processor to standard output.
c                                                  -- mjm 28 August 2001
c=======================================================================
c
c     Version 1.32:
c
c     Mode 10 prints out the occupation of every orbital in the system.
c       Note that what is written to (55) now changes a bit.
c
c                                                    -- mjm   1 Nov 2001
c=======================================================================
c
c     Version 1.33:
c
c     Change bande and its associated files to only allocate eigenvalue
c      arrays for the number of k-points which will be used on one
c      processor, rather than dimensioning for all k-points on each
c      processor.  Note that this will require recompiling the code
c      whenever we DECREASE the number of processors in the calculation.
c      The relevant parameter is "minproc", which is dimensioned in P3.
c
c     This involves replacing some of the nkd dimensions by nkdproc,
c      also defined in P3, and adding auxillary local pointers to local
c      do-loop variables.
c
c                                                    -- mjm  15 Aug 2002
c
c     When lpress is true, if setparv fails it puts its error code
c      into the variable iespv.  bande detects this, sets its own
c      error code, and returns.  The check for non-zero iepsv was done
c      whether or not lpress was true or not.  This can cause problems
c      for compilers which actually leave uninitialized variables
c      uninitialized.  We now only check iespv if setparv was actually
c      called.
c                                                    -- mjm  21 Aug 2002
c=======================================================================
c
      implicit real*8 (a-h,o-z)
cMPI
      include 'mpif.h'
      logical master,multiprc
      common /MPIcomm/ MPISIZE,MPIRANK,master,multiprc
cend MPI
c
      include 'P1'
      include 'P2'
      include 'P3'
c
c     Magnetization options:
c
c     jspins = 1 for unpolarized, 2 for polarized
c
      integer jspins
c
c     flipos is described in the Version 1.20 notes, above
c
      logical flipos(matom)
c
c     lnoncol, ltorque, spinat, dspth and dspphi are described the
c      in Version 1.30 notes, above.
c
c     lspdir is described in the Version 1.31 notes
c
      logical lnoncol,ltorque,lspdir
      real*8 spinat(3,matom),dspth(3,matom),dspphi(3,matom)
c
c-----------------------------------------------------------------------
c
c     iberr is returned non-zero if there is some kind of error in
c      the calculation.  iesp is a similar variable returned by setpar
c
      integer iberr,iesp
c
c     It is still useful to know if the pressure is to be calculated:
c
      logical lpress
c
c     It is also useful to know when we need to calculate the entropy:
c
      logical lentro
c
      real*8 ak(3)
      real*8 weight(nkd)
c
c     npts is the number of k-points
c     jsover = .true.  when there is an overlap matrix,
c            = .false. for orthogonal tight-binding.
c
      logical jsover
c
      common /lia/ npts,jsover
c
c     These are the k-points, as stored before this subroutine is called
c
      common /cd1/ x(nkd),y(nkd),z(nkd)
c
c     These are, of course, the tight-binding parameters:
c
      common /params/ param(npard,mspin)
c
c     nv is the size of the secular equation for one spin:
c
      common /struc/ nv
c
c     At the moment, these common blocks are only required here
c      so that we can call kkind and kbas in diag.f.
c      Note that setham.f gets its information from the common
c      blocks themselves.
c
c     ksk(i) tells setham where the integrals for atom i start
c      in the Hamiltonian and overlap matrices
c
      common /kstuff/ ksk(matom)
      common /codes/posn(matom,3),kkind(matom),natoms
      common /parinfo/ valence(mkind),kinds,kbas(mkind)
c
c$$$      common /neigh/ dscrn(mpair)
c
c     Unperturbed Hamiltonian and overlap matrices, in packed triangular
c      form.
c
      complex*16 hmat(mh*(mh+1)/2),smat(mh*(mh+1)/2)
c
c     spndiag is the diagonal correction to the Hamiltonian for
c      non-collinear spins.  spnmix is the off-diagonal correction.
c      see spinham.f and addspins.f for more details.  Note that
c      for non-collinear calculations mh will be double the size
c      of the non-polarized secular equation, so mh/mcol is meaningful.
c
      real*8 spndiag(mh/mcol)
      complex*16 spnmix(mh/mcol)
c
c     spndth(i) is the derivative of spndiag(i) with respect to
c      theta(i)
c     spnmxth(i) is the derivative of spnmix(i) wrt theta(i)
c     spndph(i) is the derivative of spnmix(i) wrt phi(i)
c     All derivatives are calculated in degrees.
c     The derivatives are calculated by spinham, but only when
c      lnoncol AND ltorque are true.
c
      real*8 spndth(mh/mcol)
      complex*16 spnmxth(mh/mcol),spndph(mh/mcol)
c
c     dspndg and dspnmx are the versions of spndiag and spnmix used
c      to find the pressure dependence of each eigenvalue
c
      real*8 dspndg(mh/mcol)
      complex*16 dspnmx(mh/mcol)
c
c     Eigenvalues, stored by k-point and spin.  Note that (in
c      Versions 1.22 and above) the eigenvalues never leave their
c      home processor.  Thus really nkd could be replaced by
c      nkd/#nodes, if we knew the number of nodes at compile time.
c      Or we could dynamically dimension this, if we knew how to
c      do it in Fortran 77.
c
c     Version 1.33:  Well, we don't actually know the number of nodes,
c      but we're going to assume that there are at least minproc of
c      them, and reduce the dimension from nkd to
c
c     nkdproc = nkd/minproc + 1
c
      real*8 esk(mh,nkdproc,mspin)
c
c     temporary eigenvalue array needed when we output to a QLMT or
c      QAPW file
c
      real*8 en(mh)
c
c     den(i,j,ispin) is the "density" on atom i from atoms of type j.
c      This is defined in setpar.f.  It is stored here only because
c      the routine setparv.f needs this density to perform the
c      perturbation setup.  Version 1.30 change:  we need both
c      a majority (ispin = 1) and minority (ispin = 2) density
c      in order to find paros and spinat.  For other calculations
c      ispin = 1, as we don't need to keep as we loop over spins.
c
      real*8 den(mkind,matom,mspin)
c
c     We need this stuff for the perturbation calculations.  Note
c      that we can make the dimensions smaller by setting mpress = 1
c      in P1:
c
c     Eigenvectors (unperturbed) and workspace
c
      complex*16 zeig(mhe,mhe),work(mhp)
c
c     Perturbed eigenvalues (one set for each perturbation)
c     Version 1.33: Note redimensioning of nkd
c
      real*8 desk(mhp,nkdproc,mspin)
c
c     How many perturbed states are we looking at?
c     Version 1.33: Note redimensioning of nkd
c
      integer ncut(nkdproc,mspin)
c
c     We're going to use MPI_GATHERV to collect all of the eigenvalues
c      back together.  To do this, we'll need to know the offset of
c      each eigenvector, hence the array offset, which is dimensioned
c      for a large number of processors.  We'll also need the
c      count in each processor, though one implies the other, doesn't
c      it?
c     (After version 1.22, we use MPI_REDUCE, and never collect all of
c      the eigenvalues back on one processor.  Does this mean that
c      we don't need offset, and count?  Think about deleting this
c      in future versions of the code.)
c
c     Version 1.33: yes, I think we can get rid of all this.  But
c      keep as comments for now, just to be safe:
c$$$c
c$$$      integer offset(256),count(256)
c$$$c
c$$$c     The pressure calculation needs similar offsets, both for the
c$$$c      eigenvalue perturbations and for ncut:
c$$$c
c$$$      integer poffst(256),pcount(256),ncoff(256),ncount(256)
c
c     Note that we do have to have some analog of count(i) for printing
c      QLMT files.  See the end of this subroutine.
c
c     character string for temporary file:
c
      character*16 qstring
c
c     Logical variables controlling printing from setpar/setham etc.
c
      logical lp_par,lp_parv,lp_klst
c-----------------------------------------------------------------------
c
c     Version 1.31 and up:
c
c     When lspdir = .true.
c
c     On output, each atom will have a spin vector associated with the
c      expectation value of the spin operator in the outgoing eigenstate.
c     We'll call this spin vector spdir(i,j), where i=1,2,3=x,y,z, and
c      j is that atom number.  Note that we will need both local
c      and global versions of spdir.  We'll add the local ones up with
c      an MPI_REDUCE operation.
c
c     We also need an array for the spins coming from the eigen
c
c     We could do an s,p,d or s,p,t2g,eg decomposition of the spin,
c      but we'll leave that out for now.
c
c     The spdir array will be printed to unit 42, which static.f has
c      opened as file SPNDIR.
c
c     Note that spdirl and lspdir have different meanings!
c
      real*8 spdir(3,matom), spdirl(3,matom)
c
c     We will also need to temporarily store the decomposition from
c      (yuck) every eigenvalue at every k-point.  Note that the
c      dimensions of this array are controlled by the the value of
c      msdir, set in P2.  When we don't want the decomposition we
c      can set msdir=0, and these dimensions will save a significant
c      amount of space.
c
c     The first argument is the direction, then the atom, then the
c      eigenvector, and finally the k-point.  Note that the dimensions
c      are written so as to minimize the space taken up by spev if
c      we are not looking at output directions.  (In P2, msdir = 1
c      or 0.)
c     Version 1.33: Note redimensioning of nkd
c
      real*8 spev(3,msdir*(matom-1)+1,msdir*(mh-1)+1,
     $     msdir*(nkdproc-1)+1)
c
c     orbocc(i) is the occupation of the ith orbital, or group of
c      orbitals on the current atom, as determined in diag and written
c      to "tape" 55.
c
      real*8 orbocc(mbas)
c
c     spsize is the total size of the spdir and spdirl arrays.  It will
c      be needed by MPI_REDUCE
c
      integer spsize
      parameter (spsize = 3*matom)
c
c     We import the label, in case we want to print things here:
c
      character*20 label
c
c-----------------------------------------------------------------------
c
c     Version 1.32:
c
c     kforb holds a format statement for QLMT output
c
      character*9 kforb
c
c-----------------------------------------------------------------------
c
c     Parameters required by the setpar/setham combination:
c
c     Parameters required by setparv/sethamv are distinguished by
c      a d-prefix
c
c
c     paros(1,iat) is the onsite "s" orbital for the iat-th atom.
c     2 = p, 3 = t2g, 4 = eg.
c
      real*8 paros(4,matom),dparos(4,matom)
c
c     spinos(i,iat) is the difference between the spin up and
c      spin down onsite terms when lnoncol is true.
c      dspnos is the obvious generalization:
c
      real*8 spinos(4,matom),dspnos(4,matom)
c
c     hop(i,j) is the ith Slater-Koster Hamiltonian parameter
c      for the j-th pair in the list produced by search2.f.
c
      real*8 hop(mpkind,mpair),dhop(mpkind,mpair)
c
c     over is the corresponding quantity for the overlap matrix.
c      note that if jsover is false these are not calculated.  We
c      may split these off into a separate subroutine at a later point.
c
      real*8 over(mpkind,mpair),dover(mpkind,mpair)
c
c=======================================================================
c
      parameter (zero = 0d0)
      parameter (half = 5d-1)
      parameter (one  = 1d0)
      parameter (two  = 2d0)
      parameter (ten  = 1d1)
      parameter (twenty = 2d1)
      parameter (hundred = 1d2)
c
      parameter (small = 1d-10)
c
c     This is the same minimum temperature as is used in fermisk.f.
c      Logically it should be read into fermisk on the input line,
c      but at the moment it is not.
c
      real*8 nearz,nearzi
      parameter (nearz = 1d-5)
      parameter (nearzi = one/nearz)
c
c     Convert from radians to degrees:
c
      parameter (pi = 3.141592653589793116d0)
      parameter (rad2deg = 180d0/pi)
c
c=======================================================================
c
c     Fermi weight function and its derivative.  Note that
c      tkbi = 1/tkb
c
      fb(xx) = one/(one+exp(tkbi*xx))
      dfb(xx) = -tkbi*fb(xx)*fb(-xx)
c
c=======================================================================
c
c     Set desired printing options.  ".true." prints diagnostics,
c      ".false." doesn't.  Note:  for versions > 1.10, this is read
c      from TWEAKS file, if any.
c
c$$$      data lp_par /.true./
c$$$      data lp_par /.false./
c$$$      data lp_parv /.false./
c
c=======================================================================
c
ctemp
c$$$      write(*,*) 'You are in bande on processor ',MPIRANK
cend temp
c
c     No errors so far:
c
      iberr = 0
c
c     Check jspins/lnoncol compatibility:
c
      if(lnoncol.and.(jspins.eq.1)) then
         write(0,*)
     $   'bande:  Calling for non-collinear magnetization calculation'
         write(0,*) 'with non-polarized parameters.'
         write(0,*) 'check your parameters.'
         iberr = 1
         return
      end if
c
c     In the noncollinear case, the size of the secular equation will
c      double, but we will not loop through the spins.  Thus define:
c
c     nvx = the size of the secular equation for the current problem,
c      noncollinear spins or not.
c
c     jsx = the terminus of loops over jspin.
c
      if(lnoncol) then
         nvx = 2*nv
         jsx = 1
      else
         nvx = nv
c
c        jspins = 2 for collinear spin-polarized, 1 for non-polarized
c
         jsx = jspins
      end if
c
c     Should we calculate a pressure?
c
      lpress = (mode.eq.2).or.(mode.eq.4)
c
c     Should we calculate the entropy?
c
      lentro = mode.lt.3
c
cMPI
c     Divide the calculation amoung processors as easily as possible
c
c     kproc is the minimum number of k-points per processor
c
      kproc = npts/MPISIZE
c
c     krem is the remainder:
c
      krem = npts - kproc*MPISIZE
c
c     Starting and stopping points for the current processor
c      To keep the load balanced, the first krem processors will
c      do kproc+1 k-points, and the remainder will do only
c      kproc k-points.  Remember that MPIRANK = 0 is a processor,
c      and the last one is MPISIZE-1
c
      if(MPIRANK.lt.krem) then
         istart = MPIRANK*(kproc+1)+1
         ifin = istart+kproc
      else
         istart = MPIRANK*kproc+krem+1
         ifin = istart+kproc-1
      end if
ctemp
c$$$      write(*,*) MPIRANK,npts,kproc,istart,ifin,krem
cend
c
c     Version 1.33:  this can be deleted, as far as I can tell:
c
c$$$c     Count the offsets.  Note that at each k-point we have
c$$$c      mh eigenvalues to ship out (mhp for the pressure),
c$$$c      so the count at each processor is mh(mhp)*number of processors
c$$$c
c$$$c     Note that when we are finished we should have
c$$$c
c$$$c     count(MPIRANK+1) == mh*(ifin-istart+1)
c$$$c
c$$$c     for each process.  Note that we need none of this if
c$$$c      we are only on a single processor:
c$$$c
c$$$      if(multiprc) then
c$$$         offset(1) = 0
c$$$         poffst(1) = 0
c$$$         ncoff(1)  = 0
c$$$c
c$$$c        Remember the remainder problem.
c$$$c
c$$$         if(krem.gt.0) then
c$$$            count(1)  = mh  * (kproc+1)
c$$$            pcount(1) = mhp * (kproc+1)
c$$$            ncount(1) = kproc+1
c$$$         else
c$$$            count(1)  = mh  * kproc
c$$$            pcount(1) = mhp * kproc
c$$$            ncount(1) = kproc
c$$$         end if
c$$$         do i = 2,MPISIZE
c$$$            if(i.gt.krem) then
c$$$               count(i)  = mh  * kproc
c$$$               pcount(i) = mhp * kproc
c$$$               ncount(i) = kproc
c$$$            else
c$$$               count(i)  = mh  * (kproc+1)
c$$$               pcount(i) = mhp * (kproc+1)
c$$$               ncount(i) = kproc + 1
c$$$            end if
c$$$            offset(i) = offset(i-1) +  count(i-1)
c$$$            poffst(i) = poffst(i-1) + pcount(i-1)
c$$$            ncoff(i)  = ncoff(i-1)  + ncount(i-1)
c$$$         end do
c$$$ctemp
c$$$c$$$         write(*,*) 'Proc. ',MPIRANK,' count = ',count(MPIRANK+1),
c$$$c$$$     $        ' should be ',mh*(ifin-istart+1)
c$$$cend temp
c$$$      else
c$$$c
c$$$c        Even on a serial machine we need count(1) if we want
c$$$c         to use the current algorithm for printing out eigenvalues:
c$$$c
c$$$         count(1) = mh * kproc
c$$$c
c$$$      end if
c
      do jspin=1,jsx
c
c        if one of the "lq" options is true, we're going to open up a
c         binary file to write the eigenvalue information for the
c         current processor.  At the end of all of this we'll use
c         processor zero to concatenate the files into the traditional
c         QLMT or QAPW file.  Note that the header of this file (unit
c         45) has already been written in static.f.
c
c        The calling program should have eliminated mode=7, so this
c         works:
c
         if(mode.gt.4) then
c
c           We need a processor descriptive name.  Note that more than
c            10^6 processors indicates that you have more money than you
c            can possibly use.
c
            write(qstring,'(''qeigtemp.'',i1,i6.6)') jspin,MPIRANK
c
c           Write over any existing file with this name:
c
            open(unit=55,file=qstring,status='unknown',
     $           form='unformatted')
         end if
c
         if(master) write(15,*)' SPIN=',jspin
c
c        Set up the Slater-Koster parameters.  These are
c         independent of k-point, and, in MPI mode, can be calculated
c         independently on each processor.  Eventually we probably
c         want to move this part of the calculation out of bande
c         and back into static, right below search2.
c        It may be efficient to parallelize this part of the code
c        For now print diagnositics as the default.
c
ctemp
c$$$         write(*,*) 'Calling setpar on processor ',MPIRANK
cend temp
c
c        Note that (version 1.20 and above) setpar and setparv
c         need the entire parameter array, not just one spin
c
         call setpar(jspin,param,jsover,den,flipos,lnoncol,
     $        paros,spinos,hop,over,master.and.lp_par,iesp)
c
c        Quit on error:
c
         if(iesp.ne.0) then
            iberr = 2
            return
         end if
c
c        Similar moves for the pressure terms, if needed.  Note
c         that in Harrison sign mode setparv needs the original hopping
c         and overlap parameters
c
         if (lpress) then
            call setparv(jspin,param,thisvol,jsover,den,
     $           flipos,lnoncol,dparos,dspnos,dhop,dover,
     $           master.and.lp_parv,iespv)
c
c           Note:  before 1.33, this was outside the lpress if block,
c            which could lead to problems if the iespv starts
c            life undefined:
c
c           Quit on error:
c
            if(iespv.ne.0) then
               write(0,*) 'bande:  setparv reported error ',iespv
               iberr = 3
               return
            end if
         end if
c
c        We now have all of the information needed to construct the
c         spin-polarization part of the Hamiltonian and overlap
c         matrices, if lnoncol is true.  Note that if ltorque is
c         true the program will also return derivative information,
c         in the vectors spndth, spnmxth, and spndph
c
         if(lnoncol) then
            call spinham(natoms,kkind,kbas,spinat,spinos,
     $           spndiag,spnmix,ltorque,dspth,dspphi,
     $           spndth,spnmxth,spndph)
c
c           We can use the same routine to calculate the pressure
c            derivative terms:
c
            if(lpress) call spinham(natoms,kkind,kbas,spinat,dspnos,
     $           dspndg,dspnmx,.false.,dspth,dspphi,
     $           spndth,spnmxth,spndph)
         end if
c
c
c        Now start the do loop.  Note that if MPISIZE = 1 this will
c         reduce to the serial case.
c
c        Version 1.33: note that we must keep track of the local
c         k-point, iloc, as well as the global k-point i.  In
c         a future world we'll fix this up:
c
         iloc = 0
c
         do i = istart,ifin
c
            iloc = iloc + 1
c
c           Printed if TWEAKS file has PRINT_PROGRESS keyword set:
c
            if(lp_klst) write(*,*)
     $           'Processor ',MPIRANK,' is working on k-point',i,
     $           ' spin ',jspin
c
c           The k-point coordinates are stored globally:
c
            ak(1)=x(i)
            ak(2)=y(i)
            ak(3)=z(i)
c
c           Set up the Hamiltonian and Overlap matrices.
c            Output diagnostics only when the final argument is .true.
c             note that means that you can designate a particular
c             k-point or range of k-points for diagnostic treatment.  At
c             the moment all we have is i<=kptpr, where kptpr is set in
c             tweaks.f.  This could be the place for a creative use of
c             tweaks, if you want to, say, print out the eigenvalues
c             from the 10th, 17th, and 32nd k-points.
c
c           Version 1.30 and up:  When lnoncol = .true., setham only
c            calculates the averaged single-spin Hamiltonian.
c
c           Note that the Hamiltonian is printed out in the appropriate
c            place
c
c           Version 1.33:
c            Note that the "i" here is the global k-point count:
c
            call setham(jsover,ak,paros,hop,over,hmat,smat,
     $           (.not.lnoncol).and.master.and.(i.le.kptpr))
c
c           Now add in the spin corrections for non-collinear spins.
c            Note that this will double the occupied size of hmat and
c            smat.  Also note that we really mean nv, not nvx.
c
            if(lnoncol) call addspins(nv,spndiag,spnmix,hmat,smat,
     $           master.and.(i.le.kptpr))
c
c
c           Diagonalize
c
c           Note that there is a lot of wasted space in esk if
c            we are using more than one processor.  For now,
c            however, we'll just put the variables in where they
c            would go if this was all done on one processor.
c
c           Version 1.33:  we fix this problem by working in
c            the local processor space, except for the weight,
c            which is stored globally (and is only needed to
c            output to QLMT/QAPW files).
c
            call diag(mode,lnoncol,hmat,smat,esk(1,iloc,jspin),
     $           lspdir,spev(1,1,1,iloc),
     $           nvx,zeig,jsover,ak,eigcut,weight(i),ioverr,
     $           natoms,kkind,kbas,ksk)
c
c           if ioverr is > 0, the overlap matrix for this k-point
c            is not physical (negative eigenvalues)
c           if ioverr is > 0, the overlap matrix for this k-point
c            is not physical (negative eigenvalues)
c
            if(ioverr.lt.0) then
               write(*,*)
     $              'bande:  logical inconsistency in calling diag'
               iberr = 4
c
c              No point in continuing
c
               return
            end if
c
            if(ioverr.ne.0) then
               write( *,*) 'Found non-physical overlap at k-point',i
c
c              No point in continuing with this k-point.
c
               go to 1000
            end if
c              
ctemp
c$$$            write(*,*) MPIRANK,' k = ',i,esk(1,i,jspin)
cend temp
ctemp
c$$$            if(i.le.2) then
c$$$               write(15,4113) ak(1),ak(2),ak(3),nvx
c$$$ 4113          format(3f10.6,i10)
c$$$               write(15,4114) (esk(j,i,jspin),j=1,nvx)
c$$$ 4114          format(8f10.6)
c$$$            end if
cend temp
c
            if(lpress) then
c
c              Perturbation setup (Volume Perturbation)
c               Note that all we need from the original calculation
c               are the eig and zeig.  Thus we can use hmat and smat
c               to hold the perturbation matricies.  Note that
c               we now pass the volume along here:
c
               call sethamv(jsover,ak,dparos,dhop,dover,
     $              hmat,smat,master.and.(i.eq.0))
c
c              Now add in the spin corrections for non-collinear spins.
c               Again note that we can use the addspins routine
c               for these changes, and that hmat and smat come
c               out twice as big as they started.
c
               if(lnoncol) call addspins(nv,dspndg,dspnmx,hmat,smat,
     $              master.and.(i.le.kptpr))
c
c              Calculate the perturbation of each eigenvector.
c              Again note the wasted space.  But at least we aren't
c               moving the eigenvalues from processor to processor.
c
c              Version 1.33: again reduce space:
c
               call perturb (mhp,nvx,hmat,smat,
     $              zeig,esk(1,iloc,jspin),eigcut,ncut(iloc,jspin),
     $              desk(1,iloc,jspin),work)
            end if
 1000       continue
         end do
cMPI
c        Close the current eigenvalue collector, if it is open.
c
         if(mode.gt.4) close(55)
c
c        End the jspin loop:
c
      end do
c
c     Now use the rewritten fermi.f routine to calculate the
c      Fermi level.
c
c     Use T broading, find E_Fermi using the new Fermi function.  Note
c      that the new Fermi routine is written in parallel.  There will be
c      some computaional inefficiency here, but at least we didn't have
c      to move all of the eigenvalues and derivatives around.
c
c     sfac tells us how many electrons fill each band.  In the
c      paramagnetic case (jspins=1,lnoncol=.false.) there are
c      two electrons in each band.  Otherwise, there is only one.
c
      if((jspins.eq.2).or.lnoncol) then
         sfac = one
      else
         sfac = two
      end if
ctemp
c$$$      write(*,*) 'Finding the Fermi energy on ',MPIRANK
cend temp
c
c     Version 1.33: Note we must pass both nkd and nkdproc
c
      efermi = fermi(mh,nkd,nkdproc,istart,ifin,nvx,jsx,sfac,
     $     weight,esk,tkb,zfill)
ctemp
c$$$      write(*,*) 'E_Fermi = ',efermi,' on ',MPIRANK
cend temp
c
c     We need a little broadening so that the program won't crash.
c      Use nearz as the effective limit of tkb:
c
      if(tkb.lt.nearz) then
         tkbi = nearzi
      else
         tkbi = one/tkb
      end if
c
c     Find the energy.
c     In addition, calculate the change in the sum and the Fermi
c      level due to a volume change for pressure calculation
c
c-----------------------------------------------------------------------
c     New as of version 1.10:
c
c     M. J. Gillan (J. Phys. Condensed Matter 1, 689-711 (1989))
c      showed that the conventional "Fermi broadening" technique can
c      be expressed as a minimization of the free energy:
c
c     A(T) = E(T) - T S(T)
c
c     Where E(T) is the Fermi broadened energy:
c
c     E(T) = Sum_i [w_i f((e_i-mu)/T) w_i]
c
c     and S(T) is the entropy:
c
c     S(T) = - Sum_i [w_i { f((e_i-mu)/T) ln f((e_i-mu)/T) +
c            ( 1 - f((e_i-mu)/T) ) ln ( 1 - f((e_i-mu)/T) )}]
c
c     here:
c
c     the sums are over all eigenstates of the Hamiltonian.  The
c     eigenvalue of the ith state is e_i.
c
c     w_i is the weight given to the k-point at which the ith state
c     is calculated.  For spin-restricted calculations w_i also includes
c     the factor of 2 to represent both electrons.
c
c     f(z) = [1 + exp(z)]^-1 is the Fermi weight function.  Note that
c
c     f(-z) = 1 - f(z)
c
c     f'(z) = -f(z) f(-z)
c
c     f(-z) = exp(z) f(z)
c
c     mu is the chemical potential, defined so that
c
c     Sum_i [w_i f((e_i-mu)/T)] = zfill, the number of electrons
c      in the system.
c
c     Note that e_i depends on the structure of the solid, and, in
c      particular, upon the volume.
c
c     Also, mu depends upon the temperature, T, directly, and upon
c      volume through its relationship with e_i, since zfill is
c      independent of the volume and the temperature.
c
c     Gillan then stated without proof that
c
c     E(T) = E(0) + 1/2 gamma T**2 + O[T**3]
c
c     and
c
c     A(T) = E(0) - 1/2 gamma T**2 + O[T**3]
c
c     it follows that
c
c     E(0) = 1/2 [A(T) + E(T)] + O[T**3]
c          = E(T) - 1/2 T S(T) + O[T**3]
c
c     O. Grotheer and M Fahnle (Phys. Rev. B 58, 13459-64) showed that
c      in fact we can replace O[T**3] by O[T**4].  This makes the method
c      rather powerful, since it means that we can extrapolate from
c      a relatively high temperature to T = 0.  This <em>might</em> save
c      on k-points, since high temperature results converge faster
c      than low temperature results, but this has to be checked.
c
c-----------------------------------------------------------------------
c
c     Other parameters independent of jspins:
c
      etol = twenty*tkb
      efp = efermi+etol
      efm = efermi-etol
c
      if(lpress) then
c
c        Note that these sums are performed on each processor, after
c         which we will us MPI_REDUCE to figure out the global sums
c
         esum = zero
         desum3 = zero
         dfbsum = zero
         desum1 = zero
         desum2 = zero
         dferm = zero
c
c        Note that states far below and far above efermi do
c         not contribute to the entropy
c
         if (lentro) then
            entropy = zero
            dent1 = zero
            dent2 = zero
         end if
c
         do jspin = 1,jsx
c
c           Version 1.33: keep track of the local k-point
c
            kloc = 0
c
            do k = istart,ifin
c
               kloc = kloc + 1
c
c
c              Weights are constant across a single k-point
c              Version 1.33: weights stored globally
c
               wt = sfac*weight(k)
c
c              Remember that only a certain subset of eigenvalues
c               were used.  We hope that these include all energies
c               less than efp
c
c              Version 1.33: ncut,esk,desk stored locally:
c
               do i = 1,ncut(kloc,jspin)
                  if(esk(i,kloc,jspin).lt.efp) then
                     ewt = wt*esk(i,kloc,jspin)
                     if(esk(i,kloc,jspin).lt.efm) then
                        esum = esum + ewt
                        desum1 = desum1 + wt*desk(i,kloc,jspin)
                     else
                        fbe = fb(esk(i,kloc,jspin)-efermi)
                        esum = esum + ewt*fbe
                        dfbe = dfb(esk(i,kloc,jspin)-efermi)
                        dfbsum = dfbsum + dfbe*wt
                        desum3 = desum3 + dfbe*ewt
                        dewt = wt*desk(i,kloc,jspin)
                        desum1 = desum1 + dewt*fbe
                        dferm = dferm + dfbe*dewt
                        desum2 = desum2 + dfbe*ewt*desk(i,kloc,jspin)
                        if (lentro) then
                           entropy = entropy - wt*
     $                          (fbe*log(fbe)+(one-fbe)*log(one-fbe))
c
c                          Note that dS/df = log(f/(1-f)) = -(e-ferm)/T
c
                           con1 = dfbe*wt*(esk(i,kloc,jspin)-efermi)
                           dent1 = dent1 + con1*desk(i,kloc,jspin)
                           dent2 = dent2 + con1
                        end if
                     end if
                  end if
               end do
            end do
         end do
c
c        Transfer everything needed back to the master processor.
c        Note that there is no reason to do everything on
c         every processor.
c        Only the pressure-related reductions are done here.
c        The rest are done outside of this if then/else block.
c
         if(multiprc) then
            call MPI_REDUCE (desum1,desum1g,1,MPI_DOUBLE_PRECISION,
     $           MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE (desum2,desum2g,1,MPI_DOUBLE_PRECISION,
     $           MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE (desum3,desum3g,1,MPI_DOUBLE_PRECISION,
     $           MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE (dfbsum,dfbsumg,1,MPI_DOUBLE_PRECISION,
     $           MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE (dferm,dfermg,1,MPI_DOUBLE_PRECISION,
     $           MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE (dent1,dent1g,1,MPI_DOUBLE_PRECISION,
     $           MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE (dent2,dent2g,1,MPI_DOUBLE_PRECISION,
     $           MPI_SUM,0,MPI_COMM_WORLD,ierr)
            if(master) then
c
c              Put things back where they belong:
c
               desum1 = desum1g
               desum2 = desum2g
               desum3 = desum3g
               dfbsum = dfbsumg
               dferm  = dfermg
               dent1  = dent1g
               dent2  = dent2g
            end if
         end if
      else
c
c        It's easier if we don't need the pressure:
c
         esum = zero
c
c        Note that states far below and far above efermi do
c         not contribute to the entropy
c
         if (lentro) entropy = zero
c
         do jspin = 1,jsx
c
c           Version 1.33:  again, keep local k-point for the
c            eigenvalues, global for the weights:
c
            kloc = 0
c
            do k = istart,ifin
c
               kloc = kloc + 1
c
               do i = 1,nvx
                  if(esk(i,kloc,jspin).lt.efp) then
                     wt = sfac*weight(k)
                     ewt = wt*esk(i,kloc,jspin)
                     if(esk(i,kloc,jspin).lt.efm) then
                        esum = esum + ewt
                     else
                        fbe = fb(esk(i,kloc,jspin)-efermi)
                        if (lentro) entropy = entropy - wt*
     $                       (fbe*log(fbe)+(one-fbe)*log(one-fbe))
                           esum = esum + ewt*fbe
                     end if
                  end if
               end do
            end do
         end do
      end if
      if(multiprc) then
         call MPI_REDUCE (esum,esumg,1,MPI_DOUBLE_PRECISION,
     $        MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE (entropy,entropyg,1,MPI_DOUBLE_PRECISION,
     $        MPI_SUM,0,MPI_COMM_WORLD,ierr)
         if(master) then
            esum = esumg
            entropy = entropyg
         end if
      end if
c
c     Sum up the rest of the variables:
c
      
      if (jsx.eq.2) then
c
c        Determine the spin up and spin down population
c
c        Count:
c
         do jspin = 1,jsx
c
c           Initialize
c
            sptemp = zero
c
c           Scan k-points and bands.  Note that we only get here
c            when jsx=2, so sfac = 1 and is omitted.
c
c           Version 1.33: again with the local/global differentiation:
c
            kloc = 0
c
            do k = istart,ifin
c
               kloc = kloc + 1
c
               do i = 1,nvx
                  if(esk(i,kloc,jspin).lt.efm) then
                     sptemp = sptemp + weight(k)
                  else if(esk(i,kloc,jspin).lt.efp) then
                     sptemp = sptemp +
     $                    weight(k)*fb(esk(i,kloc,jspin)-efermi)
                  end if
               end do
            end do
c
c           Note that sptemp has to be output inside the jspin loop.
c
            if(multiprc) then
               call MPI_REDUCE (sptemp,sptempg,1,
     $              MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
               if(master) sptemp = sptempg
            end if
c
c           Ideally this would refer to an array, but the old
c            Absoft compiler has serious dimension problems when
c            initializing by data-statement, and I don't think
c            you can initialize an array by parameter statement
c
c           This part only needs to be done by the master processor:
c
            if(master) then
               if (jspin.eq.1) then
                  write(15,'(//''The system has '',f10.5,
     $                 '' spin-up electrons'')') sptemp
                  smag = sptemp
               else
                  write(15,'(''and '',f10.5, '' spin-down electrons'')')
     $                 sptemp
                  smag = smag - sptemp
                  write(15,'(''for a total magnetization of '',
     $                 f10.5,'' electrons per unit cell.''//)') smag
               end if
            end if
         end do
      else
c
c        No magnetic moment:
c
         smag = zero
      end if
c
      if(master) then
         write(15,*) 'Temperature = ',tkb
         write(15,*) 'EFermi = ',efermi
         write(15,*) 'E(T) = ',esum
         if(lentro) then
            write(15,*) 'S(T) = ',entropy
c
c           Calculate the free energy:
c
            free = esum - tkb*entropy
            write(15,*) 'A(T) = ',free
c
c           Extrapolate
c
            earg1 = half*(esum+free)
            write(15,*) 'E(0) = ',earg1,' + O[T**4]'
c
c           If the pressure calculation is turned off, than
c            earg2 will be E(T):
c
            earg2 = esum
         else
            earg1 = esum
         end if
c
         if (lpress) then
c
c           Fermi level shift:
c
c           In an insulator, dfbsum will be zero, as will dferm.
c            Then dfermi is undetermined, but in an insulator
c            there will be no contribution from this term anyway:
c
            if(abs(dfbsum).gt.small) then
               dfermi = dferm/dfbsum
            else
               dfermi = zero
            end if
c
c           Pressure is - dE/dV:
c
            press = -(desum1 + (desum2 - dfermi*desum3))
c
            write(15,*) 'Fermi level shift: ',dfermi
            write(15,*) 'P(T) = ',press
c
            if(lentro) then
c
c              Calculate the entropy shift with pressure.  Remember
c               that pressure has a (-) sign.  Actually, this is
c               -T S'(V)
c
               pent = dfermi*dent2 - dent1
c
               write(15,*) 'Entropy shift = ',pent
c
c              free energy pressure:
c
               fpress = press - pent
               write(15,*) 'P_A(T) = ',fpress
c
c              extrapolated pressure:
c
               earg2 = half*(press+fpress)
               write(15,*) 'P(0) = ',earg2,' + O[T**4]'
            else
c
c              The second output is the pressure:
c
               earg2 = press
            end if
         end if
      end if
c-----------------------------------------------------------------------
c
c     When lspdir is true, add up the contributions to the
c      magnetization, find the spin direction and magnitude, and print
c      out everything to the file SPNDIR (unit 42).
c
      if(lspdir) then
c
c        Zero the accumulation arrays:
c
         do i = 1,natoms
            do j = 1,3
               spdirl(j,i) = zero
            end do
         end do
c
c        Now add up the contribution from each eigenvalue at each
c         k-point, along with the appropriate weights.  Note that there
c         is no sum over spins, as this calculation can only be done
c         in the non-collinear spin case, when jspin = 1.
c
c        Version 1.33: one more time with the local k-points:
c
         kloc = 0
c
         do k = istart,ifin
c
            kloc = kloc + 1
c
            wt = sfac*weight(k)
            do i = 1,nvx
               if(esk(i,kloc,1).lt.efp) then
                  fbe = wt
                  if(esk(i,kloc,1).gt.efm) fbe =
     $                 wt*fb(esk(i,kloc,1)-efermi)
                  do j = 1,natoms
                     do l = 1,3
                        spdirl(l,j) = spdirl(l,j) + fbe*spev(l,j,i,kloc)
ctemp
c$$$                        write(*,*) l,j,i,k,spev(l,j,i,kloc),fbe,
c$$$     $                       spdirl(l,j)
cend temp
                     end do
                  end do
               end if
            end do
         end do
c
c        Accumulate everything onto the master processor:
c
         if (multiprc) then
            call MPI_REDUCE (spdirl,spdir,spsize,MPI_DOUBLE_PRECISION,
     $           MPI_SUM,0,MPI_COMM_WORLD,ierr)
         else
            do i = 1,natoms
               do j = 1,3
                  spdir(j,i) = spdirl(j,i)
ctemp
c$$$                  write(*,*) j,i,spdirl(j,i)
cend temp
               end do
            end do
         end if
c
c        Write out the three spin components, the magnitude of the
c         vector, and the (theta,phi) direction of the vector
c         (in degrees)
c
         if(master) then
            do i = 1,natoms
               smagl = sqrt(spdir(1,i)*spdir(1,i)+
     $              spdir(2,i)*spdir(2,i)+spdir(3,i)*spdir(3,i))
               theta = rad2deg*acos(spdir(3,i)/smagl)
c
c              I've always wanted to use atan2 for SOMETHING:
c
               phi = rad2deg*atan2(spdir(2,i),spdir(1,i))
c
               write(42,'(a20,i5,4f10.6,2f12.6)')
     $              label,i,(spdir(j,i),j=1,3),
     $              smagl,theta,phi
            end do
         end if
      end if
c-----------------------------------------------------------------------
c     
c
c     Now, if one of the "lq*" modes is on, and if this is the master root,
c      collect the eigenvalues and write them out in the appropriate format:
c
      if(master.and.(mode.gt.4)) then
c
         do jspin = 1,jsx
            do i = 1,MPISIZE
c
c              Open the appropriate file.  Note that we need to get
c               the indexing right:
c
               write(qstring,'(''qeigtemp.'',i1,i6.6)') jspin,i-1
               open(unit=55,file=qstring,status='old',
     $              form='unformatted')
c
c              There should be count(i)/mh k-points in this file
c
c
c              For mode=5 or 6, we'll have the s, p, and d decomposition
c               of the eigenstates.  For mode=10 we'll have the orbital
c               by orbital decomposition, and for 9 we'll have no
c               decomposition at all.
c
c              For mode 5, 6, and 10, the line with the occupation
c               is written (nbas,orbocc(1),orbocc(2),...) where nbas is the
c               number of entries, and orbocc holds the occupation
c               of an orbital or group of orbitals.
c
c              Note that, as it says above,
c
c$$$      if(MPIRANK.lt.krem) then
c$$$         istart = MPIRANK*(kproc+1)+1
c$$$         ifin = istart+kproc
c$$$      else
c$$$         istart = MPIRANK*kproc+krem+1
c$$$         ifin = istart+kproc-1
c$$$      end if
c
c              whence the needed change seems obvious
c
c              Noting that i is the processor number:
c
c              Remembering that MPIRANK goes from 0 to MPISIZE-1:
c
               if(i.lt.krem+1) then
                  jkpts = kproc+1
               else
                  jkpts = kproc
               end if
c
ctemp
c$$$               write(*,*) count(i)/mh,' k-points for file ',i
c$$$               write(*,*) jkpts,' k-points for file ',i
cend temp
c$$$               do j = 1,count(i)/mh
               do j = 1,jkpts
c
c                 what we read depends upon the mode:
c
                  if(mode.eq.5) then
c
c                    Everything remains unformatted.  This is just
c                     a big copy routine.
c
 43   format(3f10.5,2X,f10.5,i10,d20.10,' K')
                     read(55) ak(1),ak(2),ak(3),natm,ntot,wt0
                     write(45,43) ak(1),ak(2),ak(3),wt0,ntot,wt0
                     do k = 1,ntot
                        read(55) eign
                        write(45) eign,zero
                        do l = 1,natm
                           read(55) norb,s,p,d
c
c                          norb had better be 3
c
                           write(45) s,p,d,zero
                        end do
                     end do
                  else if(mode.eq.6) then
                     read(55) ak(1),ak(2),ak(3),natm,ntot,wt0
                     write(45,'(3f10.6,2X,f10.5,i10,d20.10,'' K'')')
     $                    ak(1),ak(2),ak(3),wt0,ntot,wt0
                     do k = 1,ntot
                        read(55) eign
c
c                    If the eigenvalue is large then change the formatting a bit:
c
                        if((eign.le.-ten).or.(eign.ge.hundred)) then
                           write(45,'(f10.6,e15.6)') eign,zero
                        else
                           write(45,'(f10.7,e15.6)') eign,zero
                        end if
                        do l = 1,natm
c
c                          norb had better be 3
c
                           read(55) norb,s,p,d
                           write(45,'(4f15.9)') s,p,d,zero
                        end do
                     end do
                  else if(mode.eq.10) then
                     read(55) ak(1),ak(2),ak(3),natm,ntot,wt0
                     write(45,'(3f10.6,2X,f10.5,i10,d20.10,'' K'')')
     $                    ak(1),ak(2),ak(3),wt0,ntot,wt0
                     do k = 1,ntot
                        read(55) eign
                        write(45,'(f10.7,e15.6)') eign,zero
                        do l = 1,natm
                           read(55) norb,(orbocc(k1),k1=1,norb)
c
c                          Set up the format string
c
                           write(kforb,'(''('',i2.2,''f15.9)'')') norb
                           write(45,kforb) (orbocc(k1),k1=1,norb)
                        end do
                     end do
                  else if(mode.eq.9) then
                     read(55) ak(1),ak(2),ak(3),ntot,wt0,
     $                    (en(l),l=1,ntot)
                     write(45,'(3f10.6,2X,f10.5,i10,d20.10,'' K'')')
     $                    ak(1),ak(2),ak(3),wt0,ntot,wt0
                     write(45,'(f12.7,e15.6)')
     $                    (en(l),zero,l=1,ntot)
                  else
c
c                    Must be the QAPW mode:
c
                     read(55) ntot,ak(1),ak(2),ak(3),(en(l),l=1,ntot)
                     write(45,'(1x,3f9.6,10x,4f9.5/(38x,4f9.5))')
     $                    ak(1),ak(2),ak(3),(en(l),l=1,ntot)
                  end if
               end do
c
c              Close and delete the temporary file:
c
               close(55,status='delete')
            end do
         end do
      end if
      return
      end
