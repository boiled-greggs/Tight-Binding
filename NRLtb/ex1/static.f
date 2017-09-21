      PROGRAM static
c     
c     Author:Papaconstantopoulos
c     Modifications by Cohen and Mehl 
c     Performs energy and eigenvalue evaluations of S-K fit
c     Two center approximation
c     with d-orbitals
c     
c=======================================================================
c     
c     This software and any accompanying documentation are released "as
c     is." The U.S. Government makes no warranty of any kind, expressed
c     or implied, concerning this software and any accompanying
c     documentation, including without limitation, any warranties of
c     merchantability or fitness for a particular purpose. In no event
c     will the U.S. Government be liable for any damages, including any
c     lost profits, lost savings, or other incidental or consequential
c     damages arising out of the use, or inability of use, of this
c     software or any accompanying documentation, even if informed in
c     advance of the possibility of such damages.
c     
c=======================================================================
c     
c     REVISION HISTORY:
c     
c-----------------------------------------------------------------------
c     
c     Modified by REC 1/20/93 to simultaneously fit multiple structures
c     
c-----------------------------------------------------------------------
c     
c     This now is an "energy evaluation package", i.e., it does no fitting
c     
c     Current supported modes:
c     
c     3  -- Total Energy, no Pressure
c     4  -- Total Energy and Pressure
c     6  -- Total Energy, no Pressure, and produces an LAPW like QLMT file
c     7  -- Total Energy, dE/d(c/a) (Omitted)
c     8  -- Total Energy, and no Pressure, and produces an APW like
c             eigenvalue file (called QAPW).  (added 9/12/96)
c     
c     4 is the default mode
c     
c                                                    --  mjm 17 Oct 1994
c-----------------------------------------------------------------------
c     
c     Added possible calls to "asymp", which reads the asymptotic
c      values of the onsite parameters and decides how to interpolate
c     
c                                                    --  mjm  5 Dec 1994
c-----------------------------------------------------------------------
c     
c     Eliminated the "jstruc" dependent variables.  Now we can do
c      as many structures as we please in a single run.
c                                                    --  mjm 15 Dec 1994
c-----------------------------------------------------------------------
c     
c     See setup:  split onsite d terms, quadratic S-K parameters
c                                                    --  mjm 31 Oct 1994
c-----------------------------------------------------------------------
c     
c     Since we only need one triangle of the Hamiltonian and overlap
c      matrix, only define one such triangle.
c                                                    --  mjm  3 Feb 1997
c-----------------------------------------------------------------------
c     
c     Eliminated the "asymptotic" on-site parameters, since we never
c      use them anyway.
c     
c     Rearranged the dimensions so that the user can use the same
c      executable for single and multiple species parameter files
c                                                   --  mjm 23 Sept 1997
c-----------------------------------------------------------------------
c     
c     Add the option to constrain like atom overlap matricies to have
c      the correct behavior as R -> 0, i.e.
c      S_{l,l',m}(R) -> delta_{l,l'}
c     
c     For these terms (and only these terms), the parametrization
c      is of the form
c     
c     S(R) = [delta_{l,l'} + R( A + B R + C R^2)] Exp[-D^2 R]
c     
c     This is signaled from the parameter file if the first parameter
c      (which is always evaluated as lamda^2) is negative.
c      The flag for this is the logical variable realov in common
c      block overlap.  See "input.f and setup.f" for details.
c                                                   --  mjm 26 Sept 1997
c-----------------------------------------------------------------------
c     
c     Added mode 9:  QLMT like file with "0" atoms
c                                                   --  mjm 28 Oct  1997
c-----------------------------------------------------------------------
c     
c     This is a port of the serial code for an f77/MPI based
c      multiprocessor system.  Note that the number of processors is
c      arbitrary, and that the true parallelization takes place only
c      in the bande.f subroutine, where the procedure is assign a
c      set of k-points to each processor.  Given the setup times
c      involved, this is probably highly scalable until the number
c      of processors gets to about 1/2 the number of k-points.
c     
c     When possible, MPI commands will be set apart from the rest of
c      the code by "cMPI/cend MPI" comments.
c                                                    -- mjm  28 Jan 1998
c=======================================================================
c     
c     Version 1.05:
c     
c     This version integrates the serial and parallel codes.  Single
c      processor routines call a set of fake MPI libraries, and bande.f
c      uses a slightly different algorithm, avoiding massive data
c      transfer, when there is only one processor in the system.  In
c      addition, the form of the parameter file has been modified to
c      make it more compatible with Florian Kirchhoff's TBMD routines,
c      and to, in general, move parameter dependent information from
c      the SKIN file into the parameter file.
c                                                    -- mjm   7 Apr 1998
c-----------------------------------------------------------------------
c     
c     Continuing with the above, move the JSPIN and JSOVER parameters
c      into the parameter file, cleverly hidden in line one as
c      "NN" or something like that (see input.f for details)
c                                                    -- mjm  20 Apr 1998
c=======================================================================
c     
c     Version 1.06:
c     
c     jsover has been changed to a logical variable:
c        jsover = .true.  -> Non-orthogonal Hamiltonian (S <> identity)
c        jsover = .false. -> Orthogonal Hamiltonian (S = identity)
c                                                   -- mjm  6 July  1998
c-----------------------------------------------------------------------
c     Separated the eigenvector calculation parameters into a parameter
c      file 'P2'.
c     
c     Moved k-point generation parameters into 'P3'
c     
c                                                   -- mjm  6 July  1998
c-----------------------------------------------------------------------
c     
c     Got rid of the useless "iwrite" variable in /lia/
c
c     Renamed /overlap/ to /parcom/, now it exports the parameterization
c      type nltype.
c
c     Added a new nltype, 90000, for use in the parametrization of the
c      H atom.  See setup.f
c
c                                                   -- mjm  4 Aug   1998
c-----------------------------------------------------------------------
c     
c     The fourth field in the first line of a QLMT file is supposed to
c      indicate the number of representative atoms in the structure.
c      In our case this is just the number of atoms in the file.
c      Correct the input to reflect this.
c                                                   -- mjm 18 Sept  1998
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
c     Mode 1 -- Prints the extrapolated E(0) and E(T) in the SKENG file
c     Mode 2 -- Prints extrapolated E(0) and P(0) in the SKENG file
c     
c     The other modes remain the same.
c     
c     Change the default mode to 1.
c     
c     Eliminated the logical variables lpress,lqlmt,lqlmt0,lqapw, which
c      tell bande.f how to handle the modes.  Now the mode number is
c      passed to bande.f, which will figure out what to do on its own.
c      (Note that this requires changes to diag.f as well.)
c     
c     Eliminated the common block
c     
c     common /fermcom/ thisfe,bandes,press
c     
c     which passed the main results of bande.f back to the main program.
c      Note that the values "bandes" and "press" in bande.f have been
c      replaced by generic arguments, as they may represent different
c      quantities in different modes.
c     
c                                                    -- mjm  26 Jan 1999
c-----------------------------------------------------------------------
c     
c     Mode 5 added.  This is identical to mode 6 except that the QLMT
c      file is written unformatted.  This increases the accuracy of the
c      storage of k-points and eigenvalues, making it possible to
c      accurately evaluate the extrapolated E(0) from an arbitrary
c      temperature.
c                                                    -- mjm  11 Feb 1999
c=======================================================================
c     
c     Version 1.11:
c     
c     Broke up setup.f into two parts.  The initialization of the
c      Slater-Koster parameters and on-site terms is now in setpar.f,
c      which has been moved outside the k-point loop in bande.f.  It
c      could, with minor changes, be called from static.f, just below
c      the call to search2.  This may be useful if we ever parallelize
c      search2.
c
c     The k-point part of the Hamiltonian setup now takes place in the
c      routine setham.f.  This is called from bande.f inside the k-point
c      loop and calls rotate.f to set up the Hamiltonian matrix.
c                                                    -- mjm  30 Apr 1999
c-----------------------------------------------------------------------
c     
c     Broke setvol.f into setparv.f and sethamv.f, corresponding to
c      setpar.f and setham.f for setup.f.
c     Print the version number at the head of standard output and
c      SKOUT.                                        -- mjm  27 Jul 1999
c-----------------------------------------------------------------------
c     
c     search2.f now outputs the minimum pair bond length, dismin,
c      so that we can limit write statements to SKOUT to the "master"
c      processor.                                    -- mjm   4 Aug 1999
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
c      the input changes.  (Note -- as of 12 Jan 2000, the interpretation
c      of "density" changes.  See setpar.f for details.)
c     
c     To help in this effort, the variable "jspins" has been moved from
c      to common block /lia/ and moved into the calling sequence of
c      input.f and bande.f.
c     
c     Also note the new logical array, "flipos".  "flipos(i)" is false
c      in the default case, i.e., the spin of atom i points in the spin
c      direction we agree to be called "up."  If flipos(i) is false,
c      atom i is "pointing down," and so we flip its onsite spins.  See
c      "input1.f" and "setpar.f" to see how this happens.
c     
c     WARNINGS:
c     
c     Note that input1.f defaults "flipos(i) = .false." for all atoms.
c      Thus for unpolarized or "standard" calculations we follow the
c      default magnetic case.  Because of this, it is important that
c      FLIPOS IS CHANGED ONLY IN INPUT1.F.
c     
c     Also note that if two atoms have opposite spins, THEY ARE NO
c      LONGER EQUIVALENT.  That is, the space group of the system is now
c      different (and lower) than in the all-spin up or unpolarized case.
c     
c     The entire procedure only works if the "hopping" and "overlap"
c      parameters remain the same for up-up, up-down, and down-down
c      atoms, with only the spin-up and spin-down on-site terms
c      differing.
c                                                     -- mjm 21 Dec 1999
c-----------------------------------------------------------------------
c     
c     bande.f now reports back the magnetic moment, if any, in the
c      real*8 variable "smag".  If the system is magnetic this will be
c      printed in the SKENG file.
c                                                     -- mjm  8 Feb 2000
c-----------------------------------------------------------------------
c     
c     Let A-A, B-B, C-C, etc. interactions have Harrison's canonical
c      sign.  i.e.:
c     
c                 H    S
c     
c     ss sigma    -    +
c     sp sigma    +    -
c     pp sigma    +    -
c     pp pi       -    +
c     sd sigma    -    +
c     pd sigma    -    +
c     pd pi       +    -
c     dd sigma    -    +
c     dd pi       +    -
c     dd delta    Arbitrary (Harrison has zero)
c     
c     This is signaled by the logical parameter "harrison", stored in
c      common block parcom.
c                                                     -- mjm  7 Mar 2000
c=======================================================================
c     
c     Version 1.21:
c     
c     When running static, if a file named TWEAKS exists in the same
c      directory as SKIN, then this file is examined for certain
c      keywords, as outlined in the subroutine tweaks.f.  These keywords
c      set certain output options, such as
c     
c     (1) printing out values of the on-site, hopping, and overlap
c      matrix elements;
c     
c     (2) printing out the Hamiltonian and Overlap matrices at certain
c      k-points.
c     
c     The subroutine tweaks.f determines if the file TWEAKS exists, and,
c      if so, whether it does anything or not.
c                                                  -- mjm  7 August 2000
c=======================================================================
c     
c     Version 1.22:
c     
c     Fixed static.f so that the code will properly run in parallel
c      under MPI.  Also removed those annoying "cMPI/cend MPI" comments
c      in favor of proper indenting.
c     
c     Sustantially rewrote bande.f and fermi.f so that eigenvalues and
c      their derivatives do not have be schleped across processors.
c
c                                                 -- mjm 19 October 2000
c=======================================================================
c
c     Version 1.30:
c
c     Impliment non-collinear magnetization for spin polarized parameter
c      sets, as outlined by Pickett (1996).
c
c     The subroutine input1.f provides a logical parameter, lnoncol.
c      If true, non-collinear magnetization is enabled.  Note that we
c      must have a spin-polarized set of parameters, which is signaled
c      by input.f returning a value jspins=2.  If jspins=1 and
c      lnoncol=.true. the program will print out an error message and
c      continue.
c     If lnoncol=.true. then input1 also exports a spin direction
c      unit vector, spinat(i=1,2,3;iatom), for each atom.
c                                                  -- mjm  9 August 2000
c
c     Implement total energy calculations for arbitrary spin directions.
c                                                  -- mjm 11 August 2000
c
c     Print out spin and angular momentum decomposed eigenvector
c      information in the QLMT file.
c
c     Implement pressure calculations for arbitrary spin directions.
c                                                  -- mjm 14 August 2000
c
c     Use perturbation theory to calculate the derivative of the energy
c      with respect to the spin directions.  That is, we can determine
c      to torque on each atom.  Note that for this to work properly,
c      each atom must be considered to be independent when it comes
c      to calculating symmetry.  In general this will mean a large
c      increase in the number of k-points.
c     Derivative information will be stored in the file SKDERIV, unit
c      number 41.
c                                                  -- mjm 14 August 2000
c-----------------------------------------------------------------------
c
c     If bande returns with a non-zero value for iberr, a major error
c      has occured, so stop the calculation.
c                                                  -- mjm  9 August 2000
c=======================================================================
c
c     Version 1.31:
c
c     This is a merger of the 1.21->1.30 line of development with the
c      revised 1.22 version.
c                                                 -- mjm 23 October 2000
c
c     We specify the direction of the spin of each atom in the input.
c      But this does not guarantee that we know the direction of each
c      spin in the output.  Here we turn on an option to determine the
c      "output spin direction" of each atom.  Note that at
c      self-consistency the input and output spin directions should be
c      aligned, so this is one way of determining the directions of
c      spins at a domain wall, or in a frustrated anti-ferromagnet.
c      See input1.f for directions on how to turn on this calculation.
c      Results are in the file SPNDIR (unit 42).
c                                                -- mjm  26 October 2000
c     Updated "tweaks" to provide more information on how the program
c      is progressing through the k-point list.
c                                                  -- mjm 28 August 2001
c=======================================================================
c
c     Version 1.32:
c
c     Merge g77 (with spin directions) and Absoft (with tweaks) versions
c      of 1.31
c
c     Add a new mode (10) which prints out orbital-by-orbital
c      decompostion of the wave function in the QLMT file.
c
c     This necessitates a small change in the format of the line in
c      static which reads the mode.  Previously it was:
c
c              read(10,10) mode
c      10      format(5x,i1)
c
c     We change it to
c
c     	 read(10,'(5x,i)') mode
c
c     In Absoft and g77, this returns "mode=10" as 10, "mode= 1" or
c      "mode=1" as 1, etc.  This may vary in other compilers.
c
c                                                      -- mjm  1 November 2001
c
c     Add /fazon/ common block, to make sure that the nearly
c      ready for depreciation variables posq are propogated to
c      all nodes
c
c                                                    -- mjm  13 May 2002
c
c=======================================================================
c
c     Version 1.33:
c
c     Correctly parse the "Mode=" on the first line of the SKIN file.
c      This is accomplished via a new "modeparse" subroutine
c
c                                                    -- mjm  14 Aug 2002
c
c     Change bande and its associated files to only allocate eigenvalue
c      arrays for the number of k-points which will be used on one
c      processor, rather than dimensioning for all k-points on each
c      processor.  Note that this will require recompiling the code
c      whenever we DECREASE the number of processors in the calculation.
c      The relevant parameter is "minproc", which is dimensioned in P3.
c      Note that static.f will check to see if minproc is consistent
c      with the number of processors in the code, and abort if the
c      dimension needs to be fixed.
c
c                                                    -- mjm  15 Aug 2002
c
c=======================================================================
c     
c     Current supported modes:
c     
c     1 -- Prints Gillan's extrapolated E(0) and E(T) in the SKENG file
c     2 -- Prints extrapolated E(0) and P(0) in the SKENG file
c     3 -- Total Energy, no Pressure
c     4 -- Total Energy and Pressure
c     5 -- Mode 6 with QLMT output in unformatted file
c     6 -- Total Energy, no Pressure, and produces an LAPW like QLMT file
c     7 -- Total Energy, dE/d(c/a) (Omitted)
c     8 -- Total Energy, and no Pressure, and produces an APW like
c             eigenvalue file (called QAPW).  (added 9/12/96)
c     9 -- QLMT like file with "0" atoms
c    10 -- Mode 6, but we print out the occupation of each orbital
c             rather than p = x+y+z and d = yz+zx+xy+x^2-y^2+3z^2-r^2
c     
c     The default mode is now "1".  (26 Jan 1999)
c     
c=======================================================================
c     
      implicit real*8 (a-h,o-z)
c     MPI
      include 'mpif.h'
c=======================================================================
c     
c     Version Number
c     
      integer v_main,v_sub
      parameter (v_main =  1)
      parameter (v_sub  = 33)
c=======================================================================
c     
      logical master,multiprc
c
      common /MPIcomm/ MPISIZE,MPIRANK,master,multiprc
c
      include 'P1'
      include 'P2'
      include 'P3'
      common /types/ npkind(mptype)
      common /lattyp/ avec(3,3),bvec(3,3),bij(3,3),wsvol,lattic
      common /kmesh1/ qkpt(3,nkd,nwdd),wghtk(nkd,nwdd),nkv(nwdd)
      common/codes/posn(matom,3),kkind(matom),natoms
      common /parinfo/ valence(mkind),kinds,kbas(mkind)
      common /cutcom/ dnnsq(mint),dnn(mint),
     $     screenl(mint),scrninv(mint)
      common/latt1/plv(3,3)
      common/struc/ nv
      common/cd1/x(nkd),y(nkd),z(nkd)
      common/params/ param(npard,mspin)
      common /kstuff/ ksk(matom)
      logical jsover
      common/lia/ npts,jsover
c
c     posq stores relative phases for atomic wave functions.
c      Normally these are set to 0 (the "0 0 0" you see in the
c      SKIN file), but it's possible to make a choice of phase so
c      that, e.g., you can construct real Hamiltonian for the
c      diamond structure.  Almost useless in complex structures.
c      We may depreciate this in the next version of the code.
c      The common block /fazon/ allows for the communication of
c      the values of posq from input1, where it is generated,
c      to rotate, where it is used.
c
      common /fazon/ posq(matom,3)
c
c     spinat(i,j) is the ith Cartesian component of the spin direction
c      of the jth atom, expressed as a 3 component unit vector.  This
c      is defined in input1.f.
c
      real*8 spinat(3,matom)
c     
c     Internally, "lpress" is useful to indicate when we are
c      calculating pressure
c     
      logical lpress
c     
c     Temporary storage
c     
      real*8 temp(3)
c     
      real*8 zfill,tkb
c     
c     Flag for using standard or restricted like-atom overlap
c      parameters
c     
      logical realov,harrison
      common /parcom/ nltype,realov,harrison
c     
c     These are the "spin flip" parameters (see discussion in the
c      Version 1.20 revision notes, above).  Note that input1.f will by
c      default set all of these false.
c     
      logical flipos(matom)
c     
c
c     Version 1.30 adds non-collinear spins to a spin-polarized
c      calculation, signaled by lnoncol=.true. coming from
c      input1.f
c
      logical lnoncol
c
c     Also in 1.30, calculate the "torque" on each atom if desired.
c      Note that the value of ltorque is meaningless if lnoncol
c      is false.
c
c     In version 1.31 and above, lspdir is true if we want to
c      calculate the magnitude and direction of the OUTPUT spin
c      on each atom.  Note that lnoncol must be true if this is
c      turned on, and that, in principle (but not yet in practice)
c      ltorque and lspdir can be computed in the same run.
c
      logical ltorque,lspdir
c
c     dspth = d spinat/d theta, and dspphi = d spinat/d phi,
c      both in units of energy/degree.
      real*8 dspth(3,matom),dspphi(3,matom)
c
c     This is a trick:  sometimes we may want to define mptype=2,
c      if we're only covering monatomic systems and want to save
c      space
c     
      integer npkint(3)
c     
c     This is an array of search information passed from input1 to
c      search2 through this common block:
c     
      integer nkdarry(6)
      common/disply/ nkdarry
c     
c     lp_par is true if setpar.f is to print out the on-site, hopping,
c      and overlap Slater-Koster parameters.
c     
c     lp_parv is true if setparv.f is to print out the derivative
c      on-site, hopping, and overlap Slater-Koster parameters.
c
c     lp_klst is true if we are to print the current k-point onto
c      standard output
c     
      logical lp_par,lp_parv,lp_klst
c     
c-----------------------------------------------------------------------
c     
c     Character strings:
c     
c     Title for the SKOUT file:
c
      character*70 title
c
c     Identifier for the SKENG file:
c
      character*20 label
c
c     Input mode string:
c
      character*2 modest
c     
c=======================================================================
c     
      parameter (zero = 0d0)
c$$$      parameter (one  = 1d0)
c$$$      parameter (two  = 2d0)
c$$$      parameter (four = 4d0)
c     
c$$$      parameter (pi = 3.14159265358979323846d0)
c$$$      parameter (tpi = two*pi)
c$$$      parameter (fpi = four*pi)
c=======================================================================
c     
c     Set for splitting d on-site into t2g and eg
c     
      data npkint/4,10,14/
c=======================================================================
c     
c     Start up MPI:
c     
      call MPI_INIT (ierr)
c     
c     Find the number of processors and the rank of the current processor
c     
c     Note that we should essentially recover the scalar version
c      of this code by setting MPISIZE=1 and MPIRANK=0
c     
      call MPI_COMM_SIZE (MPI_COMM_WORLD, MPISIZE , ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD, MPIRANK , ierr)
c     
c     This is the master processor if
c     
      master = MPIRANK.eq.0
c
c     This is a multi processor calculation if
c
      multiprc = MPISIZE.gt.1
c     
c     Note that master will be true on a single processor machine.
c      See mpifake.f MPI_COMM_RANK for details.
c     
c     Print out the version number on standard output
c     
      if (master) write(*,'(''static Version '',i2,''.'',i2.2)')
     $     v_main,v_sub
c
c     New in version 1.33:  check that MPISIZE is consistent with
c      minproc:
c
      if(minproc.gt.MPISIZE) then
c
         write(0,*) 'minproc dimensioned to ',minproc
         write(0,*) 'This calculation is running with ',
     $        MPISIZE,' processors'
         write(0,*) 'minproc must be less than or equal to ',
     $        'the number of processors used'
         write(0,*) 'Redimension minproc or change the ',
     $        'number of processors'
c
         stop 'minproc > MPISIZE'
      end if
c     
c     This little ditty was originally in input0.f, but we're going
c      to eliminate that particular subroutine
c     
c     Define npkind.  1 is the number of onsite terms, 2 the number
c      of hopping terms for like elements, 3 the number of hopping
c      terms for unlike elements.
c      Assume s,p,d two center parametrization
c     
cMPInote
c     This is done on all processors, so we need not export npkind
cend MPInote
c     
      do i = 1,mptype
         npkind(i) = npkint(i)
      end do
c     
c     Input and output will only take place on processor zero:
c     
      if (master) then
c     
ctemp
c$$$         write(*,*) meigen,mpress,mh,mhe,mhp
c$$$         write(0,*) matom,mkind,mh,nkd,mpair
cend temp
c     
         open (10,file='SKIN',status='unknown',blank='zero')
         open (15,file='SKOUT',status='unknown',blank='zero')
c     
c     Print out the version number in SKOUT:
c     
         if (master) write(15,'(''static Version '',i2,''.'',i2.2)')
     $        v_main,v_sub
c-----------------------------------------------------------------------
c     
c     Version 1.21 update:
c     
c     If the file TWEAKS exists in the same directory, search it for
c     keywords.  Note that this will have to be done on the master
c     node in MPI.
c     
         call tweaks (lp_par,lp_parv,lp_klst,kptpr)
c     
c-----------------------------------------------------------------------
c     
c$$$         read(10,10) mode
c$$$ 10      format(5x,i1)
c
c        We hope that this will work:
c
         read(10,'(5x,a2)') modest
c
         call modeparse(modest,mode)
c
         write(15,*) 'mode=', mode
c
      end if
c     
c-----------------------------------------------------------------------
c     
c     Version 1.22 fix of version 1.21 update: MPI_BCAST calls MUST be
c     made on all processors.  At least that's all I know how to do.
c     
c     Broadcast these parameters to the world at large:
c     
      call MPI_BCAST(lp_par,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lp_parv,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lp_klst,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(kptpr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mode,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c     
c-----------------------------------------------------------------------
c     
c     lpress is still useful within static for evaluation of various flags
c     
      lpress = (mode.eq.2).or.(mode.eq.4)
      if(mode.eq.7) then
c
c     Information written only by processor 0:
c     
         if(master) then
c
            write(*,*) 'Calculating derivative wrt c/a'
            write(*,*) 'mode not implemented'
            write(0,*) 'no c/a derivatives yet'
         end if
         go to 10000
      end if
c     
c     Note that we've already eliminated mode 7, so this works:
c     
      if(mode.gt.4) then
         if(master) then
            write(*,21)
 21         format(/'warning:  If you have multiple structures,'/
     $           'this routine will concatenate all of your'/
     $           'eigenvalues into one file.'/)
         end if
      end if
c     
c     More stuff best done only on processor 0:
c     
      if(master) then
         if(lpress)
     $        write(*,*) 'no pressure calculated on this run'
         if(mode.eq.5) then
c     
c        Unformatted file:
c     
            write(*,*) 'producing an unformatted QLMT file'
            open(unit=45,file='QLMT',status='unknown',
     $           form='unformatted')
         else if(mode.eq.6) then
            write(*,*) 'producing an ASCII QLMT file'
            open(unit=45,file='QLMT',status='unknown')
         else if(mode.eq.8) then
            write(*,*) 'producing an eigenvalue file called QAPW'
            open(unit=45,file='QAPW',status='unknown')
         else if(mode.eq.9) then
            write(*,*) 'producing a QLMT-like file'
            open(unit=45,file='QLMT',status='unknown')
         else if(mode.eq.10) then
            write(*,*) 'producing an ASCII QLMT file'
            write(*,*) 'with orbital occupations'
            open(unit=45,file='QLMT',status='unknown')
         end if
c     
         open(unit=31,file='SKENG',status='unknown',blank='zero')
      end if
c     
c     check to see if we have the right dimensions for various
c      scenarios:
c     
      if(lpress.and.(mpress.ne.1)) then
         if(master) then
            write(*,*) 'set mpress = 1 in P2 for pressure calculations'
            write(0,*) 'lpress and mpress are inconsistent'
         end if
         go to 10000
      end if
c
      if(lpress.and.(meigen.ne.1)) then
         if(master) then
            write(*,*) 'set meigen = 1 in P2 for pressure calculations'
            write(0,*) 'mode and meigen are inconsistent'
         end if
         go to 10000
      end if
c
c     If we are decomposing the eigenvalues we need meigen to be
c      turned on:
c
      if(((mode.eq.5).or.(mode.eq.6).or.(mode.eq.10))
     $     .and.(meigen.ne.1)) then
         if(master) then
            write(*,*) 'set meigen = 1 in P2 for qlmt calculations'
            write(0,*) 'mode=',mode,' and meigen<>1'
         end if
         go to 10000
      end if
c
      if(master) then
         read(10,*) tkb,eigcut
         write(15,*) ' temp. broad. (ryd):',tkb,
     $        ' eigenvalue cutoff = ',eigcut
      end if
c
c     Continue reading input on processor 0
c
      if(master) then
c
c        Read in tight-binding parameters
c
         call input (nparam,jspins)
c
      end if
c 
c     All processors need jspins and eigcut.
c      In Versions 1.22 (except 1.30) and up they also need tkb.
c
      call MPI_BCAST(jspins,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(eigcut,1,MPI_DOUBLE_PRECISION,
     $     0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tkb,1,MPI_DOUBLE_PRECISION,
     $     0,MPI_COMM_WORLD,ierr)
c
c     nparam, the number of parameters; 
c     jsover, kinds, and kbas are needed by all processors
c
      call MPI_BCAST(nparam,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(jsover,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(kinds,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(kbas,mkind,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c
      if(kinds.gt.mkind) then
         if(master) then
            write(0,*) 'Expected <',mkind+1,' types of atom, got',
     $           kinds
            write(0,*) 'static:  kinds > mkind'
         end if
         go to 10000
      end if
c
c     now we know how many kinds of atoms we have.  this had
c      better agree with nparam, or we are in trouble.
c
      nparam0 = kinds*(29+68*kinds)
      if(nparam.ne.nparam0) then
         if(master) then
            write(0,*) 'having ',kinds,' types of atoms requires ',
     $           nparam0,' parameters'
            write(0,*) 'you have ',nparam
            write(0,*) 'static:  nparm <> nparam0(kinds)'
         end if
         go to 10000
      end if
c
c
c     At this point we have all of the preliminary information needed
c      by the program, including the kinds of atoms, the number of
c      parameters, etc.  Pass this information to all processors:
c
c     Now we need to communicate our results from processor 0 to the
c      rest of the world.
c
c     from input.f we have
c
c     realov, which tells us the type of overlap matrix elements
c      (old or new)
c
      call MPI_BCAST(realov,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
c
c     param, the parameter array itself.  Note that we'll transmit all
c      of the parameters
c
      call MPI_BCAST(param,npard*mspin,MPI_DOUBLE_PRECISION,0,
     $     MPI_COMM_WORLD,ierr)
c
ctemp
c
c     Make sure our parameters know what's what:
c
c$$$      write(*,'(''Processor '',i2,'' finds realov = '',l1)')
c$$$     $     MPIRANK,realov
c$$$      write(*,'(''Processor '',i2,'' finds realov = '',a1)')
c$$$     $     MPIRANK,realov
c$$$      write(*,'(''Processor '',i2,'' finds p = '',1pe20.12)')
c$$$     $     MPIRANK,param(27,1)
c$$$      write(*,*) MPIRANK," npkind = ",npkind
c$$$      write(*,*) MPIRANK," jsover = ", jsover
c$$$      write(*,*) MPIRANK," kbas = ",(kbas(i),i=1,kinds)
cend temp
c
c     from input0 we have:
c
c
cend MPI
c
      jstruc = 0
c
c     this should really be a "do while no eof" loop,
c      hence the indentation
c
 1000 continue
         jstruc = jstruc + 1
ctemp
c$$$         write(*,*) 'Processor ',MPIRANK,' on jstruc = ',jstruc
cend temp
c
c        here's where we'll hit the eof:
c
c
c        Read in from processor 0, as before
c
         if(master) then
c
c           Read "title" for SKOUT and "label" for SKENG.  Neither of
c            these strings need leave the master processor.  Keep
c            checking for an eof or other SKIN file error, as this
c            is where it is most likely to happen.  (E.g., blank lines
c            entered at the end of the SKIN file.)
c
            read(10,'(a70)',iostat=ioerr) title
            if(ioerr.eq.0) read(10,'(a20)',iostat=ioerr) label
c
c
c           Read zcharge is the number of electrons ABOVE OR BELOW THE
c            CHARGE CALCULATED FROM THE PARAMETER FILE CHARGES.
c
cMPI note
c           With the parallelization of fermi.f, the other processors
c            won't need zcharge, but they will need zfill
cend MPI note
c
            if(ioerr.eq.0) read(10,*,iostat=ioerr) zcharge
c
c           The title and zfile lines are the most common places we
c            exit the program.  As such, we need to tell other
c            processors to go, too.
c
         end if
c
         call MPI_BCAST(ioerr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c
         if(ioerr.ne.0) go to 10000
c
c        Read and write from processor 0.
c
         if(master) then
            write(15,'(''-------------------------------''//
     $           ''structure '',i5//a70/)') jstruc,title
            write(15,*) zcharge, ' electrons above nominal charge'
c
c           Read in and manipulate structures.  If input1 has a read
c            error, then inpstat will return with a non-zero value.
c
c           input1 will determine the actual valence charge zfill,
c            and the volume of the current unit cell (thisvol):
c
            call input1(zcharge,zfill,thisvol,inpstat,flipos,
     $           lnoncol,ltorque,lspdir,spinat,dspth,dspphi)
c
            if(inpstat.ne.0) write(*,*) 'Error ',inpstat,
     $           ' in input1 detected.  The program will exit'
c
c           Check the compatability of lnoncol and jspins:
c
            if(lnoncol.and.(jspins.eq.1)) then
               write( *,*) 'JSPINS = 1, so no magnetization'
               write( *,*) 'Non-collinear spins are turned off'
c
               write(15,*) 'JSPINS = 1, so no magnetization'
               write(15,*) 'Non-collinear spins are turned off'
c
               lnoncol = .false.
c
c              Note that ltorque must be false if lnoncol is
c
               if(ltorque) then
                  ltorque = .false.
                  write( *,*) 'No spin derivative calculations'
                  write(15,*) 'No spin derivative calculations'
               end if
c
c              Ditto for lspdir
c
               if(lspdir) then
                  ltorque = .false.
                  write( *,*) 'No output spin vector calculations'
                  write(15,*) 'No output spin vector calculations'
               end if
            end if
c
c           Open up file 41 for writing spin derivatives.  Note that
c            only the master processor will write to this file.
c
            if(ltorque) open(unit=41,file='SKDERIV',status='unknown')
c
c           Similarly, open up file 42 if we want the output
c            spin vectors to be calculated
c
            if(lspdir) open(unit=42,file='SPNDIR',status='unknown')
c
c           This is the end of the input for the current structure
c
         end if
c
c        If inpstat is non-zero we need to quit
c
         call MPI_BCAST(inpstat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c
         if(inpstat.ne.0) go to 10000
c
c        Now communicate the necessary information generated by input1
c
         call MPI_BCAST(nv,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(natoms,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(ksk,matom,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(kkind,matom,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(nkdarry,6,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(plv,9,MPI_DOUBLE_PRECISION,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_BCAST(posn,3*matom,MPI_DOUBLE_PRECISION,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_BCAST(posq,3*matom,MPI_DOUBLE_PRECISION,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_BCAST(dnn,mint,MPI_DOUBLE_PRECISION,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_BCAST(dnnsq,mint,MPI_DOUBLE_PRECISION,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_BCAST(screenl,mint,MPI_DOUBLE_PRECISION,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_BCAST(scrninv,mint,MPI_DOUBLE_PRECISION,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_BCAST(flipos,matom,MPI_LOGICAL,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_BCAST(zfill,1,MPI_DOUBLE_PRECISION,0,
     $      MPI_COMM_WORLD,ierr)
         call MPI_BCAST(lnoncol,matom,MPI_LOGICAL,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_BCAST(ltorque,matom,MPI_LOGICAL,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_BCAST(spinat,3*matom,MPI_DOUBLE_PRECISION,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_BCAST(lspdir,matom,MPI_LOGICAL,0,
     $        MPI_COMM_WORLD,ierr)
c
c        Export the derivatives of spinat if needed:
c
         if(ltorque) then
            call MPI_BCAST(dspth,3*matom,MPI_DOUBLE_PRECISION,0,
     $           MPI_COMM_WORLD,ierr)
            call MPI_BCAST(dspphi,3*matom,MPI_DOUBLE_PRECISION,0,
     $           MPI_COMM_WORLD,ierr)
         end if
c
c        thisvol is needed by setvol to calculate the pressure:
c
         call MPI_BCAST(thisvol,1,MPI_DOUBLE_PRECISION,0,
     $      MPI_COMM_WORLD,ierr)
c
cMPI note
c        Note that these calculations might as well be done on
c         each processor, otherwise they are just sitting idle:
c
cend MPI note
c
c        The following code used to be in the "multatoms" subroutine:
c
         if(master) then
            write(15,'('' FINAL ATOMIC COORDINATES:'')')
         end if
         do iat=1,natoms
            do i=1,3
               tempi=zero
               do j=1,3
                  tempi=tempi+plv(j,i)*posn(iat,j)
               end do
               temp(i) = tempi
            end do
c
            if(master) write(15,'(3f12.7,5x,3f12.7)')
     $           (posn(iat,i),i=1,3),(temp(i),i=1,3)
c
            do i=1,3
               posn(iat,i)=temp(i)
            end do
         end do
c
c        That was the end of multatoms
c
cMPI note
c        A parallel version of search2 would probably not penalize
c         the program, and would be useful in large systems.  For now,
c         however, search2 is called by all processors.
c
cend MPI note
c
         call search2(npfound,dismin,ierr)
c
c        Output certain diagnostic information on the master processor:
c
         if (master) then
            write(15,*) 'Minimum distance this configuration = ',dismin
            write(15,'(''Number of pairs = '',i6)') npfound
         end if
c
c        In this case the error message has already been printed,
c         so we can just exit if ierr > 0:
c
         if(ierr.gt.1) go to 10000
c
ctemp
c$$$         write(*,*) MPIRANK,' found ',npfound,' pairs'
cend temp
c
         if(master) then
c
c           qlmt file header (if used).  The fourth field should be the
c            number of atoms in the unit cell, or 0 if we are using the
c            lqlmt0 option
c
            if(mode.eq.5) then
               write(45) 1,1,0,natoms
            else if((mode.eq.6).or.(mode.eq.10)) then
               write(45,45) 1,1,0,natoms
 45            format(4i5,5x,'SPINS, WIND., SC. WIND.')
            else if(mode.eq.9) then
               write(45,45) 1,1,0,0
            end if
         end if
c
c        generate k-points for each structure.  Note that kptin
c           will have a serial section where it reads in the k-points
c
cMPI note
c        I thought about making kptin work for all structures, but
c         it is easier to read the information in on one processor
c         and transfer the data.
cend MPI note
         if(master) then
            call kptin
ctemp
c$$$            write(*,*) 'Back from kptin'
cend temp
         end if
c
c        Now the k-points are calculated, so spill them back out again.
c         Note that we will have less overhead if we get the
c         (x(i),y(i),z(i)) varables changed to k(3,i) type variables
c
c        Don't forget that we need to transfer the number of k-points, too.
c
         call MPI_BCAST(npts,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         if(ierr.ne.0) go to 10000
         call MPI_BCAST(x,nkd,MPI_DOUBLE_PRECISION,0,
     $        MPI_COMM_WORLD,ierr)
         if(ierr.ne.0) go to 10000
         call MPI_BCAST(y,nkd,MPI_DOUBLE_PRECISION,0,
     $        MPI_COMM_WORLD,ierr)
         if(ierr.ne.0) go to 10000
         call MPI_BCAST(z,nkd,MPI_DOUBLE_PRECISION,0,
     $        MPI_COMM_WORLD,ierr)
         if(ierr.ne.0) go to 10000
c
c        K-point weights:
c
         call MPI_BCAST(wghtk,nkd*nwdd,MPI_DOUBLE_PRECISION,0,
     $        MPI_COMM_WORLD,ierr)
ctemp
c$$$         write(*,*) MPIRANK,' has ',npts,' k-points'
c$$$         write(*,*) 'k5(',MPIRANK,') = ',x(5),y(5),z(5)
cend temp
c
c        Final output to processor 0:
c
         if(master) then
c
c           Add a little bit of information to the qlmt file so that
c            we can later use it to perform sums over occupied states:
c
            if(mode.eq.5) then
               write(45) nkv,natoms,zfill,tkb
            else if((mode.eq.6).or.(mode.eq.9).or.(mode.eq.10)) then
               write(45,55) nkv,natoms,zfill,tkb
 55            format(i5,' K-POINTS',i6,2f10.5)
            end if
         end if
c
c        bande is rewritten in parallel.  We now pass the volume
c         this way to get it into setvol:
c      
ctemp
c$$$         write(*,*) 'Calling bande'
cend temp
c
c        Note that the quantities in earg1 and earg2 depend
c         on the mode.
c        Also note that the call to bande has been modified:
c         In 1.20, to add the flipos array, and
c         in 1.30, to add the non-collinear arrays.
c
         call bande(mode,jspins,flipos,lnoncol,spinat,
     $        ltorque,dspth,dspphi,lspdir,
     $        thisvol,zfill,tkb,eigcut,wghtk(1,1),efermi,
     $        earg1,earg2,smag,lp_par,lp_parv,lp_klst,kptpr,
     $        label,iberr)
ctemp
c$$$         write(*,*) 'Back from bande, iberr = ',iberr
cend temp
         if(iberr.ne.0) then
            write(0,*) 'bande returned with error = ',iberr
            go to 10000
         end if
c
c        Final output to processor 0:
c
         if(master) then
c
c           Add the magnetic moment smag to SKENG if this is a
c            collinear magnetic system (jspins = 2 and lnoncol=false)
c
            if((jspins.eq.1).or.lnoncol) then
               if((mode.lt.3).or.(mode.eq.4)) then
                  write(31,'(a20,f12.6,3f15.9)')
     $                 label,thisvol,efermi,earg1,earg2
               else
                  write(31,'(a20,f12.6,2f15.9)')
     $                 label,thisvol,efermi,earg1
               end if
            else
               if((mode.lt.3).or.(mode.eq.4)) then
                  write(31,'(a20,f12.6,3f15.9,f10.6)')
     $                 label,thisvol,efermi,earg1,earg2,smag
               else
                  write(31,'(a20,f12.6,2f15.9,f10.6)')
     $                 label,thisvol,efermi,earg1,smag
               end if
            end if
         end if
ctemp
c$$$         write(*,*) 'Processor ',MPIRANK,' done with jstruc = ',jstruc
cend temp
      go to 1000
c
c     Here's where we finish up:
c
10000 continue
c
c     So long as mode=7 is inactive, we can do it this way:
c
      if(mode.gt.4) close(45)
c
cmpi
c
ctemp
c$$$      write(*,*) 'Processor ',MPIRANK,' ready to exit'
cend temp
      CALL MPI_FINALIZE (ierr)
c
cend MPI
c
      end
