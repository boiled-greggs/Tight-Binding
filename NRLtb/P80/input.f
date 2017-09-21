      subroutine input (nparam,jspins)
c
c     Reads in the input file from the file named in xparfl, and
c      returns the number of parameters in the file
c
c=======================================================================
c
c     REVISION HISTORY:
c
c-----------------------------------------------------------------------
c
c     Changed to read in a file name which contains the parameters, then
c      read the parameters from that file
c
c-----------------------------------------------------------------------
c
c     Note that the results of this code do not depend on the sign
c      of param(1,1).  Thus, if param(1,1) < 0, we assume that this
c      is a flag telling us that we should use the constrained version
c      of the onsite like-atom overlap parameters.  This is signaled
c      by the flag "realov".
c                                            -- mjm  26 Sept 1997
c
c-----------------------------------------------------------------------
c
c     Version 1.05:
c
c     Now certain information which used to be read in from the SKIN
c      file is now contained in the parameter file, and so is read
c      here.  This includes the cutoff distances dnn and the cutoff
c      screening parameter screenl, which now can be different for
c      the different kinds of interactions.  The number of atom types
c      in the parameter file, kinds, is now also read in here.
c
c                                            -- mjm   7 April 1998
c
c-----------------------------------------------------------------------
c
c     The "NN" at the head of the parameter file is given new meaning.
c      See the start of the executable.
c                                            -- mjm  20 April 1998
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
c     Added a special mode "90000", for the H atom.  This invokes
c      an extended polynomial for the ss_sigma H and S parameters.
c     Changed the common block which was named /overlap/ to /parcom/.
c     Eliminated the useless variable "iwrite" from /lia/ and the code.
c
c                                                mjm --  4 Aug   1998
c-----------------------------------------------------------------------
c
c     Took "jspins" out of the common block /lia/ and moved it into
c      the calling list.  This facilitates communication with input1.f
c      for use in spin-flipped and non-colinear magnetization.
c                                                mjm -- 21 Dec   1999
c-----------------------------------------------------------------------
c  
c     Let A-A, B-B, C-C, etc. interactions have Harrison's canonical
c     sign.  i.e.:
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
      implicit real*8 (a-h,o-z)
      include 'P1'
      logical jsover
      common /params/ param(npard,mspin)
      common /lia/ npts,jsover
c
c     Basic parameter file information:
c
      common /parinfo/ valence(mkind),kinds,kbas(mkind)
c
c     Cutoff and screening arrays:
c
      common /cutcom/ dnnsq(mint),dnn(mint),
     $     screenl(mint),scrninv(mint)
c
      logical realov,harrison
      common /parcom/ nltype,realov,harrison
c
c     ovtype and mgtype control the overlap and spin characteristics
c      the parametrization file is using:
c
      character*1 ovtype,magtype
c
c     Parameter title:
c
      character*79 partitle
c
c     Parameter file path name:
c
      character*79 xparfl
c
c=======================================================================
c
      parameter (zero = 0d0)
      parameter (one  = 1d0)
c
c=======================================================================
c
c     Certain characters needed by the program:
c
      character*1 muc,mlc,nuc,nlc,ouc,olc
c
      parameter (muc = "M")
      parameter (mlc = "m")
      parameter (nuc = "N")
      parameter (nlc = "n")
      parameter (ouc = "O")
      parameter (olc = "o")
c
c=======================================================================
c
      read(10,'(a79)') xparfl
      open(unit=32,file=xparfl,status='old',err=2000)
c
c     Parameter form control.  The first line of the parameter file
c      contains several characters and numbers, revealing the type
c      of parametrization file that is contained here.
c
c     The first of these tells us what kind of basis functions we
c      are using:
c
c     "N" indicates non-orthogonal basis functions
c     "O" indicates orthogonal basis functions.  Note that in this
c         case S matrix parameters will be read in, but ignored.
c     Default is "N" with a warning message.
c
c     The second letter tells us how many spins are involved.
c
c     "N" indicates the system is non-magnetic, i.e., jspin = 1
c     "M" indicates the system is magnetic.  In this case mspin
c         in the P1 file must be dimensioned to 2 or the program
c         will bomb, and the file must contain all parameters in
c         the proper order.
c     Default is "N" with a warning message.
c
c     The next field (I5 format), contains the a revision number
c      indicating the type of parameters being used.  This may
c      or may not change the number of parameters read, depending
c      on the revision number of the file.
c
c     "    0" or "     " indicates old-style overlap parameters
c       (a + b r + c r^2)exp(-d^2 r)
c     "    1" indicates a newer style:
c       (delta_(l,l') + a r + b r^2 + c r^3) exp(-d^2 r)
c     "    2" indicates the "Harrison sign" version of the parameters
c       as noted above, with old-style overlap parameters
c     "    3" indicates the "Harrison sign" version of the parameters
c       as noted above, with new-style overlap parameters
c
c     Both of these options use the same number of parameters.
c
c     At the present time, anything else defaults to "0" with a warning.
c
      read(32,'(2a1,i5)') ovtype,magtype,nltype
ctemp
c$$$      write(*,'(''Line 1:'',2a1,i5)') ovtype,magtype,nltype
cend temp
c
c     Overlap specification:
c
      if ((ovtype.eq.nuc).or.(ovtype.eq.nlc)) then
c
c        Non-orthogonal calculation (with overlap matrix)
c
         jsover = .true.
      else if ((ovtype.eq.ouc).or.(ovtype.eq.olc)) then
c
c        Orthogonal calculation (no overlap matrix)
c
         jsover = .false.
      else
c
c        Default is non-orthogonal:
c
         write( *,*) 'INVALID BASIS FUNCTION TYPE -- ASSUMED "N"'
         write(15,*) 'INVALID BASIS FUNCTION TYPE -- ASSUMED "N"'
         jsover = .true.
      end if
      if(jsover) then
         write( *,*) 'Non-Orthogonal Basis Functions'
         write(15,*) 'Non-Orthogonal Basis Functions'
      else
         write( *,*) 'Orthogonal Basis Functions'
         write(15,*) 'Orthogonal Basis Functions'
      end if
c
c     Magnetization specification:
c
      if ((magtype.eq.nuc).or.(magtype.eq.nlc)) then
         jspins = 1
      else if ((magtype.eq.muc).or.(magtype.eq.mlc)) then
         jspins = 2
         if(mspin.lt.2) then
            write(*,*) 'MSPIN must be > 1 for magnetic calculation'
            stop 'MSPIN too small'
         end if
      else
         write( *,*) 'INVALID MAGNETIZATION TYPE -- ASSUMED "N"'
         write(15,*) 'INVALID MAGNETIZATION TYPE -- ASSUMED "N"'
         jspins = 1
      end if
      if(jspins.eq.2) then
         write( *,*) 'Spin-polarized calculation'
         write(15,*) 'Spin-polarized calculation'
      else
         write( *,*) 'Spin-restricted calculation'
         write(15,*) 'Spin-restricted calculation'
      end if
c
c     Parameter form.  nltype = 0 is default
c
      if((nltype.ne.90000).and.((nltype.gt.3).or.(nltype.lt.0))) then
         write(*,*) 'PARAMETERS ASSUMED TO BE USING DEFAULT STYLE'
         nltype = 0
      end if
c
c     nltype = (1,3) if we use the new overlap parameter style:
c
      realov = (nltype.eq.1).or.(nltype.eq.3)
      if(realov) then
         write( *,*) 'Using S(l,lp,R) -> delta(l,lp) for like atoms'
         write(15,*) 'Using S(l,lp,R) -> delta(l,lp) for like atoms'
      end if
c
c     nltype = (2,3) if we use the "Harrison" sign convention
c      (see the header file)
c
      harrison = (nltype.eq.2).or.(nltype.eq.3)
      if(harrison) then
         write( *,*) 'Using Harrison sign convention for like atoms'
         write(15,*) 'Using Harrison sign convention for like atoms'
      end if
c
c     This is just a title for the parameters:
c
      read(32,'(a79)') partitle
      write(*,'(a79)') partitle
c
c     How many atom types in this file?
c
      read(32,*) kinds
      if(kinds.gt.mkind) stop 'input.f:  kinds > mkind'
c
c        How many parameters should we read?  This depends on
c         the version number.  For now, however, it is the same
c         for both of the parameter styles defined:
c
      npards = kinds*(29+68*kinds)
c
c     This is a check on the logic in the P1 file.  If this is
c      wrong, the npard (in P1) or npards has changed:
c
      if(npards.gt.npard) then
         write(*,*) 'input.f:'
         write(*,*) 'Found kinds = ',kinds,' npards = ',npards
         write(*,*) 'while in P1'
         write(*,*) 'mkind = ',mkind,' npard = ',npard
         stop 'input.f:  npards > npard (P1)'
      end if
c
c     There are kinds*(kinds+1)/2 interaction types.  For each
c      one, read in the cutoff distance and screening parameter,
c      then manipulate these parameters as necessary.
c
      do i = 1,kinds*(kinds+1)/2
         read(32,*) dnn(i),screenl(i)
         dnnsq(i) = dnn(i)*dnn(i)
         scrninv(i) = one/screenl(i)
      end do
      write( *,'(/''Cutoff parameters: ''/(i5,f10.5,f10.5))')
     $     (i, dnn(i), screenl(i), i = 1,kinds*(kinds+1)/2)
      write(15,'(/''Cutoff parameters: ''/(i5,f10.5,f10.5))')
     $     (i, dnn(i), screenl(i), i = 1,kinds*(kinds+1)/2)
c
c     Now we need some atom-specific information:
c
      do i = 1,kinds
c
c        First is the number of orbitals we'll use for each atom kind:
c
         read(32,*) kbas(i)
c
c        Next is the atomic weight.  This is useful for MD, but not
c         for us.
c
         read(32,*) atmwght
         write(*,*) 'Atomic weight of species ',i,' is ',atmwght
c
c        Finally we read in the formal valence of each atom.  The
c         parameter file lists s, p, and d, but we only need the
c         sum:
c
         read(32,*) socc,pocc,docc
         valence(i) = socc + pocc + docc
         write(*,*) 'Atom ',i,' has ',valence(i),' electrons'
      end do
c
c     Read off the parameters for each spin:
c
      do jspin=1,jspins
c
         do i=1,npards
c
c           Note that there may actually be fewer parameters
c            in the calculation than we've assumed.
c
            read(32,*,iostat=ioerr) param(i,jspin)
            if(ioerr.ne.0) then
c
c              There are i-1 parameters in the file:
c
               nparam = i-1
c
               if(nparam.lt.npards)
     $              write(*,*) 'WARNING:  zeroing parameters',
     $              nparam+1,' to ',npards
c
c              Zero the rest of the parameters (just in case)
c
               do j = i,npards
                  param(j,jspin) = zero
               end do
               go to 500
            end if
         end do
c
c        Default value of nparam is npards:
c
         nparam = npards
c
 500     write(15,'(5x,i5,e20.10)') (i,param(i,jspin),i=1,nparam)
      end do
      return
 2000 write(0,2005) xparfl
 2005 format(/'Cannot find file'/a79)
      stop 'input:  Exit 2'
      end
