      subroutine tweaks (lp_par,lp_parv,lp_klst,kptpr)
c
c     Created 7 August 2000
c     Static code version 1.21
c     By Michael Mehl (mehl@dave.nrl.navy.mil)
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
c     Examines the file TWEAKS, unit number IO, for keywords, as
c      outlined below.  If the keyword exists, appropriate action is
c      taken.
c
c=======================================================================
c
c     REVISION HISTORY
c
c     30 April 2002:
c
c     Added the tweak "PRINT_DIMENSIONS" which prints all of the
c      dimensions parameters it finds in the include files.  Note
c      that this does not require any further communication to the
c      outside world, so there is no "lp_dim" variable set.
c
c     28 August 2001:
c
c     Added the tweak "PRINT_PROGRESS", which prints keeps track of where
c      we are in the program by printing the current k-point to
c      standard output.  This adds the logical variable "lp_klst" to
c      the argument list, and the keyword variable printprog.
c
c=======================================================================
c
      implicit none
c
c=======================================================================
c
c     Scalar integers
c
c     Integer keywords.  See program for full description:
c
      integer kptpr
c
c=======================================================================
c
c     Logical variables:
c
c     Logical keywords.  See program for full description:
c
      logical lp_par,lp_parv,lp_klst
c
c=======================================================================
c
c     The keywords.  See below for description.  Note that
c
c     The input keyword is:
c
      character*21 keyword
c
c     Other keywords:
c
      character*21 printdim, printpar, printparv, printmat, printprog
c
c=======================================================================
c
c     pound ('#') is the comment symbol.  Lines in TWEAKS beginning with
c      a '#' will not be acted upon.  This allows one to leave all
c      possible TWEAKS options in place but not act on them
c
      character*1 pound
c
c=======================================================================
c
c     Data for keywords and pound:
c
      data printpar  /'PRINT_PARAMETERS     '/
      data printparv /'PRINT_DERIVATIVE_PAR '/
      data printmat  /'PRINT_MATRIX         '/
      data printprog /'PRINT_PROGRESS       '/
      data printdim  /'PRINT_DIMENSIONS     '/
c
      data pound /'#'/
c=======================================================================
c
c     There are several sets of keywords.
c
c     At the moment all of the keywords pass through the subroutine
c      bande.f:
c
c-----------------------------------------------------------------------
c
c     KEYWORD: PRINT_DIMENSIONS
c     VARIABLE: None
c     ALLOWED VALUES: N/A
c
c     If found, immediately prints all of the dimensions found in
c      the include files P1, P2, and P3 onto SKOUT, then continues.
c      Useful for determining what kind of static executable you
c      actually have.
c
c-----------------------------------------------------------------------
c
c     KEYWORD: PRINT_DERIVATIVE_PAR
c     VARIABLE: lp_parv
c     ALLOWED VALUES: .true.,.false.
c
c     lp_parv is a similar flag for setparv.f, the routine which
c      calculates derivative parameters.
c
      lp_parv = .false.
c-----------------------------------------------------------------------
c
c     KEYWORD: PRINT_MATRIX
c     VARIABLE: kptpr
c     ALLOWED VALUES: integer
c
c     The subroutine bande.f tells setham.f to print the H and S
c      matrices for k-points i.le.kptpr.  The default value of
c      kptpr is, of course, zero:
c
      kptpr = 0
c-----------------------------------------------------------------------
c
c     KEYWORD: PRINT_PARAMETERS
c     VARIABLE: lp_par
c     ALLOWED VALUES: .true.,.false.
c
c     if lp_par is true, then diagnostic information is printed
c      in setpar.f.  This is mainly the values of the on-site,
c      hopping, and overlap matrix elements.  The default value
c      is false.
c
      lp_par = .false.
c-----------------------------------------------------------------------
c
c     KEYWORD: PRINT_PROGRESS
c     VARIABLE: lp_klst
c     ALLOWED VALUES: .true.,.false.
c
c     if lp_klst is true, then the current k-point will be printed
c      to standard output.  This helps track execution of the program
c      if a large number of k-points is involved.  Default is off:
c
      lp_klst = .false.
c
c=======================================================================
c
c     Does Tweaks exist?  Try opening the it:
c
      open(unit=17,file='TWEAKS',status='old',err=10000)
c
c     If we get here the file is open, so report:
c
      write(*,*) 'TWEAKS file exists'
      write(15,*) 'TWEAKS file exists'
c
c     Read a line in TWEAKS, try to match it with a keyword.
c
c     Logical keywords are toggled from the default to .not.default
c
c     Note that if an integer or real-valued keyword exists,
c      its value is set ON THE LINE AFTER the keyword.
c
 1000 read(17,'(a21)',end=9000) keyword
c
c     Check for comment.  If not, try to act:
c
      if(keyword(1:1).eq.pound) then
c$$$         write(*,'(''TWEAKS comment: '',a21)') keyword
         continue
      else if(keyword.eq.printpar) then
         lp_par = .true.
         write(*,*)
     $        'On-site, hopping, and overlap SK parameters printed'
      else if(keyword.eq.printparv) then
         lp_parv = .true.
         write(*,*)
     $        'On-site, hopping, and overlap derivatives printed'
      else if(keyword.eq.printmat) then
         read(17,*) kptpr
         write(*,*)
     $        'H and S printed for first',kptpr,' points'
      else if(keyword.eq.printprog) then
         lp_klst = .true.
         write(*,*)
     $        'k-point progression will be tracked in standard output'
      else if(keyword.eq.printdim) then
c
c        Print out all of the parameters defined in the P1, P2 and P3
c         include files.  So as not to take up too much space here,
c         we'll do this as a subroutine call.
c
         call prtdim
c
      else
         write(*,'(''Keyword'',a21,'' did not match'')') keyword
      end if
c
c     Read a new line:
c
      go to 1000
c
c     Done, so close the file an exit:
c
 9000 close(17)
c
c     Goodbye.  We go directly here if TWEAKS does not exist:
c
10000 return
      end
c_______________________________________________________________________
c
      subroutine prtdim
      implicit none
c
c     Prints all of the dimension parameters it finds in the files:
c
      include 'P1'
      include 'P2'
      include 'P3'
c
c     onto SKOUT
c
c     Note that each variable is typed in the appropriate Pi file
c
      write(15,'(//
     $     ''*********************************************************''
     $     )')
      write(15,'(//''Parameters found in P1:''//)')
      write(15,'(''mkind            = '',i10)') mkind
      write(15,'(''matom            = '',i10)') matom
      write(15,'(''mbas             = '',i10)') mbas
      write(15,'(''mppair           = '',i10)') mppair
      write(15,'(''mspin            = '',i10)') mspin
      write(15,'(''mcol             = '',i10)') mcol
      write(15,'(''npard            = '',i10)') npard
      write(15,'(''mint             = '',i10)') mint
      write(15,'(''mh               = '',i10)') mh
      write(15,'(''mptype           = '',i10)') mptype
      write(15,'(''mpair            = '',i10)') mpair
      write(15,'(//''Parameters found in P2:''//)')
      write(15,'(''mpress           = '',i10)') mpress
      write(15,'(''mtorque          = '',i10)') mtorque
      write(15,'(''msdir            = '',i10)') msdir
      write(15,'(''meigen           = '',i10)') meigen
      write(15,'(''mhe              = '',i10)') mhe
      write(15,'(''mhp              = '',i10)') mhp
      write(15,'(//''Parameters found in P3:''//)')
      write(15,'(''nkd              = '',i10)') nkd
      write(15,'(''nopd             = '',i10)') nopd
      write(15,'(''minproc          = '',i10)') minproc
      write(15,'(''nkdproc          = '',i10)') nkdproc
      write(15,'(''nwdd             = '',i10)') nwdd
      write(15,'(//
     $     ''*********************************************************''
     $     //)')
      return
      end
