      real*8 function fermi(ntbd,nkptd,nkloc,
     $     istart,ifin,nv,jspins,el,weight,eigen,tkb,zfill)
c
c     Determines the Fermi level at temperature tkb by direct
c      evaluation of the sum
c
c     zfill(mu) = Sum[
c                  Sum[weight(k)*
c                   Sum[1/(1+Exp(tkb*(eigen(i,k,j)-mu))) , {i=1,nv}]
c                                                        , {k=1,nkpt}]
c                                                        , {j=1,jspins}]
c
c     Where:
c
c     eigen(i,k,j) is the ith eigenvalue at the kth k-point for the
c      jth spin
c     el is the number of electrons in a given band (new in 1.30)
c     weight(k) is the weight for the kth k-point
c     jspins = 1 for para, 2 for ferromagnetic calculations
c     nkpt is the number of k-points.  In this multi-processor version,
c      a given node sums from istart to ifin.
c     nv is the number of eigenvalues at each k-point
c      note:  this could easily be generalized to nv(k)
c     tkb is the temperature
c     zfill is the number of electrons
c
c     Note:  We assume that the eigenvalues at each k-point are ordered,
c      lowest to highest.  Most eigenvalue finding routines will do this
c      for you.
c
c                                           -- M.J. Mehl, Sept., 1994
c
c-----------------------------------------------------------------------
c
c     Version 1.22:
c
c     This version is written in parallel.  It assumes that the
c      eigenvalues stay on the nodes on which they were generated.
c      In particular, the k-points istart to ifin are on the current
c      processor.
c-----------------------------------------------------------------------
c
c     Version 1.30:
c
c     In order to distinguish between paramagnetic, spin-polarized,
c      and non-collinear magnetism, the number of electrons/band is
c      now input in as el, in the function call.  This replaces the
c      calculation of used in the previous versions.
c                                                 -- mjm  11 August 2000
c-----------------------------------------------------------------------
c
c     Version 1.31:
c
c     Merger of 1.21->1.22 and 1.21->1.30 changes.
c                                                -- mjm  23 October 2000
c-----------------------------------------------------------------------
c
c     Version 1.33:
c
c     Takes part in the reduction of eigenvalue storage to something
c      approaching the number appropriate for the current processor.
c      See notes in REVISION_HISTORY and bande.f for more details.
c     Note that we have added nkloc to the call.  This is the k-point
c      dimension for the eigenvalues, etc., on the local processor
c
c                                                -- mjm   15 August 2002
c=======================================================================
c
      implicit none
c
c=======================================================================
c
      include 'mpif.h'
      logical master,multiprc
c
      integer MPISIZE,MPIRANK
      common /MPIcomm/MPISIZE,MPIRANK,master,multiprc
c
c=======================================================================
c
c     These set the dimensions of weight and eigen:
c
      integer ntbd,nkptd,nkloc
c
c=======================================================================
c
      integer nv,jspins,istart,ifin
c
c     Version 1.33: note that eigen is dimensioned for the maximum
c      number of k-points on the current processor.
c
      real*8 weight(nkptd),eigen(ntbd,nkloc,*),tkb,zfill
c
c=======================================================================
c
c     beta is the inverse temperature
c     mu is the trial value of the Fermi level
c     btmu = -beta*mu
c     el is the number of electrons/level (in function call for
c      Version >= 1.30)
c     mulow and muhigh are upper and lower bounds on mu
c     nsum is the Fermi sum, dnsum its derivative wrt mu
c     nksum and dnksum are the values at each k-point
c     ex = exp(beta*(e-mu))
c     f  = 1/(1+ex)
c
      real*8 beta,mu,btmu,el,mulow,muhigh,nsum,dnsum,nksum,dnksum,ex,f
c
c     variables ending in "g" are global variants of the above
c
      real*8 nsumg,dnsumg
c
c     mulowl and muhighl are the local-node values of mulow and muhigh.
c
      real*8 mulowl,muhighl
c
c     mum20 and mup20 bound the range where the fermi function is
c      neither 1 nor 0
c     mu1 is a temporary value of mu used by the Newton's method
c     stlow and sthigh are estimates of the proper Fermi level
c
c     stlowl and sthighl are local processor values of stlow and
c      sthigh
c
      real*8 mum20,mup20,mu1,stlow,sthigh,stlowl,sthighl
c
c     nfill is the next integer larger than zfill
c
      integer nfill
c
c     do loop and error variables
c
      integer i,j,k,kloc,ierr
c
c     linst is true if the system may be an insulator,
c      otherwise it must be a metal.
c
      logical linst
c
c=======================================================================
c
      real*8 zero,one,two,twenty
c
      parameter (zero = 0d0)
      parameter (one = 1d0)
      parameter (two = 2d0)
      parameter (twenty = 2d1)
c
      real*8 half
c
      parameter (half = one/two)
c
      real*8 nearz,tol
c
      parameter (nearz = 1d-5)
      parameter (tol = 1d-8)
c
c=======================================================================
ctemp
c$$$      write(*,*) 'Welcome to fermi.f, processor ',MPIRANK,istart,ifin
cend temp
c
c     Can a possible solution exist?
c
      if(zfill.le.zero) stop 'FERMI:  No solution for zfill <= 0'
c
c     Set the inverse temperature.  If TKB is small, we'll use a default
c      value so that we're not dividing by zero
c
      if(tkb.lt.nearz) then
         beta = one/nearz
      else
         beta = one/tkb
      end if
c
c     nfill is the number of levels that must be filled if the
c      crystal is an insulator.
c
      nfill = int(zfill/el+nearz)
c
c     If nfill = zfill/el, then this crystal may be an insulator
c
      linst = abs(zfill/el-dble(nfill)).lt.tol
ctemp
c$$$      if(linst) write(*,*) 'Possible insulator with ',nfill,' levels'
cend temp
c
c     Scan through the eigenvalues to find realistic values for
c      the bounds on the Fermi level, and for the starting
c      fermi level
c
c     Version 1.33: since eigenvalues are stored locally, we replace
c      istart by 1 here:
c
c$$$      muhighl = eigen(nv,istart,1)
c$$$      mulowl  = eigen(1,istart,1)
      muhighl = eigen(nv,1,1)
      mulowl  = eigen(1,1,1)
c
c$$$      stlowl = eigen(nfill,istart,1)
c$$$      sthighl = eigen(nfill+1,istart,1)
      stlowl = eigen(nfill,1,1)
      sthighl = eigen(nfill+1,1,1)
c
      do j = 1,jspins
c
c        Version 1.33: we need to keep track of the local k-point:
c
         kloc = 0
c
         do k = istart,ifin
c
            kloc = kloc+1
c
c           Note that in metals stlow may well be above sthigh
c
c           Version 1.33: remember eigenvalues are stored locally:
c
            stlowl = max(stlowl,eigen(nfill,kloc,j))
            sthighl = min(sthighl,eigen(nfill+1,kloc,j))
c
            do i = 1,nv
               muhighl = max(muhighl,eigen(i,kloc,j))
               mulowl  = min(mulowl ,eigen(i,kloc,j))
            end do
         end do
      end do
c
c     Transfer all of the stlow, sthigh, mulow, and muhigh
c      back to the master processor, figure out the starting
c      Fermi level.
c
c     Note that MPI_ALLREDUCE will do the global properties for us:
c
c     We could define a fake MPI_ALLREDUCE in mpifake.f which would
c      copy stlowl to stlow, etc., but that would have to be
c      written for double precision variables only.  Since we
c      don't, a priori, know how many versions of MPI_ALLREDUCE
c      we will need, well issue an exemption for single
c      processor runs:
c
      if(multiprc) then
c
c        Note that we want the maximum of stlowl, the minimum
c         of sthighl, the minimum of mulowl, and the maximum
c         of muhighl:
c
         call MPI_ALLREDUCE(stlowl,stlow,1,MPI_DOUBLE_PRECISION,
     $        MPI_MAX,MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(sthighl,sthigh,1,MPI_DOUBLE_PRECISION,
     $        MPI_MIN,MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(mulowl,mulow,1,MPI_DOUBLE_PRECISION,
     $        MPI_MIN,MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(muhighl,muhigh,1,MPI_DOUBLE_PRECISION,
     $        MPI_MAX,MPI_COMM_WORLD,ierr)
      else
         stlow  = stlowl
         sthigh = sthighl
         mulow  = mulowl
         muhigh = muhighl
      end if
c
c     Note that all of the nodes know everything they need
c      to know at this point:
c
c     If the gap between stlow and sthigh is large enough and
c      linst is true, we can treat this as an insulator.
c
ctemp
c$$$      if(master) write(*,*) 'High and low eigenvalues',
c$$$     $     nfill,stlow,sthigh
cend temp
      if(linst.and.(sthigh-stlow.gt.twenty*tkb)) then
ctemp
         if(master) write(*,*) 'This crystal is an insulator'
cend temp
c
c        Our work is done here:
c
         fermi = half*(stlow+sthigh)
         return
      end if
c
c     Well, it's not an insulator, so we must do some work:
c
c     Estimated starting value for the search
c
      mu = half*(stlow+sthigh)
c
c     Calculate the fermi sum and its derivative.  We'll use
c      Newton's method where possible to get a good answer.  In
c      the case of a semiconductor or insulator we may have a few
c      problems with Newton, so we'll try something else.
c
 100  nsum = zero
      dnsum = zero
c
c     The Fermi function is non-trivial only for eigen within
c      these ranges
c
      mum20 = mu - tkb*twenty
      mup20 = mu + tkb*twenty
c
c     This saves a cycle on the RS6000
c
      btmu = -beta*mu
c
      do j = 1,jspins
c
c        Version 1.33: again remembe the local k-point:
c
         kloc = 0
c
         do k = istart,ifin
c
            kloc = kloc + 1
c
            nksum = zero
            dnksum = zero
            do i = 1,nv
c
c              Only evaluate the exponentials where needed
c
c              Version 1.33: Eigenvalues are stored locally
c
               if(eigen(i,kloc,j).lt.mum20) then
                  nksum = nksum + one
               else if(eigen(i,kloc,j).lt.mup20) then
                  ex = exp(beta*eigen(i,kloc,j)+btmu)
                  f = one/(one+ex)
                  nksum = nksum + f
c
c                 Note that f'(x) = - e^x f(x)^2
c
                  dnksum = dnksum + ex*f*f
               else
c
c                 We don't need any more eigenvalues at this
c                  k-point, since all will be above mup20
c                  (see Note above)
c
                  go to 200
               end if
            end do
c
c           Add in the k-point weights
c           Version 1.33: k-points are stored globally:
c
 200        nsum  =  nsum + weight(k)*nksum
            dnsum = dnsum + weight(k)*dnksum
         end do
      end do
c
c     Multiply both functions by the electron weight:
c
      nsum = el*nsum
c
c     Remember that dnsum is the derivative wrt mu:
c
      dnsum = beta*el*dnsum
c
c     Sum up over nodes.  Use ALLREDUCE to get the word
c      back to all of the nodes
c
      if(multiprc) then
         call MPI_ALLREDUCE(nsum,nsumg,1,MPI_DOUBLE_PRECISION,
     $        MPI_SUM,MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(dnsum,dnsumg,1,MPI_DOUBLE_PRECISION,
     $        MPI_SUM,MPI_COMM_WORLD,ierr)
c
c        Set the "global" variables back to the local ones,
c        if appropriate
c
         nsum = nsumg
         dnsum = dnsumg
      end if
c
c     Are we finished?  Note that we can answer globally:
c
      if(abs(nsum-zfill).lt.tol) then
         fermi = mu
ctemp
c$$$         write(*,*) 'Found E_Fermi = ',fermi,' on ',MPIRANK
cend temp
         return
      end if
c
c     nsum is a monotonic function of mu, so we can fit the
c      bounds:
c
      if(nsum.lt.zfill) then
         mulow = mu
      else
         muhigh = mu
      end if
c
c     Now use Newton's method.  Note that in many cases the
c      derivative will be zero, so it won't work.  In that case
c      use bisection.  Note that dnsum is non-negative.
c
c     Again, everything here is done on each node.
c
      if(dnsum.lt.tol) then
         mu1 = half*(mulow+muhigh)
      else
         mu1 = mu + (zfill-nsum)/dnsum
         if(mu1.gt.muhigh) then
            mu1 = half*(mu+muhigh)
         else if(mu1.lt.mulow) then
            mu1 = half*(mu+mulow)
         end if
      end if
ctemp
c$$$      if(master) write(*,255) mu,nsum,dnsum,mu1,mulow,muhigh
c$$$      write(*,*) 'Fermi:  MPIRANK = ',MPIRANK
c$$$      write(*,255) mu,nsum,dnsum,mu1,mulow,muhigh
c$$$ 255  format(f13.9,f13.8,f13.7,3f13.9)
cend temp
c
c     Are we going anywhere?
c
      if(muhigh-mulow.gt.tol) then
         mu = mu1
         go to 100
      end if
ctemp
 300  if(master) write(*,305) mu,nsum,zfill
c$$$ 300  write(*,305) mu,nsum,zfill
cend temp
 305  format(/'FERMI:  mu = ',f10.6,' is the best we can do'/
     $     ' Get N = ',f14.8,' wanted',f14.8)
      fermi = mu
      return
      end
