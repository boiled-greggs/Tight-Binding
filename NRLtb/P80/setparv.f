      subroutine setparv(jspin,param,thisvol,jsover,den,flipos,lnoncol,
     $     dparos,dspnos,dhop,dover,lprint,ierr)
c
c     This code has been yanked from the original program setvol.f.
c     setparv.f is used to calculate the perturbed on-site, Slater-Koster,
c      and, if jsover is true, overlap matrix elements.  These results
c      are independent of the k-point.  The routine sethamv.f is used to
c      turn this data into a k-point dependent Hamiltonian and, if
c      jsover is true, overlap matrix.
c
c=======================================================================
c
c     REVISION HISTORY for setparv.f
c
c=======================================================================
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
c                                                     mjm -- 23 Jul 1999
c=======================================================================
c
c     version 1.20:
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
c     flipos(i) determines if a spin flip happens for the i^{th} atom.
c      if flipos(i) is set ".false.," then the regular on-site
c      parameters are used.  If, however, flipos(i) is set ".true.,"
c      then we use the spin-2 onsite parameters for jspin = 1, and
c      spin-1 for jspin = 2.
c
c     We now need the entire tight-binding array, so param is now
c      a function of jspin, as it is in bande.f.
c
c     Note that input1.f defaults "flipos(i) = .false." for all atoms.
c      Thus for unpolarized or "standard" calculations we follow the
c      default magnetic case.  Because of this, it is important that
c      FLIPOS IS CHANGED ONLY IN INPUT1.F.
c
c     Note that there are several logical problems with this procedure
c      as applied to the current set of on-site parameters.  We will
c      probably have to revise this later.
c
c                                                    -- mjm  21 Dec 1999
c
c     The "density" for the jspin=1 channel comes from jspin=1 for the
c      other atoms, even if the spins are flipped.   -- mjm  12 Jan 2000
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
c     version 1.30:
c
c     When the variable lnoncol is set true, the flipos information
c      above is ignored, and we determine energy derivatives for the
c      non-collinear magnetization on-site energies, calculating
c
c     dparos(i,iat)  = Average of majority and minority
c                      on-site parameters change with volume
c     dspnos(i,iat) = 1/2 of difference between majority and minority
c                     on-site parameters change with volume
c
c     This means that  dhup = dparos + dspnos,
c      and dhdown = dparos - dspnos
c
c     Note that in the non-collinear case we assume that the hopping
c      and overlap SK parameters are the same for spin up and spin
c      down.  Thus we should only execute this routine when jspin = 1.
c      We'll check this before starting the calculation.
c
c     Note that in the non-collinear case we must put spin up and
c      spin down indices on den and vden.
c
c                                                     -- mjm 14 Aug 2000
c-----------------------------------------------------------------------
c
c     Minor changes in this version:
c
c     ierr is non-zero if this code tries to execute a forbidden
c      operation.
c                                                     -- mjm 14 Aug 2000
c=======================================================================
c
c     REVISION HISTORY for setvol.f
c
c-----------------------------------------------------------------------
c
c     Changed to accomodate multi-atom systems.       -- mjm 29 Sep 1994
c-----------------------------------------------------------------------
c
c     Allows "asymptotic" form of the onsite parameters
c                                                     -- mjm  6 Dec 1994
c-----------------------------------------------------------------------
c
c     Cubic polynomials for onsite terms (4 parameters each)
c                                                     -- mjm 14 Dec 1994
c-----------------------------------------------------------------------
c
c     On-site d terms split, quadratic Slater-Koster prefactors
c                                                     -- mjm 31 Oct 1995
c-----------------------------------------------------------------------
c
c     Replace (dhr,dhi), (dsr,dsi) by dhmat and dsmat -- mjm  5 Feb 1997
c-----------------------------------------------------------------------
c
c     Changed on-site form to match that of setup.f, where different
c      atom types contribute differently to the on-site terms through
c      their densities.                               -- mjm 15 Apr 1997
c-----------------------------------------------------------------------
c
c     search2 has been changed to eliminate the double counting of pairs.
c     This means that we'll have to make sure that the density
c      terms are incremented properly.
c                                                    mjm --  2 May  1997
c-----------------------------------------------------------------------
c
c     Eliminated references to "asymptotic" on-site terms
c                                                    mjm -- 23 Sept 1997
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
c                                                   mjm -- 25 Sept. 1997
c-----------------------------------------------------------------------
c
c     Version 1.05:
c
c     As part of the general upgrade, jsover is moved into the calling
c      parameters.  Remember that jsover = 1 for non-orthogonal
c      calculations, and 0 for orthogonal calculations.
c                                                   mjm -- 20 April 1998
c-----------------------------------------------------------------------
c
c     Removed the "labels" common block, and pass the volume ("thisvol")
c      via the calling parameters.
c                                                   mjm -- 22 April 1998
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
c                                                   mjm --  4 Aug   1998
c=======================================================================
c
      implicit none
      include 'P1'
c
c=======================================================================
c
c     Arguments
c
c     jspin is the current spin of the overall lattice:
c
      integer jspin
c
c     param holds the tight-binding parameters used to generatethe
c      onsite, Slater-Koster, and overlap parameters
c
      real*8 param(npard,mspin)
c
c     thisvol is the current unit cell volume, used in calculating
c      volume changes
c
      real*8 thisvol
c
c     jsover is true for a non-orthogonal calculation
c
      logical jsover
c
c     The onsite densities, from setpar.f, passed through bande.f.
c      Non-collinear magnetization uses both channels:
c
      real*8 den(mkind,matom,mspin)
c
c     As note above, flipos(i) is true if atom i is "flipped"
c      relative to our common notion of the up spin direction.
c
      logical flipos(matom)
c
c     lnoncol is true when we are implementing the "non-collinear"
c      magnetization option.  See the discription above in the
c      Version 1.30 notes.  lnoncol = .true. overrides the
c      flipos settings.
c
      logical lnoncol
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
c     dspnos(i,iat) = (dparos(i,iat,majority - dparos(i,iat,minority))/2
c      This is used only when lnoncol is true, in which case we also
c      have
c      dparos(i,iat) = 1/2(dparos(i,iat,majority + dparos(i,iat,minority))
c
      real*8 dspnos(4,matom)
c
      real*8 dhop(mpkind,mpair),dover(mpkind,mpair)
c
c     In Harrison sign mode we need the canonical signs of the
c      Slater-Koster integrals (see data statement below)
c
      real*8 harsign(9)
c
c     lprint is true if diagnostics are to be printed
c
      logical lprint
c
c     ierr signals error exits back to the calling routine.  "0" is
c      a normal exit.
c
      integer ierr
c
c=======================================================================
c
c     Real Arrays
c
c     This is the perturbation contribution from atoms of type mkind
c      to the density at atom matom.  For non-collinear calculations
c      we will need to hold both spin up and spin down variables:
c
      real*8 vden(mkind,matom,mspin)
c
c     This array holds the R->0 behavior of the like atom overlap
c      matrix.
c
      real*8 ovlnull(10)
c
c-----------------------------------------------------------------------
c
c     Real scalars
c
c     dist is the distance of a local pair
c
      real*8 dist
c
c     vdist is the change in distance due to a volume change
c
      real*8 vdist
c
c     b is the  exponential parameter squared (non-collinear case)
c     bup is the onsite majority spin exponential parameter squared
c     bdn is the onsite minority spin exponential parameter squared
c     dent is the local contribution to the on-site "atomic density"
c     drho is the change in density at the local site
c     drhoup and drhodn are the spin generalizations of drho
c     dpup and dpdn are the changes in the majority and minority
c      onsite spin parameters
c
      real*8 b,bup,bdn,dent,drho,drhoup,drhodn,dpup,dpdn
c
c     Temporary markers for various parameter combinations
c
      real*8 p32m,par21,abr,ovl21,vm13,vm23,vm2t,vm3t
c
c     Spin polarized versions of the above:
c
      real*8 vm13up,vm23up,vm2tup,vm3tup
      real*8 vm13dn,vm23dn,vm2tdn,vm3tdn
c
c-----------------------------------------------------------------------
c
c     Integers
c
c     jsfrm1 and jsfrm2 indicate the local spin of certain atoms:
c
      integer jsfrm1,jsfrm2
c
c     Miscellaneous counters and pointers:
c
      integer jk1,jk2,iat,iop,ip,ipt,j,jat,jkm,jkx,k,kat,npk
c
c-----------------------------------------------------------------------
c
c     Counter names:
c
      integer ipair
c
c-----------------------------------------------------------------------
c
c     Logical parameters
c
c     usenew is true if we use the "new" form for on-site parameterization
c     usehar is true if we are currently using the Harrison sign
c      convention
c
      logical usenew,usehar
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
      real*8 dscrn(mpair)
      common /neigh/ dscrn
c
c     realov is the flag for the type of like-atom overlap vector
c
      integer nltype
      logical realov,harrison
      common /parcom/ nltype,realov,harrison
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
c=======================================================================
c
      real*8 zero,one,two,three,four,third,twtrd,threeh,half
c
      parameter (zero = 0d0)
      parameter (one  = 1d0)
      parameter (two  = 2d0)
      parameter (three = 3d0)
      parameter (four = 4d0)
      parameter (third = one/three)
      parameter (twtrd = two/three)
      parameter (threeh = three/two)
      parameter (half = 5d-1)
c
c=======================================================================
c
c     Data for ovlnull.  Order is
c                  sss, sps, pps, ppp, sds, pds, pdp, dds, ddp, ddd
c
      data ovlnull/1d0, 0d0, 1d0, 1d0, 0d0, 0d0, 0d0, 1d0, 1d0, 1d0/
c
c     Data for harsign
c
      data harsign/-1d0,1d0,1d0,-1d0,-1d0,-1d0,1d0,-1d0,1d0/
c
c=======================================================================
c
c     Set exit (non-)error code:
c
      ierr = 0
c
c     Check that we are not calling for non-collinear calculations
c      with jspin = 2:
c
      if(lnoncol.and.(jspin.eq.2)) then
         write(0,*) 'setparv called with lnoncol true and jspin = 2'
         ierr = 1
         return
      end if
c
c     Off-site terms
c
      do iat=1,natoms
         do kat = 1,kinds
            vden(kat,iat,1) = zero
         end do
      end do
      if (lnoncol) then
         do IAT=1,natoms
            do kat = 1,kinds
               vden(kat,iat,2)=zero
            end do
         end do
      end if
c
      kneigh=1
      do ipair=1,npairs
         jkind(1)=jkind_list(ipair,1)
         jkind(2)=jkind_list(ipair,2)
         jk1 = jkind(1)
         jk2 = jkind(2)
         jatm(1) =jatm_list(ipair,1)
         jatm(2) =jatm_list(ipair,2)
         iat=jatm(1)
         jat=jatm(2)
         dist=dist_list(ipair)
c
c        Changing the volume by V -> V + dV changes distances
c         by (1+dV/V)^(1/3) = 1 + dV/(3V) + O(dV^2)
c
         vdist = third*dist/thisvol
c
c        For on-site terms determine local "density" and its
c         perturbation and then call rotate below
c
         ip=(jkind(2)-1)*(1+npkind(1)*4)+1
c
         if(lnoncol) then
c
c           Non-collinear spin setup.  We need to store both majority
c            and minority spin densities simultaneously.  We'll denote
c            these by "up" and "dn" or "1" and "2", respectively,
c            even though the local direction might be pointing to
c            Topeka.
c
            bup = param(ip,1)*param(ip,1)
            bdn = param(ip,2)*param(ip,2)
c
c           Do not forget the lack of double-counting:
c
            vden(jk2,iat,1) = vden(jk2,iat,1) +
     $           vdist*(dscrn(ipair)-bup)*
     $           exp(-bup*dist)*screen_list(ipair)
            vden(jk2,iat,2) = vden(jk2,iat,2) +
     $           vdist*(dscrn(ipair)-bdn)*
     $           exp(-bdn*dist)*screen_list(ipair)
            if(jat.ne.iat) then
               ipt = (jkind(1)-1)*(1+npkind(1)*4)+1
               bup = param(ipt,1)*param(ipt,1)
               bdn = param(ipt,2)*param(ipt,2)
               vden(jk1,jat,1) = vden(jk1,jat,1) +
     $           vdist*(dscrn(ipair)-bup)*
     $              exp(-bup*dist)*screen_list(ipair)
               vden(jk1,jat,2) = vden(jk1,jat,2) +
     $           vdist*(dscrn(ipair)-bdn)*
     $              exp(-bdn*dist)*screen_list(ipair)
            end if
         else
c
c           Non-collinear spin setup.
c
c           What we want (he thinks) is for the jspin=1 channel of atom
c            iat to get its marching orders (e.g., the density) from the
c            jspin=1 channel of atom jat, regardless of the relative
c            orientation of the two atoms.  (Take 1->2 and/or iat<->jat
c            in the last sentence as needed.)  Thus the spin flip check
c            <<should>> be:
c
            if(flipos(iat)) then
c
c              This little ditty depends upon the fact that jspin is
c               either 1 or 2, and nothing else:
c
               jsfrm2 = 3 - jspin
            else
               jsfrm2 = jspin
            end if
c
            b=param(ip,jsfrm2)*param(ip,jsfrm2)
c
c           The this is the contribution from atom 2 (jat) to
c            the "density" at atom 1 (iat)
c
            dent=Exp(-b * dist)*screen_list(ipair)
            vden(jk2,iat,1) = vden(jk2,iat,1) +
     $           (dscrn(ipair)-b)*dent*vdist
c
c           Now for the reversed pair:
c
            if(jat.ne.iat) then
c
c              In line with the above, we should have
c
               if(flipos(jat)) then
                  jsfrm2 = 3 - jspin
               else
                  jsfrm2 = jspin
               end if
c
               ipt = (jkind(1)-1)*(1+npkind(1)*4)+1
               b = param(ipt,jsfrm2)*param(ipt,jsfrm2)
               dent = Exp(-b*dist)*screen_list(ipair)
               vden(jk1,jat,1) = vden(jk1,jat,1) +
     $              (dscrn(ipair)-b)*dent*vdist
            end if
         end if
c
         jk1=min(jkind(1),jkind(2))
         jk2=max(jkind(1),jkind(2))
c
c        See setup.f for parameter ordering and program logic
c
         if(jk2.eq.jk1) then
            ip = kinds*(1+npkind(1)*(3*kinds+1))
     $           + 4*npkind(2)*(jk1-1)
            npk = npkind(2)
         else
            ip = kinds*(1+npkind(1)*(3*kinds+1)) +
     $           4*npkind(2)*kinds +
     $           4*npkind(3)*((jk2-1)*(jk2-2)/2+jk1-1)
            npk = npkind(3)
         end if
c
c        Check if we are using the Harrison form
c
         usehar = harrison.and.(jk2.eq.jk1)
c
c        Note the special extension for H_{ss sigma}:
c
         if(nltype.ne.90000) then
            do j=1,npk
               abr = param(ip+4*j-3,jspin)+
     $              dist*(param(ip+4*j-2,jspin)+
     $              dist*param(ip+4*j-1,jspin))
               p32m = -param(ip+4*j,jspin)**2
c
c        par21 is the actual hopping term for this ME
c
               par21= abr*Exp(p32m*dist)*
     $              screen_list(ipair)
               dhop(j,ipair) = vdist*(((param(ip+4*j-2,jspin)+
     $              two*dist*param(ip+4*j-1,jspin))+p32m*abr)*
     $              Exp(p32m*dist)*
     $              screen_list(ipair)+dscrn(ipair)*par21)
c
c              In the Harrison form, get the right sign:
c
               if(usehar.and.j.lt.npk) dhop(j,ipair) = harsign(j)*
     $              dhop(j,ipair)*sign(one,abr)
c
            end do
         else
            j = 1
            abr = param(ip+4*j-3,jspin) +
     $           dist*(param(ip+4*j-2,jspin) +
     $           dist*(param(ip+4*j-1,jspin) +
     $           dist*(param(ip+4*j+1,jspin) +
     $           dist*param(ip+4*j+2,jspin))))
            p32m = -param(ip+4*j,jspin)**2
            par21= abr*Exp(p32m*dist)*
     $           screen_list(ipair)
            dhop(j,ipair) = vdist*(((param(ip+4*j-2,jspin) +
     $           dist*(two*param(ip+4*j-1,jspin) +
     $           dist*(three*param(ip+4*j+1,jspin) +
     $           dist*(four*param(ip+4*j+2,jspin))))) +
     $           p32m*abr)*Exp(p32m*dist) *
     $           screen_list(ipair)+dscrn(ipair)*par21)
            do j = 2,npk
               dhop(j,ipair) = zero
            end do
         end if
         if(lprint.and.dist.le.8d0) then
            write(15,'('' Final hopping derivatives'',
     $           f10.5,2i4,/(7f11.7))')
     $           dist,iat,jat,(dhop(j,ipair),j=1,npk)
         end if
c
c        Overlap matrices
c
         if(jsover)then
c
c           Do we use the new form for onsite integrals?
c
            usenew=realov.and.(jk2.eq.jk1)
c
            if(jk2.eq.jk1) then
               iop = kinds*(1+npkind(1)*(3*kinds+1)) +
     $              4*npkind(2)*kinds +
     $              4*npkind(3)*(kinds*(kinds-1)/2) +
     $              4*npkind(2)*(jk1-1)
               npk = npkind(2)
            else
               iop = kinds*(1+npkind(1)*(3*kinds+1)) +
     $              8*npkind(2)*kinds +
     $              4*npkind(3)*(kinds*(kinds-1)/2) +
     $              4*npkind(3)*((jk2-1)*(jk2-2)/2+jk1-1)
               npk = npkind(3)
            end if
c
c           Modify for Hydrogen ss_sigma:
c
            if(nltype.ne.9000) then
               if(usenew) then
c
c                 The new form for on-site integrals:
c
                  do j=1,npk
                     abr = ovlnull(j)+dist*(
     $                    param(iop+4*j-3,jspin)+
     $                    dist*(param(iop+4*j-2,jspin)+
     $                    dist*param(iop+4*j-1,jspin)))
                     p32m = -param(iop+4*j,jspin)**2
                     ovl21 = abr*Exp(p32m*dist)*
     $                    screen_list(ipair)
                     dover(j,ipair) = vdist*(((param(iop+4*j-3,jspin)+
     $                    two*dist*param(iop+4*j-2,jspin)+
     $                    three*dist*dist*param(iop+4*j-1,jspin))+
     $                    p32m*abr)*Exp(p32m*dist)*
     $                    screen_list(ipair)+dscrn(ipair)*ovl21)
c
c                    Note the Harrison sign for overlap is opposite
c                     that for Hopping:
c
                     if(usehar.and.j.lt.npk) dover(j,ipair) =
     $                    -harsign(j)*dover(j,ipair)*sign(one,abr)
                  end do
               else
                  do j=1,npk
                     abr = param(iop+4*j-3,jspin)+
     $                    dist*(param(iop+4*j-2,jspin)+
     $                    dist*param(iop+4*j-1,jspin))
                     p32m = -param(iop+4*j,jspin)**2
                     ovl21 = abr*Exp(p32m*dist)*
     $                    screen_list(ipair)
                     dover(j,ipair) = vdist*(((param(iop+4*j-2,jspin)+
     $                    two*dist*param(iop+4*j-1,jspin))+p32m*abr)*
     $                    Exp(p32m*dist)*
     $                    screen_list(ipair)+dscrn(ipair)*ovl21)
c
c                    Note the Harrison sign for overlap is opposite
c                     that for Hopping:
c
                     if(usehar.and.j.lt.npk) dover(j,ipair) =
     $                    -harsign(j)*dover(j,ipair)*sign(one,abr)
                  end do
               end if
            else
               j = 1
               abr = param(iop+4*j-3,jspin) +
     $              dist*(param(iop+4*j-2,jspin) +
     $              dist*(param(iop+4*j-1,jspin) +
     $              dist*(param(iop+4*j+1,jspin) +
     $              dist*param(iop+4*j+2,jspin))))
               p32m = -param(iop+4*j,jspin)**2
               par21= abr*Exp(p32m*dist)*
     $              screen_list(ipair)
               dover(j,ipair) = vdist*(((param(iop+4*j-2,jspin) +
     $              dist*(two*param(iop+4*j-1,jspin) +
     $              dist*(three*param(iop+4*j+1,jspin) +
     $              dist*(four*param(iop+4*j+2,jspin))))) +
     $              p32m*abr)*Exp(p32m*dist) *
     $              screen_list(ipair)+dscrn(ipair)*par21)
               do j = 2,npk
                  dover(j,ipair) = zero
               end do
            end if
            if(lprint.and.dist.le.8d0) then
               write(15,'('' Final overlap derivatives''/(7f11.7))')
     $              (dover(j,ipair),j=1,npk)
            end if
         end if
      end do
c
c     On-site term setup
c
      do iat=1,natoms
c
c        Describe the current atom:
c
         jkind(1)=kkind(iat)
         jkind(2)=kkind(iat)
         jk1 = jkind(1)
         jatm(1)=iat
         jatm(2)=iat
         ip=(jkind(1)-1)*(1+npkind(1)*4)+1
c
         if(lnoncol) then
c
c           On-site contributions from atoms of the same type
c
            vm13up = den(jk1,iat,1)**third
            vm23up = vm13up*vm13up
            vm2tup = two*vm23up
            vm3tup = threeh*vm23up
            vm13dn = den(jk1,iat,2)**third
            vm23dn = vm13dn*vm13dn
            vm2tdn = two*vm23dn
            vm3tdn = threeh*vm23dn
c s
            drhoup = twtrd*vden(jk1,iat,1)/vm13up
            drhodn = twtrd*vden(jk1,iat,2)/vm13dn
            do k = 1,4
               dpup = drhoup*(param(ip+2,1)+
     $              vm2tup*(param(ip+3,1)
     $              +vm3tup*param(ip+4,1)))
               dpdn = drhodn*(param(ip+2,2)+
     $              vm2tdn*(param(ip+3,2)
     $              +vm3tdn*param(ip+4,2)))
               dparos(k,iat) = half*(dpup+dpdn)
               dspnos(k,iat) = half*(dpup-dpdn)
               ip = ip + 4
            end do
c
c           On-site positions for atoms of a different type
c
            do jk2 = 1,kinds
c
c              Note that not all kinds of atoms need be represented
c               for this part of the calculation.  If they are,
c               then den(jk2,iat,1) <> 0
c
               if((jk2.ne.jk1).and.(den(jk2,iat,1).ne.zero)) then
c
c                 New density contribution
c
                  vm13up = den(jk2,iat,1)**third
                  vm23up = vm13up*vm13up
                  vm2tup = two*vm23up
                  vm3tup = threeh*vm23up
                  drhoup = twtrd*vden(jk2,iat,1)/vm13up
                  vm13dn = den(jk2,iat,2)**third
                  vm23dn = vm13dn*vm13dn
                  vm2tdn = two*vm23dn
                  vm3tdn = threeh*vm23dn
                  drhodn = twtrd*vden(jk2,iat,1)/vm13dn
c
c                 Find the right indicies
c
                  jkx = max(jk1,jk2)
                  jkm = min(jk1,jk2)
                  ip = kinds*(1+npkind(1)*4) +
     $                 6*npkind(1)*((jkx-1)*(jkx-2)/2+jkm-1)
                  if(jk2.lt.jkx) ip = ip + 3*npkind(1)
c
                  do k = 1,4
                     dpup = drhoup*(param(ip+2,1)+
     $                    vm2tup*(param(ip+3,1)
     $                    +vm3tup*param(ip+4,1)))
                     dpdn = drhodn*(param(ip+2,2)+
     $                    vm2tdn*(param(ip+3,2)
     $                    +vm3tdn*param(ip+4,2)))
                     dparos(k,iat) = half*(dpup+dpdn)
                     dspnos(k,iat) = half*(dpup-dpdn)
                     ip = ip + 3
                  end do
               end if
            end do
         else
c
c           The assumption here is that the "density" is determined for
c            the spin of the source atoms, and that the on-site
c            parameter is c determined according to the spin of the
c            on-site atom.  Obviously c this is a gross approximation.
c
c           So test to see if atom iat flips the spin, and assign
c            the appropriate fake value of jspin:
c
            if(flipos(iat)) then
c
c              Again, we rely on the observation that 0 < jspin < 3:
c
               jsfrm1 = 3 - jspin
            else
               jsfrm1 = jspin
            end if
c
c           On-site contributions from atoms of the same type
c
            vm13 = den(jk1,iat,1)**third
            vm23= vm13*vm13
            vm2t = two*vm23
c
c           vm3t = 1.5 vm23
c
            vm3t = threeh*vm23
c s
            drho = twtrd*vden(jk1,iat,1)/vm13
c
c           4 terms in the onsite expansion, now
c
c           Note that
c
c              k = 1 -> s parameters
c              k = 2 -> p parameters
c              k = 3 -> t2g parameters
c              k = 4 -> eg parameters
c
            do k = 1,4
               dparos(k,iat) = drho*(param(ip+2,jsfrm1)+
     $              vm2t*(param(ip+3,jsfrm1)
     $              +vm3t*param(ip+4,jsfrm1)))
               ip = ip + 4
            end do
c
c           On-site positions for atoms of a different type
c
            do jk2 = 1,kinds
c
c              Note that not all kinds of atoms need be represented
c               for this part of the calculation.  If they are,
c               then den(jk2,iat,1) <> 0
c
               if((jk2.ne.jk1).and.(den(jk2,iat,1).ne.zero)) then
c
c                 New density contribution
c
                  vm13 = den(jk2,iat,1)**third
                  vm23 = vm13*vm13
                  vm2t = two*vm23
                  vm3t = threeh*vm23
                  drho = twtrd*vden(jk2,iat,1)/vm13
c
c                 Find the right indicies
c
                  jkx = max(jk1,jk2)
                  jkm = min(jk1,jk2)
                  ip = kinds*(1+npkind(1)*4) +
     $                 6*npkind(1)*((jkx-1)*(jkx-2)/2+jkm-1)
                  if(jk2.lt.jkx) ip = ip + 3*npkind(1)
c
                  do k = 1,4
                     dparos(k,iat) = dparos(k,iat) +
     $                    drho*(param(ip+1,jsfrm1) +
     $                    vm2t*(param(ip+2,jsfrm1) +
     $                    vm3t*param(ip+3,jsfrm1)))
                     ip = ip + 3
                  end do
               end if
            end do
         end if
         if(lprint)then
            write(15,'('' Final onsite perturbations, atom='',
     $           i5,/,(4F11.7))')iat,(dparos(k,iat),k=1,4)
            if(lnoncol) write
     $           (15,'('' Final onsite spin perturbations, atom='',
     $           i5,/,(4F11.7))')iat,(dspnos(k,iat),k=1,4)
         end if
      end do
      return
      end
