      subroutine setpar(jspin,param,jsover,den,flipos,lnoncol,
     $     paros,spinos,hop,over,lprint,ierr)
c
c     This code has been yanked from the original program setup.f.
c     setpar.f is used to calculate the on-site, Slater-Koster, and
c      (if jsover is true) overlap matrix elements.  These results
c      are independent of the k-point.  The routine setham.f is used
c      to turn this data into a Hamiltonian and (if josver.eq.true)
c      overlap matrix.
c
c=======================================================================
c
c     REVISION HISTORY for setup.f:
c
c-----------------------------------------------------------------------
c
c REC 1/20/93
c scaling of parameters with simultaneous fitting of multiple structures
c
c-----------------------------------------------------------------------
c
c REC 6/3/93 uses list generated by search2.f
c
c-----------------------------------------------------------------------
c
c     Includes setup of perturbed Hamiltonian and overlap matrix
c      In this case, we'll use it to calculate the change in
c      eigenvalues as a function of volume.  That is, we're going
c      to calculate the pressure, as well as demonstrate that
c      we can do perturbation theory.
c                                                mjm -- 6 Sept 1994
c
c-----------------------------------------------------------------------
c
c     Note that the onsite parameters in PAR and OVL
c      should really be treated separately from the
c      hopping parameters.  For one thing, the onsites
c      never need jk1 <> jk2.
c                                                mjm -- 14 Sept 1994
c
c-----------------------------------------------------------------------
c
c     Pulled the perturbation initialization out of here
c      it will be installed in a subroutines of the form
c      "setxxx", where "xxx" describes the type of perturbation.
c      For example, setvol sets up the perturbations needed
c      to calculate the pressure.
c                                                mjm -- 19 Sept 1994
c
c-----------------------------------------------------------------------
c
c     Try to reduce the number of exponentials we need to calculate.
c      First, note that the onsite terms need to be calculated only
c      once per structure.  Then note that a similar thing applies
c      to the hopping and overlap parameters.
c                                                mjm -- 27 Sept 1994
c
c-----------------------------------------------------------------------
c
c     Modified to work with many atom types.  This includes the use
c      of a "valence" charge on each atom to help define the density.
c                                                mjm --  3 Oct  1994
c
c-----------------------------------------------------------------------
c
c     Include an "asymptotic" form of the on-site parameters for
c      rho^(2/3) < rhocut(i), i = s,p,d.  This allows the on-site
c      parameters to smoothly approach their atomic values for
c      densities smaller than those obtained in the original fit.
c                                                mjm --  6 Dec  1994
c
c-----------------------------------------------------------------------
c
c     Change from 3 parameters specifying the onsite term to
c      4 parameters                              mjm -- 14 Dec  1994
c
c-----------------------------------------------------------------------
c
c     Split d onsite terms into t2g and eg.  Not Kosher for arbitrarily
c      aligned systems, but it might work OK if we are careful.
c     Also make Slater-Koster prefactors quadratic in R.
c     Yes, we've done all of this before.
c                                                mjm -- 31 Oct  1995
c
c-----------------------------------------------------------------------
c
c     Now (hr,hi) and (sr,si) are stored in the complex arrays
c      hmat and smat.                            mjm --  5 Feb  1997
c
c-----------------------------------------------------------------------
c
c     Change the format of the on-site terms for binary systems.  Now
c      the on-site "densities" will be broken up by atom type, and
c      each atom type will have a separate contribution to the on-site
c      terms.  This eliminates the need for the "charges" zval.
c     Also reverse the indicies on par and ovl so that the mpkind
c      index is first, since this is the most rapidly varying.  This
c      will require a change in rotate.f.
c                                                mjm -- 15 Apr  1997
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
c     Finally eliminated all reference to the "asymptotic" on-site
c      parameters.
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
c     REVISION HISTORY for setpar.f
c
c=======================================================================
c
c     version 1.11:
c
c     Split off from setup.f.  Put on-site, Slater-Koster, and (optional)
c      overlap parameters into the call statement.
c
c     Rewrite the "pars,parp,part2g,pareg" onsite arrays as
c      paros(type,atom), where type=1 for s, 2 for p, etc.  This will
c      allow us to easily transform between "t2g<>eg" and "t2g=eg=d"
c      type so of onsite parametrizations.
c
c                                                     mjm -- 30 Apr 1999
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
c      above is ignored, and we determine the parameters
c
c     paros(i,iat)  = Average of majority and minority
c                     on-site parameters
c     spinos(i,iat) = 1/2 of difference between majority and minority
c                     on-site parameters
c
c     needed for calculating non-collinear magnetization.
c
c     This means that  hup = paros + spinos, and hdown = paros - spinos
c
c     Note that in the non-collinear case we assume that the hopping
c      and overlap SK parameters are the same for spin up and spin
c      down.  Thus we should only execute this routine when jspin = 1.
c      We'll check this before starting the calculation.
c
c                                                     -- mjm 10 Aug 2000
c-----------------------------------------------------------------------
c
c     Minor changes in this version:
c
c     ierr is non-zero if this code tries to execute a forbidden
c      operation.
c     We don't need the /kstuff/ ksk(matom) common block here, so it
c      has been removed.
c                                                     -- mjm 10 Aug 2000
c=======================================================================
c
      implicit real*8 (a-h,o-z)
      include 'P1'
c
c-----------------------------------------------------------------------
c
c     Arguments:
c
c     param holds the tight-binding parameters used to generate the
c      onsite, Slater-Koster Hamiltonian, and overlap parameters.
c      For version 1.20 and above param carries a spin argument
c
      real*8 param(npard,mspin)
c
c     jsover is true if these parameters are non-orthogonal, i.e.,
c      and overlap matrix needs to be generated.
c
      logical jsover
c
c     den(i,j,ispin) is the "density" of atoms of kind i at the jth
c      atom used in calculating the traditional on-site energy.  It is
c      passed out to the calling routine (bande.f) to be used in
c      perturbation parameter setups (setparv.f).  Version 1.30:
c      add spin channels for non-collinear magnetization case.
c      All other calculations can be done with ispin = 1.
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
c     paros(1,iat) is the onsite "s" orbital for the iat-th atom.
c     2 = p, 3 = t2g, 4 = eg.
c
      real*8 paros(4,matom)
c
c     spinos(i,iat) = (paros(i,iat,majority - paros(i,iat,minority))/2
c      This is used only when lnoncol is true, in which case we also
c      have
c      paros(i,iat) = 1/2(paros(i,iat,majority + paros(i,iat,minority))
c
      real*8 spinos(4,matom)
c
c     hop(i,j) is the ith Slater-Koster Hamiltonian parameter
c      for the j-th pair in the list produced by search2.f.
c
      real*8 hop(mpkind,mpair)
c
c     over is the corresponding quantity for the overlap matrix.
c      note that if jsover is false these are not calculated.  We
c      may split these off into a separate subroutine at a later point.
c
      real*8 over(mpkind,mpair)
c
c-----------------------------------------------------------------------
c
      common /codes/posn(matom,3),kkind(matom),natoms
      common /parinfo/ valence(mkind),kinds,kbas(mkind)
      common /types/ npkind(mptype)
      common /relat/dlv(3),jkind(2),jatm(2),kneigh
      common /pairs/tt_list(3,mpair),
     $     dlv_list(3,mpair),dist_list(mpair),
     $     screen_list(mpair),
     $     jkind_list(mpair,2),
     $     jatm_list(mpair,2),
     $     npairs
c
c     realov is the flag for the type of like-atom overlap vector
c     usenew is true if realov is true and the atoms are of the
c      same species
c
      logical realov,usenew,harrison
      common /parcom/ nltype,realov,harrison
c
c     This array holds the R->0 behavior of the like atom overlap
c      matrix (used only if realov is true).
c
      real*8 ovlnull(10)
c
c=======================================================================
c
c     If lprint is true we print out some diagnostic quantities
c
      logical lprint
c
c=======================================================================
c
      parameter (zero = 0d0)
      parameter (half = 5d-1)
      parameter (two  = 2d0)
      parameter (three = 3d0)
      parameter (twtrd = two/three)
c
c=======================================================================
c
c     Data for ovlnull.  Order is
c                  sss, sps, pps, ppp, sds, pds, pdp, dds, ddp, ddd
c
      data ovlnull/1d0, 0d0, 1d0, 1d0, 0d0, 0d0, 0d0, 1d0, 1d0, 1d0/
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
         write(0,*) 'setpar called with lnoncol true and jspin = 2'
         ierr = 1
         return
      end if
c
c     Off-site terms
c
c     Here is where "search" was  spliced in, now uses list
C.......................................................
C--->  DLV(3) IS THE DIRECT LATTICE VECTOR NEEDED TO
C--->  TRANSLATE THE NEIGHBOR AT R2 INTO THE UNIT CELL
C.......................................................
C--->  SEARCH PATIENTLY THROUGH ALL ATOMS IN THE UNIT
C--->  CELL AND IN ALL UNIT CELLS WHICH BORDER IT.
C.......................................................
c
c     Initialize the density matrix for the on-site terms.  Note
c      that this is also needed for the pressure calculation, hence
c      we output them through the call statement.
c
      do IAT=1,natoms
         do kat = 1,kinds
            den(kat,iat,1)=zero
         end do
      end do
      if (lnoncol) then
         do IAT=1,natoms
            do kat = 1,kinds
               den(kat,iat,2)=zero
            end do
         end do
      end if
c
c     Scroll through the current neighbor list, produced by
c      search2.f:
c
      do ipair=1,npairs
c
c        Identify the pair, and find its separation distance:
c
c        Which atoms in the unit cell?
c
         jatm(1) =jatm_list(ipair,1)
         jatm(2) =jatm_list(ipair,2)
         iat=jatm(1)
         jat=jatm(2)
c
c        What kinds of atoms are they?
c
         jkind(1)=jkind_list(ipair,1)
         jkind(2)=jkind_list(ipair,2)
         jk1 = jkind(1)
         jk2 = jkind(2)
c
c        How far apart are they?
c
         dist=dist_list(ipair)
c
c        Parameters are arranged as follows:
c         On-site, 4 for each atom type
c         Off-site, 3 parameters for each pair
c        Because like pairs (11, 22, ...) have the npkind(2) independent
c         parameters, while unlike pairs (21, 31, 32, ...) have npkind(3)
c         parameters, we arrange the parameters as follows
c
c        Onsite:  1, 2, 3, ...
c        Hopping:  11,22,33,...,21,31,32,41,42,43,...
c        Overlap:  11,22,33,...,21,31,32,41,42,43,...
c
c        Note that we have eliminated the "extra" terms which
c         deal with the orientation of the atom in the local
c         environment.  We might want to add these back someday.
c
c-----------------------------------------------------------------------
c
c     ONSITE "DENSITY" CALCULATION:
c
c        For on-site terms determine local "density" and then
c         call rotate (now in setham).
c
c        Now 4 parameters for each onsite term
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
            den(jk2,iat,1) = den(jk2,iat,1) +
     $           exp(-bup*dist)*screen_list(ipair)
            den(jk2,iat,2) = den(jk2,iat,2) +
     $           exp(-bdn*dist)*screen_list(ipair)
            if(jat.ne.iat) then
               ipt = (jkind(1)-1)*(1+npkind(1)*4)+1
               bup = param(ipt,1)*param(ipt,1)
               bdn = param(ipt,2)*param(ipt,2)
               den(jk1,jat,1) = den(jk1,jat,1) +
     $              exp(-bup*dist)*screen_list(ipair)
               den(jk1,jat,2) = den(jk1,jat,2) +
     $              exp(-bdn*dist)*screen_list(ipair)
            end if
         else
c
c           Collinear spin setup.  Note that the spins may be flipped:
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
c           The this is the contribution from atom 2 (jat), of type
c            jkind(2), to the "density" at atom 1 (iat)
c
            dent=Exp(-b * dist)*screen_list(ipair)
c
            den(jk2,iat,1)=den(jk2,iat,1)+dent
c
c           search2 has eliminated double counting, so we have to
c            put the densities in for the reversed case:
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
               den(jk1,jat,1) = den(jk1,jat,1) + dent
            end if
         end if
c
         jk1=min(jkind(1),jkind(2))
         jk2=max(jkind(1),jkind(2))
c
c-----------------------------------------------------------------------
c
c     HOPPING MATRIX ELEMENTS:
c
c        Now initialize the hopping matrix elements.
c         These are independent of the k-point, which is why we have
c         put them in this routine.
c
c        The exact path we take through this is determined by
c         the atom types
c
         if(jk2.eq.jk1) then
c
c           4 onsite terms, cross on-site terms if kinds > 1,
c            and now terms/hopping-overlap ME
c
            ip = kinds*(1+npkind(1)*(3*kinds+1))
     $           + 4*npkind(2)*(jk1-1)
c
c           Like atoms:  npkind(2) parameters
c
            npk = npkind(2)
         else
            ip = kinds*(1+npkind(1)*(3*kinds+1)) +
     $           4*npkind(2)*kinds +
     $           4*npkind(3)*((jk2-1)*(jk2-2)/2+jk1-1)
c
c           Unlike atoms:  npkind(3) parameters
c
            npk = npkind(3)
         end if
c
c        Either way, set up the parameters, which have the form
c         (a + b R + c R^2) Exp(-d^2 R)
c
ctemp
c$$$         if(lprint.and.(dist.lt.7d0)) write(*,503)
c$$$     $        jk1,jk2,ip,ipair,dist,screen_list(ipair)
c$$$ 503     format(4i5,f10.5,g18.8)
cend temp
c
c        In parameter mode 90000 we use a special parameterization for
c         the H atom, affecting only ss_sigma.  This extends
c         the parametrization polynomial by two terms.  Note
c         the futzing of the integrals:
c
         if(nltype.ne.90000) then
c
c           Standard type
c
            do j=1,npk
               abr = param(ip+4*j-3,jspin)+dist*(param(ip+4*j-2,jspin)+
     $              dist*param(ip+4*j-1,jspin))
               p32m = -param(ip+4*j,jspin)**2
               hop(j,ipair) = abr*Exp(p32m*dist)*
     $              screen_list(ipair)
            end do
c
c           Harrison sign convention mode.  So far this is only on
c            for like atoms.  Note that hop(10,...) (dd delta) is
c            not affected
c
            if(harrison.and.(jk2.eq.jk1)) then
               hop(1,ipair) = -abs(hop(1,ipair))
               hop(2,ipair) = +abs(hop(2,ipair))
               hop(3,ipair) = +abs(hop(3,ipair))
               hop(4,ipair) = -abs(hop(4,ipair))
               hop(5,ipair) = -abs(hop(5,ipair))
               hop(6,ipair) = -abs(hop(6,ipair))
               hop(7,ipair) = +abs(hop(7,ipair))
               hop(8,ipair) = -abs(hop(8,ipair))
               hop(9,ipair) = +abs(hop(9,ipair))
            end if
         else
            j = 1
            abr = param(ip+4*j-3,jspin) +
     $           dist*(param(ip+4*j-2,jspin) +
     $           dist*(param(ip+4*j-1,jspin) +
     $           dist*(param(ip+4*j+1,jspin) +
     $           dist*param(ip+4*j+2,jspin))))
            p32m = -param(ip+4*j,jspin)**2
            hop(j,ipair) = abr*Exp(p32m*dist)*
     $           screen_list(ipair)
            do j = 2,npk
               hop(j,ipair) = zero
            end do
         end if
c
         if(lprint.and.(dist.le.8d0))
     $        write(15,'('' Final hopping parameters, dist='',
     $        f10.5,i4,''('',i4,'')'',i4,''('',i4,'')''/
     $        (7f11.7))')
     $        dist,jatm(1),jkind(1),jatm(2),jkind(2),
     $        (hop(j,ipair),j=1,npk)
c
c-----------------------------------------------------------------------
c
c     OVERLAP MATRIX ELEMENTS:
c
c        Non-Orthogonal Basis -- Setup Overlap Matrix
c
         if (jsover) then
c
c           The setup is identical to the one made for the
c            hopping term, except that these parameters are lower
c            on the list
c
c           usenew is true if we use the new type of overlap
c            matrix elements
c
            usenew = realov.and.(jk2.eq.jk1)
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
c           Again, initialize on the first go-round
c
c           Make the same changes for nltype = 90000 as are done
c            in the H integral
c
            if(nltype.ne.90000) then
               do j=1,npk
                  abr = param(iop+4*j-3,jspin)+
     $                 dist*(param(iop+4*j-2,jspin)+
     $                 dist*param(iop+4*j-1,jspin))
                  p32m = -param(iop+4*j,jspin)**2
c
c                 If usenew=.true., use the new expansion of the
c                  overlap parameters to get what we hope is
c                  a better fit.  Note that if jk2=jk1 then
c                  npk = npkind(2) = 4, so this will work.
c
                  if(usenew) abr = ovlnull(j)+dist*abr
c
                  over(j,ipair) =  abr*Exp(p32m*dist)*
     $                 screen_list(ipair)
               end do
c
c              Harrison sign convention mode.  So far this is only on
c               for like atoms.  Note that hop(10,...) (dd delta) is
c               not affected
c
               if(harrison.and.(jk2.eq.jk1)) then
                  over(1,ipair) = +abs(over(1,ipair))
                  over(2,ipair) = -abs(over(2,ipair))
                  over(3,ipair) = -abs(over(3,ipair))
                  over(4,ipair) = +abs(over(4,ipair))
                  over(5,ipair) = +abs(over(5,ipair))
                  over(6,ipair) = +abs(over(6,ipair))
                  over(7,ipair) = -abs(over(7,ipair))
                  over(8,ipair) = +abs(over(8,ipair))
                  over(9,ipair) = -abs(over(9,ipair))
               end if
            else
               j = 1
               abr = param(iop+4*j-3,jspin) +
     $              dist*(param(iop+4*j-2,jspin) +
     $              dist*(param(iop+4*j-1,jspin) +
     $              dist*(param(iop+4*j+1,jspin) +
     $              dist*param(iop+4*j+2,jspin))))
               p32m = -param(iop+4*j,jspin)**2
               over(j,ipair) =  abr*Exp(p32m*dist)*
     $              screen_list(ipair)
               do j = 2,npk
                  over(j,ipair) = zero
               end do
            end if
c
         if(lprint.and.(dist.le.8d0))
     $        write(15,'('' Final overlap parameters, dist='',
     $        f10.5,i4,''('',i4,'')'',i4,''('',i4,'')''/
     $        (7f11.7))')
     $        dist,jatm(1),jkind(1),jatm(2),jkind(2),
     $        (over(j,ipair),j=1,npk)
c
         endif
      end do
c
c-----------------------------------------------------------------------
c
c     CALCUATE ONSITE PARAMETERS FROM DENSITY:
c
      do iat=1,natoms
c
c        Get the information needed about the current atom:
c
         jkind(1)=kkind(iat)
         jkind(2)=kkind(iat)
         jk1 = jkind(1)
         jatm(1)=iat
         jatm(2)=iat
         ip=(jkind(1)-1)*(1+npkind(1)*4)+1
c
         if (lnoncol) then
c
c           The majority and minority on-site terms are determined
c            from the majority and minority densities, and will be
c            rotated later.  The same thing is done in the spin-flip
c            case, below, but the only allowed angle is 180 degrees.
c
c           Remember that we need both spin up and spin down densities:
c
            vm23up = den(jk1,iat,1)**twtrd
            vm23dn = den(jk1,iat,2)**twtrd
            do k = 1,4
               pup = param(ip+1,1)+vm23up*(param(ip+2,1)+
     $              vm23up*(param(ip+3,1)+vm23up*param(ip+4,1)))
               pdn = param(ip+1,2)+vm23dn*(param(ip+2,2)+
     $              vm23dn*(param(ip+3,2)+vm23dn*param(ip+4,2)))
c
c              Find average and difference:
c
               paros (k,iat) = half*(pup + pdn)
               spinos(k,iat) = half*(pup - pdn)
               ip = ip + 4
            end do
c
c           Look at other atom types, if any:
c
            do jk2 = 1,kinds
               if((jk2.ne.jk1).and.
     $              (max(den(jk2,iat,1),den(jk2,iat,2)).gt.zero)) then
                  vm23up = den(jk2,iat,1)**twtrd
                  vm23dn = den(jk2,iat,2)**twtrd
c
                  jkx = max(jk1,jk2)
                  jkm = min(jk1,jk2)
                  ip = kinds*(1+npkind(1)*4) +
     $                 6*npkind(1)*((jkx-1)*(jkx-2)/2+jkm-1)
                  if(jk2.lt.jkx) ip = ip + 3*npkind(1)
                  do k = 1,4
                     pup = param(ip+1,1)+vm23up*(param(ip+2,1)+
     $                    vm23up*(param(ip+3,1)+vm23up*param(ip+4,1)))
                     pdn = param(ip+1,2)+vm23dn*(param(ip+2,2)+
     $                    vm23dn*(param(ip+3,2)+vm23dn*param(ip+4,2)))
c
c                    Find average and difference:
c
                     paros(k,iat)  = half*(pup+pdn)
                     spinos(k,iat) = pup - pdn
                     ip = ip + 4
                  end do
               end if
            end do
         else
c           The assumption here is that the "density" is determined for
c            the spin of the source atoms, and that the on-site parameter
c            is determined according to the spin of the on-site atom.
c            Obviously this is a gross approximation.
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
c            as iat:
c
            vm23 = den(jk1,iat,1)**twtrd
c
c           The "standard" expression is now a cubic polynomial
c            in vm23, rather than quadratic, so there are 4 parameters
c            per onsite term, as promised.
c
c            k = 1 -> s parameters
c            k = 2 -> p parameters
c            k = 3 -> t2g parameters
c            k = 4 -> eg parameters
c
            do k = 1,4
               paros(k,iat) = param(ip+1,jsfrm1)+
     $              vm23*(param(ip+2,jsfrm1)+vm23*(param(ip+3,jsfrm1)+
     $              vm23*param(ip+4,jsfrm1)))
               ip = ip + 4
            end do
            do jk2 = 1,kinds
c
c              We have already done jk2 = jk1:
c
c              Note that not all kinds of atoms need be represented
c               for this part of the calculation.  If they are,
c               then den(jk2,iat,1) <> 0
c
               if((jk2.ne.jk1).and.(den(jk2,iat,1).ne.zero)) then
c
c                 Density contribution
c
                  vm23 = den(jk2,iat,1)**twtrd
c
c                 Find the right indicies
c
                  jkx = max(jk1,jk2)
                  jkm = min(jk1,jk2)
                  ip = kinds*(1+npkind(1)*4) +
     $                 6*npkind(1)*((jkx-1)*(jkx-2)/2+jkm-1)
                  if(jk2.lt.jkx) ip = ip + 3*npkind(1)
c
c                 Now do the expansions:
c
                  do k = 1,4
                     paros(k,iat) = paros(k,iat) +
     $                    vm23*(param(ip+1,jsfrm1)+
     $                    vm23*(param(ip+2,jsfrm1)+
     $                    vm23*param(ip+3,jsfrm1)))
                     ip = ip + 3
                  end do
               end if
            end do
         end if
         if (lprint) then
            write(15,1305)
     $           iat,jkind(1),(jk2,iat,den(jk2,iat,1),jk2=1,kinds)
 1305       format('Atom ',i5,' is of type ',i5,' Densities:'/
     $           (' den(',i5,',',i5,',   up) = ',f15.9,
     $           ' den(',i5,',',i5,',   up) = ',f15.9))
            if(lnoncol) write(15,1315)
     $           (jk2,iat,den(jk2,iat,2),jk2=1,kinds)
 1315       format(' den(',i5,',',i5,', down) = ',f15.9,
     $           ' den(',i5,',',i5,', down) = ',f15.9)
            write(15,'('' Final onsite (spd) parameters, atom='',
     $           i5,/,(4F11.7))')iat,(paros(j,iat),j=1,4)
            if(lnoncol) write(15,'(''Final spin parameters, atom='',
     $           i5,/,(4F11.7))')iat,(spinos(j,iat),j=1,4)
         end if
      end do
c
c-----------------------------------------------------------------------
c
c     And away we go...
c
      return
      end
