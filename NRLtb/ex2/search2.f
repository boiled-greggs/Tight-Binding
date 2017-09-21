      SUBROUTINE search2(npfound,dismin,ierr)
c
c this routine searchs the pairs and saves all required data
c for each pair for each structure
c REC 6/3/93
c
c scaling of parameters with simultaneous fitting of multiple structures
c
c=======================================================================
c
c     REVISION HISTORY:
c
c-----------------------------------------------------------------------
c
c     Recognizing that we only need to fill the upper triangle of the
c      Hamiltonian matrix, and with appropriate changes to setup,
c      we only need to calculate the pairs for which jatom(2) >= jatom(1).
c      This not only speeds up the search, it reduces memory, as we
c      are no longer double-counting pairs.
c
c                                                      -- mjm  2 May 1997
c
c-----------------------------------------------------------------------
c
c     Corrected a minor mistake -- when exiting from the program
c      running under MPI, we should go back to the main program
c      and shut down properly.  Thus we've added ierr to the calling
c      sequence.  When ierr > 0, and error has occurred and we should
c      shut down the run.
c
c-----------------------------------------------------------------------
c
c     Version 1.05:
c
c     Eliminated the "mnn" parameter for dnnsq, and dnn.  Allow screenl
c      to be different for diffent types of interactions (A-A, B-C, etc.)
c      Also screenl and its inverse, scrninv, are no longer defined in
c      parameter statements, but instead are read in input.f
c
c                                            -- mjm   7 April 1998
c
c-----------------------------------------------------------------------
c
c     Version 1.11:
c
c     Added npfound and dismin to the argument list, where they will
c      be used to print diagnositic output from the master processor
c                                            -- mjm   4 August 1999
c=======================================================================
c
      implicit real*8 (a-h,o-z)
c
      include 'P1'
      common /kstuff/ ksk(matom)
      common /codes/posn(matom,3),kkind(matom),natoms
      common /latt1/plv(3,3)
      common /parinfo/ valence(mkind),kinds,kbas(mkind)
      common /cutcom/ dnnsq(mint),dnn(mint),
     $     screenl(mint),scrninv(mint)
      common /relat/dlv(3),jkind(2),jatm(2),kneigh
cMPInote
c     For easier communication, the indicies
c
c$$$      common /displ/nd1,nd2,nd3,kd1,kd2,kd3
c
c     Have been replaced by a single vector:
c
      integer nkdarry(6)
      common/disply/ nkdarry
c
cend MPInote
      common /pairs/tt_list(3,mpair),
     $     dlv_list(3,mpair),dist_list(mpair),
     $     screen_list(mpair),
     $     jkind_list(mpair,2),
     $     jatm_list(mpair,2),
     $     npairs
c
c     Needed for pressure calculations
c
      common /neigh/ dscrn(mpair)
c
      dimension imat(3),t(3),tt(3),rlv(3)
c
      parameter (zero = 0d0)
      parameter (one  = 1d0)
      parameter (five = 5d0)
      parameter (eps  = 1d-3)
c
c=======================================================================
c
C--->  DLV(3) IS THE DIRECT LATTICE VECTOR NEEDED TO
C--->  TRANSLATE THE NEIGHBOR AT R2 INTO THE UNIT CELL
C IT IS NOW IN LATTICE COORDINATES, SO THAT AK (i.e. k-points)
C should be in units of reciprocal lattice coordinates
c
C--->  SEARCH PATIENTLY THROUGH ALL ATOMS IN THE UNIT
C--->  CELL AND IN ALL UNIT CELLS WHICH BORDER IT.
c
c     Start with ierr = 0:
c
      ierr = 0
c
cMPInote
c
c     I suppose this might be done by equivalence statements
c
      nd1 = nkdarry(1)
      nd2 = nkdarry(2)
      nd3 = nkdarry(3)
      kd1 = nkdarry(4)
      kd2 = nkdarry(5)
      kd3 = nkdarry(6)
cend MPInote
      npairs=0
ctemp
c     dismin2 is the minimum distance^2 found in the calculation,
c      independent of atom type
      dismin2 = 1d10
cend temp
      do iat=1,natoms
         jkind(1)=kkind(iat)
         jatm(1) =iat
c$$$         nfirst=0
c$$$         nsecnd=0
c$$$         nthird=0
c
c        Here is the change.  Instead of running the index on
c         jat from 1 to natoms, run it from iat to natoms.  This
c         eliminates the double counting.
c
         do jat=iat,natoms
            jatm(2)=jat
            jkind(2)=kkind(jat)
c
c           Find the interaction type:
c
            if(jkind(1).eq.jkind(2)) then
               intknd = jkind(1)
            else
               jk1 = min(jkind(1),jkind(2))
               jk2 = max(jkind(1),jkind(2))
               intknd = kinds + (jk2-1)*(jk2-2)/2 + jk1
            end if
            screenr=dnn(intknd) - five*screenl(intknd)
C
C--->  DISPLACE ATOM JAT BY LATTICE VECTORS TO
C--->  FIND ALL EQUIVALENT ATOMS IN OTHER UNIT CELLS WHICH
C--->  ARE 0TH, 1ST, 2ND, 3RD NEIGHBORS TO ATOM IAT
c
            do ix=1,kd1
               imat(1)=ix-nd1-1
               do iy=1,kd2
                  imat(2)=iy-nd2-1
                  do 1400 iz=1,kd3
                     imat(3)=iz-nd3-1
                     dsq = zero
                     do i=1,3
                        rlv(i)=zero
                        do j=1,3
                           rlv(i)=rlv(i)+plv(j,i)*
     $                          dble( imat(j) )
                        end do
                        T(I)=POSN(JAT,I)+RLV(I)-
     $                       POSN(IAT,I)
                        dsq = dsq + t(i)*t(i)
                     end do
                     dd=dsq-eps
                     if(dd.gt.dnnsq(intknd)) go to 1400
                     if(dd.lt.dnnsq(intknd)) kneigh=1
                     if(dd.lt.eps)                        kneigh=0
ctemp
                     if(dd.gt.eps) dismin2 = min(dismin2,dsq)
cend temp
                     if(kneigh.eq.0) go to 1400
                     do j=1,3
                        dlv(j)=dble( imat(j) )
                     end do
                     npairs=npairs+1
                     if(npairs.gt.mpair) then
c
c                        Print error message and exit gracefully:
c
                         write(0,*) 'Finding more than ',mpair,
     $                     ' pairs in setup.f'
                         ierr = 1
                         return
                     end if
                     dist=sqrt(dsq)
c$$$                     if(jprndx.ne.0) write(6,9) iat,jat,
c$$$     $                    ksk(iat),ksk(jat),dlv
c$$$                     if(jprdsq.ne.0)write(6,2) iat,jat,imat,t,dist
                     tinv=one/dist
                     do j=1,3
                        tt(j)=t(j)*tinv
                     end do
C
C--->  TT(3) ARE THE DIRECTION COSINES OF THE VECTOR
C--->  SEPARATING ATOM 2 FROM ATOM 1
C
                     expdl = exp((dist-screenr)*scrninv(intknd))
                     screen_list(npairs)=one/(expdl+one)
c
c                    Actually, the derivative of screen_list
c                      wrt dist is dscrn*screen_list
c
                     dscrn(npairs)=
     $                    -expdl*screen_list(npairs)*scrninv(intknd)
                     jkind_list(npairs,1)=jkind(1)
                     jkind_list(npairs,2)=jkind(2)
                     jatm_list(npairs,1)=jatm(1)
                     jatm_list(npairs,2)=jatm(2)
                     dist_list(npairs)=dist
                     do i=1,3
                        dlv_list(i,npairs)=dlv(i)
                        tt_list(i,npairs)=tt(i)
                     end do
 1400             continue
               end do
            end do
         end do
c$$$         if(jprnbr.ne.0)write(6,1) iat,nfirst,nsecnd,nthird
      end do
c$$$      write(15,9200) npairs
c$$$ 9200 format(' Number of pairs = ',i6)
      dismin = sqrt(dismin2)
c$$$      write(15,*) 'Minimum distance this configuration = ',dismin
c
    1 format(2X,'ATOM',I3,5X,3I3,' 1ST,2ND,3RD NEIGHBORS')
    2 format(' PAIR',2I2,' DLV=',3I3,' T=',3F6.2,
     $     ' DIST =',F7.3)
    9 format(' NEIGH',I2,' IAT,JAT',2I2,' INDICES',2I2,
     $     '  DLV=',3F9.5,' DIST=',F9.5)
 26   format(2I3,14F8.5)
c
c     Set npfound = npairs so that we can report what's going on:
c
      npfound = npairs
      return
      end
