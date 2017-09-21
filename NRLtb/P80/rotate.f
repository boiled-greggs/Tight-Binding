      subroutine rotate(l,m,n,ak,par,ovl,hmat,jsover,smat,ksk,kbas)
      implicit real*8 (a-h,o-z)
c
C*******************************************************
C--->  ROTATES THE TWO-CENTER INTEGRALS GIVEN BY TABLE 1
C--->  OF SLATER=KOSTER TO THE NEW DIRECTION WITH
C--->  DIRECTION COSINES L,M,N
C*******************************************************
c
C........................................................
c
c     Note the mapping onto the spd basis:
c
C........................................................
C---> CODE: S X Y Z YZ ZX XY XX-YY 3*ZZ-RR
C---> CODE: 1 2 3 4  5  6  7   8      9
C........................................................
C---> CODE: SSS SPS PPS PPP SDS PDS PDP DDS DDP DDD PSS DSS DPS DPP
C---> CODE:  1   2   3   4   5   6   7   8   9   10   11   12   13   14
C........................................................
c
c=======================================================================
c
c     REVISION HISTORY:
c     
c-----------------------------------------------------------------------
c
c     Cleaned up as best I could -- mjm  15 Dec 1993
c
c-----------------------------------------------------------------------
c
c     Eliminate the statement functions.
c
c-----------------------------------------------------------------------
c
c     Changed par and ovl to split the onsite and
c      hopping parameters into two arrays -- mjm  14 Sept 1994
c
c-----------------------------------------------------------------------
c
c     In the interest of reducing calls to this function, the
c      onsite terms are now handled in the calling "setxxx"
c      subroutines.                       -- mjm  29 Sept 1994
c
c-----------------------------------------------------------------------
c
c     The Hamiltonian and overlap matricies are now stored
c      as upper-triangle only complex arrays.
c                                         -- mjm   5 Feb  1997
c
c-----------------------------------------------------------------------
c
c     Reorder the parameters on par and ovl for faster execution.
c                                         -- mjm  15 Apr  1997
c
c-----------------------------------------------------------------------
c
c     Version 1.05:
c
c     As part of the general upgrade, jsover is moved into the calling
c      parameters.  Remember that jsover = 1 for non-orthogonal
c      calculations, and 0 for orthogonal calculations.
c
c                                         -- mjm  20 Apr  1998
c-----------------------------------------------------------------------
c
c     Version 1.06:
c
c     jsover has been changed to a logical variable:
c        jsover = .true.  -> Non-orthogonal Hamiltonian (S <> identity)
c        jsover = .false. -> Orthogonal Hamiltonian (S = identity)
c
c                                                mjm --  6 July  1998
c=======================================================================
c
      include 'P1'
c      PARAMETER (MNN=3,MNN1=MNN+1)
c      PARAMETER (MH=9,MATOM=7,MKIND=5,MNN=3,MNN1=MNN+1)
c
c     Yo Ho Ho!  Look out for these:
c
      real*8 l,m,n,ll,mm,nn,lmn,lm,mn,nl
c$$$      common/struc/ nv
c
      logical jsover
c
c     The k-point
c
      real*8 ak(3)
c
c     Onsite and Hopping parameters:
c
      real*8 par(mpkind,mkind,mkind),ovl(mpkind,mkind,mkind)
c
c     Hamiltonian and Overlap matrices
c
      complex*16 hmat(mh*(mh+1)/2),smat(mh*(mh+1)/2)
c
      integer ksk(matom),kbas(mkind)
      common /relat/dlv(3),jkind(2),jatm(2),kneigh
      common /fazon/posq(matom,3)
      dimension p(14),e(9,9),parity(9)
c
c     phase is the complex representation of the phase angle, (cs,sn)
c
      complex*16 phase
c
c     Change these numbers to DP constants
c
c     This is what they should be:
c
      parameter (tpi = 6.283185307179586d0)
      parameter (rt3 = 1.732050807568877d0)
c
      parameter (half  = 5d-1  )
c
      parameter (rt3h = half*rt3)
      parameter (rt3q = half*rt3h)
      parameter (rt32 = 2d0*rt3)
c
      DATA PARITY/1.d0,3*-1.d0,5*1.d0/
c
c_______________________________________________________________________
c
      kk1=ksk(jatm(1))
      kk2=ksk(jatm(2))
c
c     If kk2 < kk1 then we're not going to fill in the matrix
c      anyway: (this is now taken care of in setup and setvol)
c
c$$$      if(kk2.lt.kk1) return
c
      kb1=kbas(jkind(1))
      kb2=kbas(jkind(2))
      lover=0
c
      ll=l*l
      mm=m*m
      nn=n*n
      lm=l*m
      mn=m*n
      nl=n*l
      lmn=lm*n
      dot=ak(1)*(dlv(1)-posq(jatm(2),1)+posq(jatm(1),1))
     $  + ak(2)*(dlv(2)-posq(jatm(2),2)+posq(jatm(1),2))
     $  + ak(3)*(dlv(3)-posq(jatm(2),3)+posq(jatm(1),3))
      cs=cos(tpi*dot)
      sn=sin(tpi*dot)
      phase = dcmplx(cs,sn)
ctemp
c$$$      write(6,3) kk1,kk2,kb1,kb2,sn,cs,dot
cend temp
c
      do j = 1,mpkind
         p(j) = par(j,jkind(2),jkind(1))
      end do
ctemp
c$$$      write(6,2) jatm,jkind,l,m,n
cend temp
C...S,S
  150 E(1,1)=P(1)
C...S,P
      e(1,2)=l*p(2)
      e(1,3)=m*p(2)
      e(1,4)=n*p(2)
C...S,D
      rp5 = rt3*p(5)
      e(1,5)=mn*rp5
      e(1,6)=nl*rp5
      e(1,7)=lm*rp5
      e(1,8)=half*rp5*(ll-mm)
      e(1,9)=(nn-half*(ll+mm))*p(5)
C...P,P
      p34 = p(3)-p(4)
      e(2,2)=p(4)+ll*p34
      e(3,3)=p(4)+mm*p34
      e(4,4)=p(4)+nn*p34
      e(2,3)=lm*p34
      e(2,4)=nl*p34
      e(3,4)=mn*p34
C...P,D
      p67 = rt3*p(6)-2d0*p(7)
      e(2,5)=lmn*p67
      e(2,6)=n*(ll*p67+p(7))
      e(2,7)=m*(ll*p67+p(7))
      e(3,5)=n*(mm*p67+p(7))
      E(3,6)=E(2,5)
      e(3,7)=l*(mm*p67+p(7))
      e(4,5)=m*(nn*p67+p(7))
      e(4,6)=l*(nn*p67+p(7))
      E(4,7)=E(2,5)
      e(2,8) = l*(half*(ll-mm)*p67+p(7))
      p67a = p(6)-rt3*p(7)
      e(2,9) = l*(nn*p67a-(half*(ll+mm))*p(6))
      e(3,8)=m*(half*(ll-mm)*p67-p(7))
      e(3,9) = m*(nn*p67a-(half*(ll+mm))*p(6))
      e(4,8)=n*(half*(ll-mm)*p67)
      e(4,9) = n*(nn*p( 6)-half*(ll+mm)*(p( 6)-rt32*p( 7)))
C...D,D
      p8910 = (3d0*p(8)+p(10))-4d0*p(9)
      e(5,5) = ll*p(10) + mm*(nn*p8910) + (mm+nn)*p(9)
      p910m = p(9) - p(10)
      e(5,6) = lm*(p910m + nn*p8910)
      e(5,7) = nl*(p910m + mm*p8910)
      e(6,6) = mm*p(10) + ll*(nn*p8910) + (ll+nn)*p(9)
      e(6,7) = mn*(p910m + ll*p8910)
      e(7,7) = nn*p(10) + ll*(mm*p8910) + (ll+mm)*p(9)
      e(5,8) = mn*((half*(ll-mm))*p8910-p910m)
      e(6,8) = nl*((half*(ll-mm))*p8910+p910m)
      e(7,8) = lm*(half*(ll-mm))*p8910
      p89 = p(8) - 2d0*p(9)
      p8910a = p89 + p(10)
      p8910b = p8910a + p89
      p89m = p(8) - p(9)
      e(5,9) = rt3*mn*(nn*p89m - (half*(ll+mm))*p8910a)
      e(6,9) = rt3*nl*(nn*p89m - (half*(ll+mm))*p8910a)
      e(7,9) = rt3h*lm*(p(10)-(ll+mm)*p(8) + nn*p8910b)
      e(8,8) = (ll+mm)*p(9) + nn*p(10) + ((half*(ll-mm))**2)*p8910
      e(8,9) = rt3q*(ll-mm)*(p(10) -(ll+mm)*p(8) + nn*p8910b)
      e(9,9) = nn*nn*p(8) + (ll+mm)*(2.5d-1*(ll+mm)*(3d0*p(10)+p(8))
     $     - nn*(p(8)-3d0*p(9)))
c
      DO I=2,9
         DO J=1,I-1
            E(I,J)=parity(i)*(E(J,I)*PARITY(J))
         end do
      end do
      IF (JKIND(1).ne.JKIND(2)) then
C...P,S
         e(2,1)=l*p(11)
         e(3,1)=m*p(11)
         e(4,1)=n*p(11)
C...D,S
         rp12 = rt3*p(12)
         e(5,1)=mn*rp12
         e(6,1)=nl*rp12
         e(7,1)=lm*rp12
         e(8,1)=half*rp12*(ll-mm)
         e(9,1)=(nn-half*(ll+mm))*p(12)
C...D,P
         p134 = rt3*p(13)-2d0*p(14)
         e(5,2)=lmn*p134
         e(6,2)=n*(ll*p134+p(14))
         e(7,2)=m*(ll*p134+p(14))
         e(5,3)=n*(mm*p134+p(14))
         E(6,3)=E(5,2)
         e(7,3)=l*(mm*p134+p(14))
         e(5,4)=m*(nn*p134+p(14))
         e(6,4)=l*(nn*p134+p(14))
         E(7,4)=E(5,2)
         e(8,2)=l*(half*(ll-mm)*p134+p(14))
         p134a = p(13)-rt3*p(14)
         e(9,2) = l*(nn*p134a-(half*(ll+mm))*p(13))
         e(8,3)=m*(half*(ll-mm)*p134-p(14))
         e(9,3) = m*(nn*p134a-(half*(ll+mm))*p(13))
         e(8,4)=n*(half*(ll-mm)*p134)
         e(9,4) = n*(nn*p(13)-half*(ll+mm)*(p(13)-rt32*p(14)))
      end if
C........................................................
      IF(LOVER.ne.1) then
         kkt = (kk2+1)*kk2/2
c
c           There are two cases:  kk1 = kk2 and kk1 < kk2,
c            (kk1 > kk2 is never stored).  If kk1 < kk2, then
c            all the matrix elements will be needed, otherwise
c            only those for i <= j will be needed:
c
         if(kk2.gt.kk1) then
            do j = 1,kb2
               kk2j = kk2+j
               mij = kkt + kk1
               do i = 1,kb1
                  mij = mij + 1
                  hmat(mij) = hmat(mij) + e(i,j)*phase
               end do
               kkt = kkt + kk2j
            end do
         else
c
c           In this case kb1 = kb2
c
            do j = 1,kb2
               kk2j = kk2+j
               mij = kkt + kk1
               do i = 1,j
                  mij = mij + 1
                  hmat(mij) = hmat(mij) + e(i,j)*phase
               end do
               kkt = kkt + kk2j
            end do
         end if
C........................................................
C---> THIS CONTRIBUTION TO S-K MATRIX IS COMPLETE
C........................................................
         IF(.not.JSOVER) RETURN
c
c        Now setup the overlap matrix
c
         lover=1
         do j = 1,mpkind
            p(j) = ovl(j,jkind(2),jkind(1))
         end do
         go to 150
      end if
      kkt = (kk2+1)*kk2/2
      if(kk2.gt.kk1) then
         do j = 1,kb2
            kk2j = kk2+j
            mij = kkt + kk1
            do i = 1,kb1
               mij = mij + 1
               smat(mij) = smat(mij) + e(i,j)*phase
            end do
            kkt = kkt + kk2j
         end do
      else
         do j = 1,kb2
            kk2j = kk2+j
            mij = kkt + kk1
            do i = 1,j
               mij = mij + 1
               smat(mij) = smat(mij) + e(i,j)*phase
            end do
            kkt = kkt + kk2j
         end do
      end if
C........................................................
C---> THIS CONTRIBUTION TO OVERLAP MATRIX IS COMPLETE
C........................................................
      RETURN
 1    FORMAT(' DIAG. TERMS   KIND =',I2,'  BASIS =',I2,
     $     '   SK INDEX =',I3)
 2    FORMAT(' JATM =',2I3,'  JKIND =',2I3,'  L,M,N =',
     $     3F8.5)
 3    FORMAT(' SUB. ROTATE...KK1,KK2,KB1,KB2 =',4I5
     $     ,3X,'SIN,COS=',2F8.4,3X,'DOT=',F8.4)
      END
