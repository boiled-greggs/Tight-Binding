      subroutine kptin
c
c=======================================================================
c
c     REVISION HISTORY:
c
c-----------------------------------------------------------------------
c
c     Version 1.06:
c
c     jsover has been changed to a logical variable:
c        jsover = .true.  -> Non-orthogonal Hamiltonian (S <> identity)
c        jsover = .false. -> Orthogonal Hamiltonian (S = identity)
c
c     Moved k-point generation parameters into 'P3'
c
c                                                mjm --  6 July  1998
c
c     Got rid of iwrite in /lia/.                  mjm --  4 Aug   1998
c
c     Version 1.20:
c
c     Got rid of jspins in /lia/.                  mjm -- 21 Dec   1999
c
c     Fixed bug which printed statement 9005 twice.   mjm -- 16 Aug 2000
c
c=======================================================================
c
      IMPLICIT REAL*8 (A-H,O-Z)
      save nk
      LOGICAL LSPEC
      logical ILAT,NEWSYM,LINFOK
      INCLUDE 'P1'
      include 'P3'
      parameter (EPS=1d-6)
      parameter (ONE=1d0)
      parameter (ZERO=0d0)
C
      CHARACTER*255 SGRPFIL
      common/cd1/ x(nkd),y(nkd),z(nkd)
c$$$      common/codes/posn(matom,3),kkind(matom),natoms
c
c     See note about kinds below:
c
c$$$      common /parinfo/ valence(mkind),kinds,kbas(mkind)
c
      common/latt1/plv(3,3)
      logical jsover
      common /lia/ npts,jsover
      COMMON /LATTYP/ AVEC(3,3),BVEC(3,3),BIJ(3,3),WSVOL,LATTIC
      COMMON /SPACE / RSPC(3,3,NOPD),RSPCT(3,3,NOPD),TNPSPC(3,NOPD),NOP
      common /kmesh1/ qkpt(3,nkd,nwdd),wghtk(nkd,nwdd),nkv(nwdd)
      COMMON /KMESH2/ KPDIV1(NWDD),KPDIV2(NWDD),KPDIV3(NWDD),
     $     LSPEC(NWDD)
c
c -- Begin mjm --
c
c     If the number of k-points, NKV(M), = -1313, then read the rest
c      of the k-point information from a file
c
      character*255 kfile
c
c -- End mjm --
c
      real*8 A(9),B(9),EL(4)
C
      parameter (pi = 3.14159265358979323846d0)
      parameter (tpi = 2d0*pi)
c
      WGHTS=0d0
c
c     There is only one window in the calculation
c
      M=1
c
c NEW SYMMETRY TYPE?
ctemp
c$$$      write(0,*) 'Ready to read NEWSYM'
cend temp
      READ(10,'(7x,L1)')NEWSYM
ctemp
c$$$      write(0,'(''NEWSYM = '',l1)') newsym
cend temp
c
c     Defining kio keeps us from trying to close an unopened file
c      if newysm is false:
c
      kio = 0
c
      if(newsym)then
         read(10,'(7x,I1)')LATTIC
C     COMPUTE RECIPROCAL LATTICE GENERATORS:
C
         L=0
         DO J=1,3
            DO I=1,3
               L=L+1
               AVEC(I,J)=PLV(J,I)
               A(L)=AVEC(I,J)
            end do
         end do
         B(1)=A(5)*A(9)-A(6)*A(8)
         B(2)=A(6)*A(7)-A(4)*A(9)
         B(3)=A(4)*A(8)-A(5)*A(7)
         B(4)=A(8)*A(3)-A(9)*A(2)
         B(5)=A(9)*A(1)-A(7)*A(3)
         B(6)=A(7)*A(2)-A(8)*A(1)
         B(7)=A(2)*A(6)-A(3)*A(5)
         B(8)=A(3)*A(4)-A(1)*A(6)
         B(9)=A(1)*A(5)-A(2)*A(4)
         WSVOL=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
         L=0
         prefact = tpi/wsvol
         DO J=1,3
            DO I=1,3
               L=L+1
               BVEC(I,J)=prefact*B(L)
            end do
         end do
         WSVOL=ABS(WSVOL)
         WRITE (6,9004) WSVOL
         WRITE(15,9004) WSVOL
         WRITE (6,9001)
         write(15,9001)
         WRITE (6,9002) (A(I),I=1,9)
         write(15,9002) (A(I),I=1,9)
         WRITE (6,9005)
         write(15,9005)
         WRITE (6,9002) ((BVEC(I,J),I=1,3),J=1,3)
         write(15,9002) ((BVEC(I,J),I=1,3),J=1,3)
 9001    FORMAT (/2X,'PRIMITIVE VECTORS OF THE DIRECT LATTICE:'/)
 9002    FORMAT (3F10.6)
 9003    FORMAT (27X,3D18.8)
 9004    FORMAT (/2X,'VOLUME OF THE UNIT CELL=',F16.8)
 9005    FORMAT (/2X,'PRIMITIVE VECTORS OF THE RECIPROCAL LATTICE:'/)
 9006    FORMAT (//30X,'RECIPROCAL SPACE METRIC:'/)
 9007    FORMAT (I3,A77)
 9008    FORMAT (//10X,'LATTIC = ',I2,2X,A77/)
 9009    FORMAT (/10X,'HEX. LATT. CONSTANTS: A= ',F10.6,'  C= ',F10.6/)
C
C     CALCULATE THE METRIC IN RECIPROCAL SPACE:
C
         DO I=1,3
            DO J=1,3
               BIJ(J,I)=ZERO
            end do
         end do
         DO J=1,3
            DO I=1,3
               DO K=1,3
                  BIJ(I,J)=BVEC(K,I)*BVEC(K,J)+BIJ(I,J)
               end do
            end do
         end do
         WRITE (6,9006)
         write(15,9006)
         WRITE (6,9003) ((BIJ(I,J),I=1,3),J=1,3)
         write(15,9003) ((BIJ(I,J),I=1,3),J=1,3)
c$$$      do i=1,NATOMS
c$$$        do J=1,3
c$$$           POS(J,I)=POSN(I,J)
c$$$        end do
c$$$      end do
c
c     this is just to fake out rdspgp so it generates atoms with
c     symmetry.  you NEED to fix this so that the different types get
c     listed correctely in KKIND you need to fix POSQ also.
c     Alternatively, just don't use the generated atoms in kptgen
c     right now until you start doing complicated structures
c
c$$$      NTYPE=1
c$$$      do 930 i=1,kinds
c$$$         NEQ(i)=1
c$$$ 930  continue
C SPACE GROUP FILE
         READ(10,'(A)')SGRPFIL
         READ(10,'(5x,L1)')ILAT
         OPEN(25,FILE=SGRPFIL,STATUS='OLD')
c$$$      CALL RDSPGP(.TRUE.,ILAT,25,NAT)
         CALL RDSPGP(.TRUE.,ILAT,25)
         CLOSE(25)
c     DO 940 I=1,NATOMS
c     DO 935 J=1,3
c            POSN(I,J)=POS(J,I)
c 935        CONTINUE
c 940        CONTINUE
c      CALL MULTATOMS
         READ  (10,1000) NKV(M),SFK
c
c -- Begin mjm --
c
         if(nkv(m).eq.-1313) then
c
c        IF NKV(M) has the "magic" number -1313, then read all of the
c         k-point information from the file named in the next line
c
            kio = 33
            read(10,1313) kfile
 1313       format(a255)
            open(unit=kio,file=kfile,status='old',err=13000)
            read (kio,1000) nkv(m),sfk
         else
c
c        Otherwise read the k-point information directly from this file:
c
c
c        Note that "READ(10" is replace by "READ(KIO" beyond this point
c
            kio = 10
         end if
c
c -- End mjm --
c
         IF(SFK.LE.EPS) SFK=ONE
C
c$$$      IF (TRIA .AND. NKV(M).GT.0)
c$$$     *STOP 'PTIOE2: TRIA = .TRUE., SO NKV(M) MUST BE .LE. 0'
C
         IF (NKV(M).GT.0) THEN
C
c$$$      LKPTGN=.FALSE.
C
            IF(NKV(M).EQ.999) THEN
C
               READ(kio,1000) KPDIR
               NK=0
               DO KP=1,KPDIR
                  READ(kio,*)  KD,(QKPT(I,NK+1,M),I=1,3),
     *                 (QKPT(I,NK+KD,M),I=1,3)
                  DO I=1,3
                     EL(I)=(QKPT(I,NK+KD,M)-QKPT(I,NK+1,M))/DBLE(KD-1)
                  end do
                  NK=NK+1
                  WGHTK(NK,M)=ONE
c
c              Note that the loop will not execute if kd <=1, so
c               no problem ensues:
c
c$$$               IF (KD.LE.1) GOTO 49
                  DO KK=2,KD
                     NK=NK+1
                     WGHTK(NK,M)=ONE
                     DO I=1,3
                        QKPT(I,NK,M)=QKPT(I,NK-1,M)+EL(I)
                     end do
                  end do
ctemp
                  write(*,*) nk,' K-points after direction ',kp
cend temp
               end do
               WGHTS=NK
               NKV(M)=NK
               IF (NKV(M).GT.NKD) STOP
     $              'PTIOE2: NKD DIMENSION EXCEEDED'
               WRITE (6,1010) NKV(M)
               WRITE(15,1010) NKV(M)
C
            ELSE
C
               IF (NKV(M).GT.NKD) STOP
     $              'PTIOE2: NKD DIMENSION EXCEEDED'
               WRITE (6,1010) NKV(M)
               WRITE(15,1010) NKV(M)
               IF (NKV(M).NE.999) THEN
                  NK=NKV(M)
                  DO I=1,NK
                     READ (kio,*) (QKPT(J,I,M),J=1,3),WGHTK(I,M)
                     DO J=1,3
                        QKPT(J,I,M)=QKPT(J,I,M)/SFK
                     end do
                     WGHTS=WGHTS+WGHTK(I,M)
                  end do
                  WGHTS=ONE/WGHTS
               END IF
            END IF
         ELSE
c$$$      LKPTGN=.TRUE.
            READ (kio,1013) LSPEC(M),LINFOK
            READ (kio,1014)  KPDIV1(M),KPDIV2(M),KPDIV3(M),KPST,KPEND
            WRITE ( 6,1015) LSPEC(M)
            WRITE ( 6,1016) KPDIV1(M),KPDIV2(M),KPDIV3(M),KPST,KPEND
            WRITE(15,1015) LSPEC(M)
            WRITE(15,1016) KPDIV1(M),KPDIV2(M),KPDIV3(M),KPST,KPEND
C
            IF (LINFOK) NKINFO=-NKV(M)
C
C     NOTE THAT NKV(M) WILL BE OVERWRITTEN IN KPTGEN, SO THAT
C     IF LINFOK.EQ.TRUE THEN NKINFO IS THE NUMBER OF INFORMATION-ONLY
C     K-POINTS (I.E. THEIR WEIGHTS ARE SET TO ZERO).
C     NOW GENERATE THE K-POINT MESH:
C
c$$$      CALL KPTGEN (25,IOK,M)
            CALL KPTGEN (M)
            NK=NKV(M)
C
C     READ IN INFORMATION-ONLY POINTS IF ANY:
C
            IF (LINFOK) THEN
               DO I=1,NKINFO
                  NK=NK+1
                  IF (NK.GT.NKD)
     $                 STOP 'PTIOE2: INFO PTS. EXCEED NKD DIMENSION'
                  READ (kio,*) (QKPT(J,NK,M),J=1,3)
                  WGHTK(NK,M)=ZERO
                  DO J=1,3
                     QKPT(J,NK,M)=QKPT(J,NK,M)/SFK
                  end do
               end do
               NKV(M)=NK
            END IF
C
            DO I=1,NK
               WGHTS=WGHTS+WGHTK(I,M)
            end do
C
            IF (ABS(WGHTS-ONE).GT.EPS) THEN
               WRITE ( 6,*) 'PTIOE2: SUM OF WEIGHTS .NE. 1, SUM=',WGHTS
               WRITE(15,*) 'PTIOE2: SUM OF WEIGHTS .NE. 1, SUM=',WGHTS
               wghtsi = one/wghts
               do i=1,nk
                  wghtk(i,m)=wghtk(i,m)*wghtsi
               end do
            END IF
C
            IF (KPST.LE.0) KPST=1
            IF (KPEND.EQ.0) KPEND=NKV(M)
C
            IF ((KPST.GT.1) .OR. (KPEND.LT.NKV(M))) THEN
               WGHTS=ZERO
               NK=0
               DO I=KPST,KPEND
                  NK=NK+1
                  WGHTK(NK,M)=WGHTK(I,M)
                  WGHTS=WGHTS+WGHTK(NK,M)
                  DO J=1,3
                     QKPT(J,NK,M)=QKPT(J,I,M)
                  end do
               end do
               NKV(M)=NK
               WGHTS=ONE/WGHTS
            END IF
C
            WRITE ( 6,1010) NKV(M)
            WRITE(15,1010) NKV(M)
C
         END IF
C
         write(6,'(i5)') nk
         DO I=1,NK
            WGHTK(I,M)=WGHTK(I,M)*WGHTS
            CALL MAPK(EL,QKPT(1,I,M),-1)
            WRITE(6,'(4f19.15)') (QKPT(J,I,M),J=1,3),WGHTK(I,M)
            WRITE(15,1011) (QKPT(J,I,M),J=1,3),WGHTK(I,M),
     $        EL(1),EL(2),EL(3)
         end do
C NEWSYM
      endif
ctemp
c$$$      write(*,*) 'Assigning k-points'
cend temp
      do i=1,nk
         x(i) = qkpt(1,i,m)
         y(i) = qkpt(2,i,m)
         z(i) = qkpt(3,i,m)
      end do
      NPTS=NK
c      IF(.NOT.NEWSYM)THEN
C IF NOT NEWSYM
c      do 120 i=1,NATOMS
c        do 120 J=1,3
c         POSN(I,J)=POS(J,I)
c 120  continue
c      endif
c
c     Close k-point file, if we must
c
      if(kio.eq.33) close(kio)
ctemp
c$$$      write(*,*) 'Exiting kptin'
cend temp
      RETURN
c
c -- Begin mjm --
c
13000 write(0,13005) kfile
13005 format(/'Cannot find k-point file'/a255)
c$$$      call exit(3)
      stop 'kptin:  Exit 3'
c
c -- End mjm --
c
 1000 FORMAT(I5,2F10.6)
 1001 FORMAT (/10X,'NUMBER OF    ENERGY WINDOWS=',I3/
     *         10X,'NUMBER OF SEMI-CORE WINDOWS=',I3/
     +        10X,'NUMBER OF VALENCE ELECTRONS=',F10.6/
     +        10X,'TKB=',F10.6)
 1002 FORMAT(5F10.6)
 1003 FORMAT(I4,4F10.6)
CLPS93
 1005 FORMAT(10X,I2,15F8.2)
CLPS93
 1006 FORMAT(12X,10F8.2)
 1007 FORMAT(/10X,'VACUUM:',2F8.2)
 1008 FORMAT(10X,'RKM = ',F8.3)
 1010 FORMAT(/,10X,'NUMBER OF K-POINTS=',I5,/,
     *       10X,'BVEC CORD',14X,'WT',14X,'CARTESIAN CORD',/)
 1011 FORMAT(' (',3F10.5,' )',F13.8,'  (',3F10.5,' )')
 1012 FORMAT (2I5,2F10.6)
 1013 FORMAT (6X,L1,8X,L1)
 1014 FORMAT (5I5)
 1015 FORMAT (//10X,' THE K-POINT GENERATOR IS INVOKED:  LSPEC=',L1//
     *          10X,'   .T. MEANS SPECIAL K-POINTS ARE GENERATED'/
     *          10X,'   .F. MEANS UNIFORM K-POINT MESH GENERATED')
 1016 FORMAT (//10X,'BRILLOUIN ZONE DIVISIONS:',/
     *          10X,'      KDIV1,KDIV2,KDIV3=',3I5/,
     *          10X,'      KPST,KPEND=',2I5,/)
 1017 FORMAT(I5,6F10.5)
C
      END
