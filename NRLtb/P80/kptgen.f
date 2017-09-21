      SUBROUTINE KPTGEN  (NW)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C           GENERATES THE K-PT MESH IN THE IRREDUCIBLE BRILLOUIN ZONE
C           AND PROVIDES A MAPPING TO EQUIVALENT POINTS IN THE FULL BZ.
C
C           IF LSPEC=.TRUE.  SPECIAL K-POINTS ARE GENERATED
C           IF LSPEC=.FALSE.  CONVENTIONAL UNIFORM MESH
C                             K-POINTS ARE GENERATED
C
C                      BY H. KRAKAUER JULY 1984
C
C                      FIXED BY D. SINGH 1987
C
c=======================================================================
c
c     REVISION HISTORY:
c
c-----------------------------------------------------------------------
c
c     Version 1.06:
c
c     Moved k-point generation parameters into 'P3'
c
c                                                mjm --  6 July  1998
c=======================================================================
c
      LOGICAL LSPEC
      logical LBCC,LHEX,LTET
      include 'P1'
      include 'P3'
c$$$      PARAMETER (NKBZD=96*NKD)
c$$$      parameter (nkbzd = nopd*nkd)
      parameter (nkbzd = 2*nkd)
C
      COMMON /LATTYP/ AVEC(3,3),BVEC(3,3),BIJ(3,3),WSVOL,LATTIC
      COMMON /SPACE / RSPC(3,3,NOPD),RSPCT(3,3,NOPD),TNPSPC(3,NOPD),NOP
      common /kmesh1/ qkpt(3,nkd,nwdd),wghtk(nkd,nwdd),nkv(nwdd)
      COMMON /KMESH2/ KPDIV1(NWDD),KPDIV2(NWDD),KPDIV3(NWDD),
     $     LSPEC(NWDD)
      real*8 RECIPK(3,NKBZD),REALK(3,NKBZD),RADK(NKBZD),TEM1(NKBZD)
      integer IST(NKBZD+1),INDEX(NKBZD)
      logical LIN(NKBZD)
      DIMENSION U(3),UU(3),UK(3)
c
      parameter (toler = 1d-4)
C
      parameter (pi = 3.14159265358979323846d0)
      parameter (fpi = 4d0*pi)
c
c      IF(NW.EQ.1)THEN
c      OPEN(90,FILE='AFOP',STATUS='UNKNOWN')
c      READ(90,'(3F20.0)',END=400)((AFCART(I,J),I=1,3),J=1,3)
c      READ(90,'(3F20.0)',END=400)AFSH(1),AFSH(2),AFSH(3)
c      DO 402 I=1,3
c      DO 402 J=1,3
c      T(I,J)=BVEC(J,I)/TPI
c  402 TINV(I,J)=AVEC(I,J)
c      CALL MATMLT(AFCART,TINV,WORK,3)
c      CALL MATMLT(T,WORK,AFDR,3)
c      DO 403 I=1,3
c      DO 403 J=1,3
c      T(I,J)=AVEC(J,I)/TPI
c  403 TINV(I,J)=BVEC(I,J)
c      CALL MATMLT(AFCART,TINV,WORK,3)
c      CALL MATMLT(T,WORK,AFREC,3)
c  400 CONTINUE
c      CLOSE(90)
c      ENDIF
C
      LBCC=LATTIC.EQ.2
      LHEX=LATTIC.EQ.4
      LTET=(LATTIC.EQ.7).OR.(LATTIC.EQ.8).OR.(LATTIC.EQ.9)
      IF(LTET)THEN
         WRITE( 6,*)' SPECIAL TREATMENT OF BCT LATTICE'
         WRITE( *,*)' SPECIAL TREATMENT OF BCT LATTICE'
      ENDIF
C
      IF (LBCC) THEN
C
C         SPECIAL TREATMENT FOR BCC LATTICE USING SIMPLE
C         CUBIC CELL:
C
          WRITE ( 6,*) '    BCC LATTICE IN KPTGEN'
          WRITE ( *,*) '    BCC LATTICE IN KPTGEN'
          ALAT=(2.D0*WSVOL)**0.33333333333333333D0
          WRITE ( 6,*) 'ALAT= ',ALAT
          WRITE ( *,*) 'ALAT= ',ALAT
          ALAT=FPI/ALAT
      ENDIF
C
      KDIV1=KPDIV1(NW)
      KDIV2=KPDIV2(NW)
      KDIV3=KPDIV3(NW)
      IF(LTET)THEN
         KDIV1=2*KDIV1
         KDIV2=2*KDIV2
      ENDIF
      IF(LSPEC(NW))THEN
         IF(LHEX)THEN
            MAX1=KDIV1+MOD(KDIV1,2)
            MAX2=MAX1
            IF (KDIV1.NE.KDIV2) STOP 'KPTGEN: KDIV ERROR FOR HEX LAT.'
         ELSE
            MAX1=KDIV1+1
            MAX2=KDIV2+1
         ENDIF
         MAX3=KDIV3+1
      ELSE
         MAX1=KDIV1
         MAX2=KDIV2
         MAX3=KDIV3
      ENDIF
      NK=0
      KIN=0
      F1=0.5D0/DBLE(KDIV1)
      F2=0.5D0/DBLE(KDIV2)
      F3=0.5D0/DBLE(KDIV3)
      DO II=-MAX1,MAX1,2
         UK(1)=F1*DBLE(II)
         DO JJ=-MAX2,MAX2,2
            UK(2)=F2*DBLE(JJ)
            DO KK=-MAX3,MAX3,2
               UK(3)=F3*DBLE(KK)
               IF(LATTIC.EQ.7)THEN
                  UK(1)=DBLE(II)/DBLE(KDIV1)
                  UK(2)=DBLE(JJ)/DBLE(KDIV2)
                  UK(3)=UK(3)+0.5D0*(UK(1)-UK(2))
               ENDIF
               IF(LATTIC.EQ.8)THEN
                  UK(1)=DBLE(II)/DBLE(KDIV1)
                  UK(2)=DBLE(JJ)/DBLE(KDIV2)
                  UK(3)=UK(3)+0.5D0*UK(1)
               ENDIF
               IF(LATTIC.EQ.9)THEN
                  UK(1)=DBLE(II)/DBLE(KDIV1)
                  UK(2)=DBLE(JJ)/DBLE(KDIV2)
                  UK(3)=UK(3)+0.5D0*(UK(1)+UK(2))
               ENDIF
               NK=NK+1
               IF(NK.GT.NKBZD)THEN
                  WRITE( 6,*)'KPTGEN: NKBZD DIM. EXCEED., NK= ',NK
                  WRITE( *,*)'KPTGEN: NKBZD DIM. EXCEED., NK= ',NK
                  STOP'KPTGEN: NKBZD DIMENSION EXCEEDED'
               ENDIF
               U(1)=UK(1)
               U(2)=UK(2)
               U(3)=UK(3)
               LIN(NK)=II.LE.KDIV1.AND.II.GT.-KDIV1
     *              .AND.JJ.LE.KDIV2.AND.JJ.GT.-KDIV2
     *              .AND.KK.LE.KDIV3.AND.KK.GT.-KDIV3
               IF(LIN(NK))KIN=KIN+1
C
               IF(LBCC)THEN
                  DO J=1,3
                     UU(J)=UK(J)*ALAT
                     REALK(J,NK)=UU(J)
                  end do
                  CALL MAPK (UU,U,1)
               ELSE
                  CALL MAPK (REALK(1,NK),UK,-1)
               ENDIF
C
C     SHIFT TO THE WIGNER SEITZ CELL BY
C     FINDING THE SMALLEST MAGNITUDE EQUIVALENT K-VECTOR
C     BY SHIFTING WITH SOME NEAR SHELL G-VECTORS:
C
               R2=1.D50
               DO I1=-3,3
                  DO I2=-3,3
                     DO I3=-3,3
                        UU(1)=U(1)+I1
                        UU(2)=U(2)+I2
                        UU(3)=U(3)+I3
                        R=GDOT(UU,UU,BIJ)
                        IF(R.LE.R2)THEN
                           R2=R
                           UT1=UU(1)
                           UT2=UU(2)
                           UT3=UU(3)
                        ENDIF
                     end do
                  end do
               end do
               RECIPK(1,NK)=UT1
               RECIPK(2,NK)=UT2
               RECIPK(3,NK)=UT3
               RADK(NK)=R2
            end do
         end do
      end do
C
      CALL HSORT(NK,RADK,INDEX)
      F=1.D0/DBLE(KIN)
      DO J=1,NK
         TEM1(J)=RADK(INDEX(J))
      end do
      DO J=1,NK
         RADK(J)=TEM1(J)
      end do
      DO I=1,3
         DO J=1,NK
            TEM1(J)=RECIPK(I,INDEX(J))
         end do
         DO J=1,NK
            RECIPK(I,J)=TEM1(J)
         end do
      end do
      R2=RADK(1)+TOLER
      NST=1
      IST(1)=1
      DO J=2,NK
         IF(RADK(J).GT.R2)THEN
            R2=RADK(J)+TOLER
            NST=NST+1
            IST(NST)=J
         ENDIF
      end do
      IST(NST+1)=NK+1
      NOVK=0
      DO I=1,NST
         K1=NOVK+1
         IS1=IST(I)
         IS2=IST(I+1)-1
         DO J=IS1,IS2
            IEQ=0
            XX=RECIPK(1,J)
            YY=RECIPK(2,J)
            ZZ=RECIPK(3,J)
            DO K=K1,NOVK
               U(1)=QKPT(1,K,NW)
               U(2)=QKPT(2,K,NW)
               U(3)=QKPT(3,K,NW)
               DO IOX=1,NOP
                  DO I1=-1,1
                     X=RSPCT(1,1,IOX)*U(1)+RSPCT(1,2,IOX)*U(2)+
     *                    RSPCT(1,3,IOX)*U(3)+DBLE(I1)
                     DO I2=-1,1
                        Y=RSPCT(2,1,IOX)*U(1)+RSPCT(2,2,IOX)*U(2)+
     *                       RSPCT(2,3,IOX)*U(3)+DBLE(I2)
                        DO I3=-1,1
                           Z=RSPCT(3,1,IOX)*U(1)+RSPCT(3,2,IOX)*U(2)+
     *                          RSPCT(3,3,IOX)*U(3)+DBLE(I3)
                           IF(ABS(X-XX)+ABS(Y-YY)+ABS(Z-ZZ).LE.1.D-4)
     $                          IEQ=K
                           IF(IEQ.EQ.0.AND.(ABS(X+XX)+ABS(Y+YY)+
     $                          ABS(Z+ZZ)).LE.1.D-4) IEQ=K
                        end do
                     end do
                  end do
               end do
            end do
            IF(IEQ.EQ.0)THEN
               NOVK=NOVK+1
               QKPT(1,NOVK,NW)=RECIPK(1,J)
               QKPT(2,NOVK,NW)=RECIPK(2,J)
               QKPT(3,NOVK,NW)=RECIPK(3,J)
               IF(LIN(INDEX(J)))THEN
                  WGHTK(NOVK,NW)=F
               ELSE
                  WGHTK(NOVK,NW)=0.D0
               ENDIF
            ELSE
               IF(LIN(INDEX(J)))WGHTK(IEQ,NW)=WGHTK(IEQ,NW)+F
            ENDIF
         end do
      end do
      NKV(NW)=NOVK
C
c     Set nwndw = 1
c
c$$$c$$$      IF(TRIA.AND.NW.EQ.NWNDW)THEN
c$$$      IF(TRIA.AND.(NW.EQ.1)) THEN
c$$$      OPEN (IOBZP,FILE='BZKPT',STATUS='UNKNOWN',FORM='UNFORMATTED')
c$$$C TEST
c$$$C     WRITE( *,'(I5)')NK
c$$$C     WRITE( *,'(3F12.7,I7)')((REALK(J,I),J=1,3),IRDUK(I),I=1,NK)
c$$$C
c$$$      REWIND IOBZP
c$$$      WRITE (IOBZP) NK
c$$$      WRITE (IOBZP) (IRDUK(I),I=1,NK)
c$$$      WRITE (IOBZP) ((REALK(J,I),J=1,3),I=1,NK)
c$$$      CLOSE (IOBZP)
c$$$      ENDIF
      RETURN
      END
