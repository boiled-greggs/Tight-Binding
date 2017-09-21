      SUBROUTINE RDSPGP(IPRT,ILAT,IO)
      IMPLICIT REAL*8 (A-H,O-Z)
C     ******************************************************************
C     *****  THIS SUBROUTINE READS THE SPACE GROUP OPERATIONS (R/TAU).
C     *****  THE NON-PRIMITIVE TRANSLATION VECTOR, TAU, IS READ IN
C     *****  DIRECT-LATTICE COORDINATES AND STORED IN ARRAY TNPSPC.
C     *****  THE ROTATION MATRIX, R, IS READ IN EITHER CARTESIAN
C     *****  COORDINATES OR DIRECT-LATTICE COORDINATES (A SIMILARITY
C     *****  TRANSFORMATION IS PERFORMED TO OBTAIN ALL REPRESENTATIONS
C     *****  OF THE OPERATOR, R, REGARDLESS OF THE FORM THAT IS INPUT):
C     *****
C     *****      ILAT .EQ. TRUE   THE MATRIX, R, IS READ IN DIR. LAT.
C     *****                       COORDINATES
C     *****      ILAT .EQ. FALSE  THE MATRIX, R, IS READ IN CARTESIAN
C     *****                       COORDINATES
C     *****      RSPC             STORES THE MATRIX, R, IN CART. COOR.
C     *****      RSPDIR           STORES R IN DIR. LAT. COOR.
C     *****      RSPCT            STORES R IN RECIP. LAT. COOR.
C     *****
C     *****  OPTION TO REDEFINE THE SPACE GROUP OPERATIONS RELATIVE
C     *****  TO A NEW ORIGIN (PARALLEL AXES, HOWEVER) FROM THAT FOR
C     *****  WHICH THEY ARE DEFINED ON INPUT:
C     *****    SHIFT(3)  IS THE VECTOR (IN LATTICE COORDINATES) FROM
C     *****    THE OLD ORIGIN TO THE NEW ORIGIN.  IF THIS VECTOR
C     *****    IS NOT ZERO ON INPUT, RDSPGP WILL CALCULATE THE NEW
C     *****    NON-PRIMITIVE TRANSLATION VECTORS.
C     *****
C     *****  NOTE:  LARGER INPUT FORMAT IS USED FOR INPUT OF ROTATION
C     *****         MATRICES IN CARTESIAN COORDINATES.
C     *****
C     *****  INPUT  1) GPNAME,NOP,SHIFT    (A10,I5,5X,3F10.5)
C     *****         2) RSPCT OR RSPC       (18F4.0) OR (3D20.10)
C     *****         3) TNPSPC              (9F8.5)
C     *****
C     *****                            H. KRAKAUER   FEB 1983
C     *****
C     *****  MODIFIED TO GENERATE ISYE AND EQUIVALENT ATOM POSITIONS
C     *****                          C. S. WANG.   JAN. 1985
C     ******************************************************************
c
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
c-----------------------------------------------------------------------
c
c     Added a line to check that the number of space group operations,
c      NOP, is no larger than the allowed dimension, NOPD
c
c                                                mjm -- 18 Dec   1998
c
c     Fixed bug wherein statements 100[45] were printed to standard
c      output twice, rather than once to s.o. and once to unit 15.
c                                                     mjm -- 16 Aug 2000
c=======================================================================
c
      LOGICAL IPRT,ILAT
      CHARACTER*10 GPNAME
      INCLUDE 'P1'
      INCLUDE 'P3'
C
      DIMENSION T(3,3),TINV(3,3),WORK(3,3),SHIFT(3),R(3)
      COMMON /LATTYP/ AVEC(3,3),BVEC(3,3),BIJ(3,3),WSVOL,LATTIC
      COMMON /SPACE / RSPC(3,3,NOPD),RSPCT(3,3,NOPD),TNPSPC(3,NOPD),NOP
      COMMON /SPACED/ RSPDIR(3,3,NOPD)
C
      parameter (zero  = 0d0)
      parameter (one   = 1d0)
      parameter (two   = 2d0)
      parameter (three = 3d0)
c
      parameter (pi = 3.14159265358979323846d0)
      parameter (tpi = 2d0*pi)
c
      parameter (tpiinv = one/tpi)
C
 1001 FORMAT (A10,I5,5X,3F10.5)
 1002 FORMAT (18F4.0)
 1003 FORMAT (9F8.5)
 1004 FORMAT (//'        ON INPUT:',2X,
     * ' THE ROTATION MATRICES ARE IN DIRECT-LATTICE COORDINATES'/)
 1005 FORMAT (//'        ON INPUT:',2X,
     *       ' THE ROTATION MATRICES ARE IN CARTESIAN COORDINATES'/)
 1006 FORMAT (3D20.10)
 1007 FORMAT (//30X,'THERE ARE',I3,1X,'SPACE GROUP OPERATIONS'/
     *          30X,'THE POINT GROUP IS ',A10/)
 1008 FORMAT (4X,'OPER. NO.',1X,'NON-PRIMITIVE TRANSLATION',10X,
     * 'CARTESIAN ROTATION MATRIX',8X,
     * 'RECIPROCAL-LATTICE-COOR. ROTATION MATRIX'//)
 1009 FORMAT (/5X,I5,3F10.4)
 1010 FORMAT (45X,3F10.5,10X,3F10.5/)
 1011 FORMAT (/'        NEW ORIGIN IS SHIFTED FROM THE OLD ORIGIN BY:'//
     *         '        R=(',3(F8.5,2X),')  IN DIRECT-LATTICE UNITS.'//)
 1012 FORMAT (10X,F10.8,41H OF THE UNIT CELL VOLUME IS INTERSTITIAL.  )
 1021 FORMAT(' ATOM NO=',I5,4X,'ISYE=',I5)

 1031 FORMAT(' TYPE=',I4,'  ATOM=',I4,'  POS=',3F10.5)
C
C     SET UP THE TRANSFORMATION MATRICES FOR THE SIMILARITY
C     TRANSFORMATIONS BETWEEN CARTESIAN AND DIRECT-LATTICE
C     REPRESENTATIONS OF THE ROTATION OPERATOR.
C
      DO I=1,3
         DO J=1,3
            T(I,J)   =BVEC(J,I)*tpiinv
            TINV(I,J)=AVEC(I,J)
         end do
      end do
C
      READ (IO,1001) GPNAME,NOP,SHIFT(1),SHIFT(2),SHIFT(3)
c
c     If NOP > NOPD then bad things will happen
c
      if(nop.gt.nopd) stop 'RDSPGP:  NOP > NOPD -- REDIMENSION IN P3'
C
      IF (ILAT) THEN
C
C        ROTATION OPERATORS, R, ARE READ IN DIR. LAT. COOR.:
C
         WRITE( *,1004)
         WRITE(15,1004)
         READ (IO,1002) (((RSPDIR(I,J,K),J=1,3),I=1,3),K=1,NOP)
         READ (IO,1003) ((TNPSPC(I,K),I=1,3),K=1,NOP)
C
C        OBTAIN THE ROTATION MATRICES IN CARTESIAN COORDINATES BY APPLYING
C         THE SIMILARITY TRANSFORMATION RSPC=TINV*RSPDIR*T
C
         DO K=1,NOP
            CALL MATMLT (RSPDIR(1,1,K),T,WORK,3)
            CALL MATMLT (TINV,WORK,RSPC(1,1,K),3)
         end do
C
      ELSE
C
C        ROTATION MATRICES ARE INPUT IN CARTESIAN COORDINATES:
C
         WRITE( *,1005)
         WRITE(15,1005)
         DO K=1,NOP
            READ (IO,1006) ((RSPC(I,J,K),J=1,3),I=1,3)
         end do
         READ (IO,1003) ((TNPSPC(I,K),I=1,3),K=1,NOP)
         IF (LATTIC.EQ.4) THEN
            CS3D2=SQRT(THREE)/TWO
            DO K=1,NOP
               DO I=1,3
                  DO J=1,3
                     IF (ABS(RSPC(I,J,K)-CS3D2).LT.0.001)
     $                    RSPC(I,J,K)=CS3D2
                     IF (ABS(RSPC(I,J,K)+CS3D2).LT.0.001)
     $                    RSPC(I,J,K)=-CS3D2
                  end do
               end do
            end do
         END IF
C
C        OBTAIN THE ROTATION MATRICES IN DIRECT-LATTICE COORDINATES
C        BY APPLYING THE SIMILARITY TRANSFORMATION RSPDIR=T*RSPC*TINV:
C
         DO K=1,NOP
            CALL MATMLT (RSPC(1,1,K),TINV,WORK,3)
            CALL MATMLT (T,WORK,RSPDIR(1,1,K),3)
         end do
C
      END IF
C
C     NOW OBTAIN THE ROTATION OPERATORS, R, IN RECIPROCAL LATTICE
C     COORDINATES.  FIRST SET UP THE SIMILARITY TRANSFORMATION
C     FROM CARTESIAN COORDINATES.
C
      DO 100 I=1,3
      DO 100 J=1,3
      T(I,J)   =AVEC(J,I)/TPI
  100 TINV(I,J)=BVEC(I,J)
C
C     PERFORM RSPCT=T*RSPC*TINV
C
      DO 105 K=1,NOP
      CALL MATMLT (RSPC(1,1,K),TINV,WORK,3)
  105 CALL MATMLT (T,WORK,RSPCT(1,1,K),3)
C
C      GENERATE THE EQUIVALENT ATOM POSITONS.
C temporarily disabled
c$$$      NEXT=1
c$$$      IA2=0
c$$$      DO 360 III=1,NTYPE
c$$$      IA1=IA2+1
c$$$      IA2=IA2+NEQ(III)
c$$$      NBEGIN=NEXT
c$$$      DO 350 IA=IA1,IA2
c$$$      DO 340 IOP=1,NOP
c$$$      DO 330 I=1,3
c$$$      AT2=TNPSPC(I,IOP)
c$$$      DO 320 J=1,3
c$$$      AT2=AT2+RSPDIR(I,J,IOP)*POS(J,IA)
c$$$  320 CONTINUE
c$$$      TEMP(I,NEXT)=MOD(AT2+5.D0,1.D0)
c$$$  330 IF(ABS(TEMP(I,NEXT)-1.D0).LT.SMALL)TEMP(I,NEXT)=0.D0
c$$$      IF (NEXT.EQ.NBEGIN) THEN
c$$$      NEXT=NEXT+1
c$$$      GOTO 340
c$$$      ENDIF
c$$$      DO 332 JA=NBEGIN,NEXT-1
c$$$      DO 331 I=1,3
c$$$      IF(ABS(TEMP(I,JA)-TEMP(I,NEXT)).GT.SMALL)GOTO 332
c$$$  331 CONTINUE
c$$$      GOTO 340
c$$$  332 CONTINUE
c$$$      WRITE (6,9998) IOP,NEXT
c$$$ 9998 FORMAT(' OP=',I5,' ATOM=',I5)
c$$$      WRITE (6,9999) (TNPSPC(I,IOP),I=1,3)
c$$$      DO 339 I=1,3
c$$$  339 WRITE (6,9999) (RSPDIR(J,I,IOP),J=1,3)
c$$$      WRITE (6,9999) (TEMP(I,NEXT),I=1,3)
c$$$ 9999 FORMAT(' TS=',3F10.5)
c$$$      NEXT=NEXT+1
c$$$  340 CONTINUE
c$$$  350 CONTINUE
c$$$      NEQ(III)=NEXT-NBEGIN
c$$$  360 CONTINUE
c$$$      NAT=NEXT-1
c$$$      IF(NAT.GT.NATD) THEN
c$$$      IA2=0
c$$$      DO 363 N=1,NTYPE
c$$$      IA1=IA2+1
c$$$      IA2=IA2+NEQ(N)
c$$$      DO 363 I=IA1,IA2
c$$$      WRITE(6,1031) N,I,(TEMP(J,I),J=1,3)
c$$$      WRITE(*,1031) N,I,(TEMP(J,I),J=1,3)
c$$$  363 CONTINUE
c$$$      PRINT *,' NAT  .LT. NATD'
c$$$      PRINT *,NAT,NATD
c$$$      STOP
c$$$      ENDIF
c$$$      DO 361 IA=1,NAT
c$$$      DO 361 I=1,3
c$$$  361 POS(I,IA)=TEMP(I,IA)
c$$$
c$$$
c$$$
c$$$C
c$$$C     CALCULATE THE INTERSTITIAL VOLUME
c$$$C
c$$$      VINT=ZERO
c$$$      DO 365 JATOM=1,NTYPE
c$$$  365 VINT=VINT+NEQ(JATOM)*VOL(JATOM)
c$$$      VINT=ONE-VINT
c$$$C
c$$$      WRITE (6,1012) VINT
c$$$
c
C end temporary
C
C     CHECK TO SEE IF THE ORIGIN IS TO BE SHIFTED:
C
      S=ZERO
      DO 90 J=1,3
   90 S=S+ABS(SHIFT(J))
C
      IF (S.GT.0.001)
     *THEN
C
C     SHIFT IS NOT EQUAL TO (0,0,0) -- FIND SPACE GROUP OPERATIONS
C     RELATIVE TO THE NEW ORIGIN:
C
      WRITE (6,1011) (SHIFT(J),J=1,3)
      WRITE( *,1011) (SHIFT(J),J=1,3)
C
      DO 110 K=1,NOP
C
C     ROTATE THE SHIFT-VECTOR WITH RSPDIR
C
      DO 120 I=1,3
      R(I)=ZERO
      DO 120 J=1,3
  120 R(I)=R(I)+RSPDIR(I,J,K)*SHIFT(J)
C
C     DEFINE NEW NON-PRIMITIVE VECTORS:
C
      DO 130 J=1,3
  130 TNPSPC(J,K)=R(J)+TNPSPC(J,K)-SHIFT(J)
  110 CONTINUE
C
C     OBTAIN THE ATOMIC POSITION VECTORS WITH RESPECT TO THE
C     NEW ORIGIN:
C
c$$$      ISITE=0
c$$$      DO 200 N=1,NTYPE
c$$$      NAT=NEQ(N)
c$$$      DO 200 NQ=1,NAT
c$$$      ISITE=ISITE+1
c$$$      DO 210 J=1,3
c$$$  210 POS(J,ISITE)=POS(J,ISITE)-SHIFT(J)
c$$$  200 CONTINUE
c$$$C
c$$$C     ALL SPACE GROUP OPERATIONS ARE NOW DEFINED RELATIVE TO THE
c$$$C     NEW COORDINATE SYSTEM.
c$$$C
          ENDIF
c$$$C
c$$$      IA2=0
c$$$      DO 380 III=1,NTYPE
c$$$      IA1=IA2+1
c$$$      IA2=IA2+NEQ(III)
c$$$      DO 370 IA=IA1,IA2
c$$$      WRITE (6,1031)III,IA,(POS(I,IA),I=1,3)
c$$$      WRITE (*,1031)III,IA,(POS(I,IA),I=1,3)
c$$$  370 CONTINUE
c$$$  380 CONTINUE

      IF (.NOT.IPRT) GOTO 209
C     OUTPUT RESULTS:
C
      WRITE (6,1007) NOP,GPNAME
      WRITE (6,1008)
      DO 30 K=1,NOP
      WRITE (6,1009) K,(TNPSPC(I,K),I=1,3)
      DO 31 I=1,3
   31 WRITE (6,1010) (RSPC(I,J,K),J=1,3),(RSPCT(I,J,K),J=1,3)
   30 CONTINUE
      WRITE ( *,1007) NOP,GPNAME
C
C     ADD MECHANICS TO FIND ISYE INTERNALLY
C
C
  209 CONTINUE
c$$$      IA2=0
c$$$      DO 260 III=1,NTYPE
c$$$      IA1=IA2+1
c$$$      IA2=IA2+NEQ(III)
c$$$      ISYE(IA1)=1
c$$$      WRITE (6,1021) IA1,ISYE(IA1)
c$$$      WRITE (*,1021) IA1,ISYE(IA1)
c$$$      IF (IA2.LE.IA1) GOTO 260
c$$$      DO 250 IA=IA1+1,IA2
c$$$      DO 240 IOP=1,NOP
c$$$      DO 230 I=1,3
c$$$      AT2=TNPSPC(I,IOP)
c$$$      DO 220 J=1,3
c$$$      AT2=AT2+RSPDIR(I,J,IOP)*POS(J,IA)
c$$$  220 CONTINUE
c$$$      AT2=MOD(AT2+TEN,ONE)
c$$$      AT1=MOD(POS(I,IA1)+TEN,ONE)
c$$$      IF(ABS(AT1-1.D0).LT.SMALL)AT1=0.D0
c$$$      IF(ABS(AT2-1.D0).LT.SMALL)AT2=0.D0
c$$$      IF(ABS(AT1-AT2).GT.SMALL)GOTO240
c$$$  230 CONTINUE
c$$$      ISYE(IA)=IOP
c$$$      WRITE (6,1021) IA,IOP
c$$$      WRITE (*,1021) IA,IOP
c$$$      GOTO 250
c$$$  240 CONTINUE
c$$$      WRITE(*,'(A80)')' ISYE CANNOT BE FOUND'
c$$$      STOP'ISYE NOT FOUND'
c$$$  250 CONTINUE
  260 CONTINUE
      RETURN
      END
