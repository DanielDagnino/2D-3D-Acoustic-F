C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2011 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                     FFTPACK  version 5.1                      *
C     *                                                               *
C     *                 A Fortran Package of Fast Fourier             *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *               Paul Swarztrauber and Dick Valent               *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE DRFFT2F (LDIM, L, M, R, WSAVE, LENSAV, WORK,
     1  LENWRK, IER)
      INTEGER LDIM, L, M, LENSAV, LENWRK, IER, IDX, MODL, MODM,
     1                      IDH, IDW
      REAL    R(LDIM,M), WSAVE(LENSAV), WORK(LENWRK)
C
C
C INITIALIZE IER
C
      IER = 0
C
C VERIFY LENSAV
C
      LWSAV =   L+INT(LOG(REAL(L))/LOG(2.))+4
      MWSAV =   2*M+INT(LOG(REAL(M))/LOG(2.))+4
      MMSAV =   M+INT(LOG(REAL(M))/LOG(2.))+4
C
      IF (LENSAV .LT. LWSAV+MWSAV+MMSAV) THEN
        IER = 2
        CALL DXERFFT ('RFFT2F', 6)
        GO TO 100
      ENDIF
C
C VERIFY LENWRK
C
      IF (LENWRK .LT. (L+1)*M) THEN
        IER = 3
        CALL DXERFFT ('RFFT2F', 8)
        GO TO 100
      ENDIF
C
C VERIFY LDIM IS AS BIG AS L
C
      IF (LDIM .LT. L) THEN
        IER = 5
        CALL DXERFFT ('RFFT2F', -6)
        GO TO 100
      ENDIF
C
C TRANSFORM FIRST DIMENSION OF ARRAY
C
      CALL DRFFTMF(M,LDIM,L,1,R,M*LDIM,WSAVE(1),
     .     L+INT(LOG(REAL(L))/LOG(2.))+4,WORK,LENWRK,IER1)
C
      IF(IER1.NE.0) THEN
         IER=20
         CALL DXERFFT('RFFT2F',-5)
         GO TO 100
      ENDIF
C
      LDX = 2*INT((L+1)/2)-1
      DO I=2,LDX
      DO J=1,M
      R(I,J) = .5*R(I,J)
      END DO
      END DO
      DO J=1,M
      DO I=3,LDX,2
      R(I,J) = -R(I,J)
      END DO
      END DO
C
C     PRINT*, 'FORWARD TRANSFORM IN THE I DIRECTION'
C     DO I=1,L
C       PRINT*, (R(I,J),J=1,M)
C     END DO
C
C RESHUFFLE TO ADD IN NYQUIST IMAGINARY COMPONENTS
C
      MODL = MOD(L,2)
      MODM = MOD(M,2)
C
C TRANSFORM SECOND DIMENSION OF ARRAY
C
      CALL DRFFTMF(1,1,M,LDIM,R,M*LDIM,
     1     WSAVE(LWSAV+MWSAV+1),MMSAV,WORK,LENWRK,IER1)
      DO J=2,2*((M+1)/2)-1
      R(1,J) = .5*R(1,J)
      END DO
      DO J=3,M,2
      R(1,J) = -R(1,J)
      END DO
      LDH = INT((L+1)/2)
      IF(LDH.GT.1) THEN
      LDW = LDH+LDH
C
C     R AND WORK ARE SWITCHED BECAUSE THE THE FIRST DIMENSION
C     OF THE INPUT TO COMPLEX CFFTMF MUST BE EVEN.
C
      CALL DR2W(LDIM,LDW,L,M,R,WORK)
      CALL DCFFTMF(LDH-1,1,M,LDH,WORK(2),LDH*M,
     1     WSAVE(LWSAV+1),MWSAV,R,L*M, IER1)
      IF(IER1.NE.0) THEN
         IER=20
         CALL DXERFFT('RFFT2F',-5)
         GO TO 100
      ENDIF
      CALL DW2R(LDIM,LDW,L,M,R,WORK)
      END IF
C
      IF(MODL.EQ.0) THEN
      CALL DRFFTMF(1,1,M,LDIM,R(L,1),M*LDIM,
     1     WSAVE(LWSAV+MWSAV+1),MMSAV,WORK,LENWRK,IER1)
      DO J=2,2*((M+1)/2)-1
      R(L,J) = .5*R(L,J)
      END DO
      DO J=3,M,2
      R(L,J) = -R(L,J)
      END DO
      END IF
C
C     PRINT*, 'FORWARD TRANSFORM IN THE J DIRECTION'
C     DO I=1,L
C       PRINT*, (R(I,J),J=1,M)
C     END DO
C
      IF(IER1.NE.0) THEN
         IER=20
         CALL DXERFFT('RFFT2F',-5)
         GO TO 100
      ENDIF
C
C
  100 CONTINUE
      RETURN
      END
