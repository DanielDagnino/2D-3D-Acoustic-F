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
      SUBROUTINE DSINT1B ( N, INC, X, LENX, WSAVE, LENSAV, 
     1                   WORK, LENWRK, IER)
      INTEGER    N, INC, LENX, LENSAV, LENWRK, IER
      REAL       X(INC,*), WSAVE(LENSAV), WORK(LENWRK)
C
      IER = 0
C
      IF (LENX .LT. INC*(N-1) + 1) THEN
        IER = 1
        CALL DXERFFT ('SINT1B', 6)
        GO TO 100
      ELSEIF (LENSAV .LT. N/2 + N + INT(LOG(REAL(N))/LOG(2.)) +4) THEN
        IER = 2
        CALL DXERFFT ('SINT1B', 8)
        GO TO 100
      ELSEIF (LENWRK .LT. (2*N+2)) THEN
        IER = 3
        CALL DXERFFT ('SINT1B', 10)
        GO TO 100
      ENDIF
C
      CALL DSINTB1(N,INC,X,WSAVE,WORK,WORK(N+2),IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL DXERFFT ('SINT1B',-5)
      ENDIF
C
  100 CONTINUE
      RETURN
      END
