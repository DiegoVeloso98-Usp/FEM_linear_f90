* *******************************************************************
* COPYRIGHT (c) 1993 Council for the Central Laboratory
*                    of the Research Councils
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of the ACADEMIC
* SOFTWARE LICENCE, see http://hsl.rl.ac.uk/hsl2007/cou/academic.html
*
* Please note that for an ACADEMIC Licence:
*
* 1. The Packages may only be used for academic research or teaching
*    purposes by the Licensee, and must not be copied by the Licensee for
*    use by any other persons. Use of the Packages in any commercial
*    application shall be subject to prior written agreement between
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory (STFC) of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither STFC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
* Original date 15 March 1993

* 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC30AD(N,NE,A,IRN,ICN,S,W,LP,IFAIL)
      INTEGER N,NE
      DOUBLE PRECISION A(NE)
      INTEGER IRN(NE),ICN(NE)
      DOUBLE PRECISION S(N),W(N,4)
      INTEGER LP,IFAIL
      INTRINSIC LOG,ABS,MAX,MIN
      INTEGER M,MAXIT,R,P,MP
      PARAMETER (M=1,MAXIT=10,MP=4,P=3,R=2)
      DOUBLE PRECISION ONE,RMIN,ZERO
      PARAMETER (ONE=1D0,RMIN=0.1,ZERO=0D0)
      DOUBLE PRECISION AK,BK
      INTEGER I,ITER,J,K
      DOUBLE PRECISION PP,RM,RR,RRL,U
      IFAIL = 0
      IF (N.LT.1) THEN
         IFAIL = -1
         GO TO 130
      ELSE IF (NE.LE.0) THEN
         IFAIL = -2
         GO TO 130
      END IF
      DO 10 I = 1,N
         S(I) = ZERO
         W(I,M) = ZERO
         W(I,R) = ZERO
   10 CONTINUE
      DO 40 K = 1,NE
         U = ABS(A(K))
         IF (U.EQ.ZERO) GO TO 40
         I = IRN(K)
         J = ICN(K)
         IF ( MIN(I,J).LT.1 .OR. MAX(I,J).GT.N ) GO TO 40
         U = LOG(U)
         W(I,M) = W(I,M) + ONE
         W(I,R) = W(I,R) - U
         W(J,M) = W(J,M) + ONE
         IF(I.EQ.J) GO TO 40
         W(J,R) = W(J,R) - U
   40 CONTINUE
      RR = ZERO
      DO 50 I = 1,N
         IF (W(I,M).EQ.ZERO) W(I,M) = ONE
         W(I,P) = W(I,R)/W(I,M)
         W(I,MP) = W(I,R)
         RR = RR + W(I,R)**2/W(I,M)
   50 CONTINUE
      RM = RMIN*NE
      IF (RR.LE.RM) RETURN
      DO 120 ITER = 1,MAXIT
         DO 80 K = 1,NE
            IF (A(K).EQ.ZERO) GO TO 80
            I = ICN(K)
            J = IRN(K)
            IF(I.EQ.J) GO TO 80
            IF ( MIN(I,J).LT.1 .OR. MAX(I,J).GT.N ) GO TO 80
            W(I,MP) = W(I,MP) + W(J,P)
            W(J,MP) = W(J,MP) + W(I,P)
   80    CONTINUE
         PP = ZERO
         DO 90 I = 1,N
            PP = PP + W(I,P)*W(I,MP)
   90    CONTINUE
         AK = RR/PP
         RRL = RR
         RR = ZERO
         DO 100 I = 1,N
            S(I) = S(I) + AK*W(I,P)
            W(I,R) = W(I,R) - AK*W(I,MP)
            RR = RR + W(I,R)**2/W(I,M)
  100    CONTINUE
         IF (RR.LE.RM) RETURN
         BK = RR/RRL
         DO 110 I = 1,N
            W(I,P) = W(I,R)/W(I,M) + BK*W(I,P)
            W(I,MP) = W(I,P)*W(I,M)
  110    CONTINUE
  120 CONTINUE
  130 IF (LP.GT.0) WRITE (LP,'(/A/A,I3)')
     +    ' **** Error return from MC30AD ****',' IFAIL =',IFAIL
      END
