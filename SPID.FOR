!          PROGRAM SPID.FOR   FIFTH VERSION 18.02.96
!
! *** SPID= 'Simple Plastic Incremental Deformation'***
! *****************************************************
! *                                                   *  
! * METAL FORMING AND THE FINITE ELEMENT METHOD       *
! *                                                   * 
! *****************************************************
! 
!  THIS IS MAIN PROGRAM OF RIGID VISCOPLASTIC FINITE
!  ELEMENT METHOD FOR SIMPLE FORMING PROCESS.
      IMPLICIT INTEGER*4 (I-N),REAL*8 (A-H, O-Z)
      CHARACTER TITLE*70
      COMMON /TITL/ TITLE
      COMMON /RVA2/ EPS(5,100), STS(5,100), TEPS(100)
      COMMON /INOT/ INPT,MSSG,IUNIT,IUNI2,zz
      COMMON /MSTR/NUMNP,NUMEL,IPLNAX
      COMMON /INVR/ NOD(4,100),NBCD(2,100)
      COMMON /TSTP/ NINI,NCUR,NSEND,NITR,DTMAX
      COMMON /ITRC/ ITYP,ICONV
      INPT = 5
      MSSG = 6
      IUNIT= 3
      IUNI2= 4

!      
! READ INPUT
!
      CALL INPRED
      OPEN(IUNIT,FILE='SPID.OUT', STATUS='UNKNOWN',FORM='FORMATTED')
      OPEN(MSSG, FILE='SPID.MSG', STATUS='UNKNOWN',FORM='FORMATTED') 
!	OPEN(zz, FILE='zz.txt', STATUS='UNKNOWN',FORM='FORMATTED') 
      WRITE(MSSG,1020) TITLE
      CALL PRTINP
      CALL BAND(NOD,NUMEL,NUMNP) 
!
!   STEP SOLUTION
!
      NINI = NINI+1
!
      DO 300 N=NINI,NSEND
      NCUR = N
!
      WRITE(MSSG,1050) N
      IF(N .NE. NINI) GO TO 80
      ICOUNT = 0
   50 ITYP = 2
      CALL NONLIN
      ICOUNT =ICOUNT+1
   80 ITYP =1
      CALL NONLIN
      IF(ICONV .EQ. 2 .AND. ICOUNT .GT. 3) GO TO 900
      IF(ICONV .EQ. 2) GO TO 50
!
      CALL POTSOL
      CALL PRTSOL
      CALL RSTFIL
 300  CONTINUE
!
      CLOSE(IUNIT) 
      CLOSE(MSSG)
!	close(zz)
      STOP
!
 900  CONTINUE
      WRITE(MSSG,1070)
!
      STOP
 1020 FORMAT(1H1,//,5X,'OUTPUT OF S P I D',//,5X,' MESSAGE FILE FOR '/,5X,A,//)
 1050 FORMAT(///,' ITERATION PROCESS FOR STEP ',I5,//)
 1070 FORMAT(/,' STOP BECAUSE SOLUTION DOES NOT CONVERGE. ')
      END
! ***************************************************************
      SUBROUTINE ADDBAN(B,A,NQ,LM,QQ,PP)
!
!     ---- CHAPTER 7.2,
!     ---- EQUATION(7.8)
!
! ASSEMBLE GLOBAL STIFFNESS MATRIX FROM ELEMENTAL STIFFNESS MATRIX
!
      IMPLICIT REAL*8 (A-H, O-Z), INTEGER*4 (I-N)
      DIMENSION B(1), A(NQ,1), QQ(8), PP(8,8), LM(1)
!
      DO 100 I = 1,8
      II =LM(I)
      DO 50 J = 1,8
      JJ = LM(J)-LM(I)+1
      IF(JJ.LE. 0) GO TO 50
      A(II,JJ) =A(II,JJ)+PP(I,J)
   50 CONTINUE
      B(II) = B(II) + QQ(I)
  100 CONTINUE
      RETURN
      END
! *********************************************************************
      SUBROUTINE BAND(NOD,NUMEL,NUMNP)
!
! ---- CHAPTER 7.2
! DETERMINE MAXIMUM HALF BANDWIDTH, MBAND AND
!           TOTAL NUMBER OF EQUATIONS, NEQ 
!
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
      COMMON /CNEQ/ NEQ,MBAND
      DIMENSION NOD(4,1)
!
      MBAND = 0
      DO 100 N= 1, NUMEL
      NMIN = NOD(1,N)
      NMAX = NOD(1,N)
      DO 50 I= 2,4
      IF(NMIN .GT. NOD(I,N)) NMIN = NOD(I,N)
      IF(NMAX .LT. NOD(I,N)) NMAX = NOD(I,N)
   50 CONTINUE
      MB = (NMAX-NMIN+1)*2
      IF(MBAND.LT.MB) MBAND=MB
  100 CONTINUE
      NEQ = NUMNP*2
!	write(zz,*)mband,neq
      RETURN
      END
! ***********************************************************
      SUBROUTINE BANSOL(B,A,NQ,MM)
!
!     --------CHAPTER 7.2
! THIS SUBROUTINE SOLVES THE BANDED SYMMETRIC MATRIX EQUATIONS BY
! THE GAUSSIAN ELIMINATION
! B:  LOAD VECTOR
! A:  SYMMETRIC MATRIX IN BANDED FORM
! MM: HALF BANDWIDTH
! NQ: NUMBER OF EQUATIONS
!
      IMPLICIT REAL*8 (A-H, O-Z), INTEGER*4 (I-N)
      COMMON /INOT/INOT,MSSG,IUNIT,IUNI2
      DIMENSION B(1),A(NQ,1)
!
      DO 200 N= 1,NQ
      IF(A(N,1) .LE. 0.) GO TO 800
      DO 150 L =2,MM
      IF(A(N,L) .EQ. 0.) GO TO 150
      C = A(N,L) / A(N,1)
      I = N+L-1
      J=0
      DO 100 K = L,MM
      J = J+1
 100  A(I,J) = A(I,J)-C*A(N,K)
      A(N,L) = C
  150 CONTINUE
  200 CONTINUE
!
! LOAD VECTOR REDUCTION
!
      DO 300 N= 1,NQ
      DO 250 L= 2,MM
      I = N+L-1
      IF(I .GT. NQ) GO TO 250
      B(I)=B(I)-A(N,L)*B(N)
  250 CONTINUE
      B(N) = B(N) / A(N,1)
  300 CONTINUE
!
! BACK SUBSITUATION
!
      DO 400 M = 1,NQ
      N = NQ+1-M
      DO 350 K= 2,MM
      L =N+K-1
  350 B(N) = B(N)-A(N,K)*B(L)
  400 CONTINUE
      RETURN
!
  800 CONTINUE
      WRITE(MSSG,1020) N
      STOP
 1020 FORMAT(/,'NEGATIVE PIVOT AT EQUATION NO. ',I5)
      END
! *********************************************************
      SUBROUTINE DISBDY(URZ,NBCD,B,A,NEQ,MBAND,ITYP)
!
!   ----CHAPTER 7.3
! APPLY DISPLACEMENT BOUNDARY CONDITION
!
      IMPLICIT REAL*8 (A-H, O-Z), INTEGER*4 (I-N)
      DIMENSION B(1),A(NEQ,1),NBCD(1),URZ(1)
!
      IF(ITYP .EQ. 2) GO TO 120
      DO 100 N = 1,NEQ
      IF(NBCD(N) .EQ. 0) GO TO 100
      DO 70 I = 2, MBAND
      II = N-I+1
      IF(II .LE. 0) GO TO 50
      A(II,I)=0.
   50 CONTINUE
      II=N+I-1
      IF(II .GT. NEQ) GO TO 70
      A(N,I) =0.
   70 CONTINUE
      B(N) = 0.
      A(N,1) =1.
  100 CONTINUE
      RETURN
!
  120 CONTINUE
      DO 200 N = 1,NEQ
      IF(NBCD(N) .EQ. 0) GO TO 200
      DO 170 I = 2,MBAND
      II= N-I+1
      IF(II .LE. 0) GO TO 150 
      B(II) = B(II)-A(II,I)*URZ(N)
      A(II,I) = 0.
  150 CONTINUE
      II = N+I -1
      IF(II .GT. NEQ) GO TO 170
      B(II) = B(II) -A(N,I)*URZ(N)
      A(N,I) = 0.
  170 CONTINUE
      B(N) = URZ(N)
      A(N,1) =1.
  200 CONTINUE
      END
! ******************************************************
      SUBROUTINE ELSHLF(PP,QQ,RZ,URZ,EPS,TEPS,IPLNAX,IDREC,NEL) 
!     
!  -------CHAPTER 6.5
!  -------EQUATIONS (6.43), AND (6.44)
!  EVALUATION OF ELEMENTAL STIFFNESS MATRIX
!
!    IDREC: IF =1, NEWTON-RAPHSON ITERATION
!              =2, DIRECT ITERATION          
!
      IMPLICIT REAL*8 (A-H, O-Z), INTEGER*4 (I-N)
      DIMENSION RZ(2,1), URZ(2,1), B(4,8), EPS(1), TEPS(100)
      DIMENSION QQ(8), PP(8,8),S2(2),W2(2)
      DATA S2/-0.57735026918963D0,0.57735026918963D0/,W2/2*1.0D0/
!
      DO 10  I= 1, 8
      QQ(I) = 0.
      DO 10  J= 1, 8
      PP(I,J) = 0.
   10 CONTINUE
!
!  CARRY OUT ONE POINT INTEGRATION
!
      S = 0.
      T = 0.
      CALL STRMTX(RZ,B,WDXJ,S,T,IPLNAX,NEL)
      WDXJ = WDXJ*4
      CALL VSPLON(QQ,PP,B,URZ,EPS,WDXJ,IDREC)
!
!    REGULAR INTEGRATION
!
      DO 100 I = 1, 2
      S = S2(I)
      DO 50 J = 1, 2
      T= S2(J)
      CALL STRMTX(RZ,B,WDXJ,S,T,IPLNAX,NEL)
      WDXJ = WDXJ*W2(I)*W2(J)
      CALL VSPLST(QQ,PP,B,URZ,TEPS(1),WDXJ,IDREC)
   50 CONTINUE
  100 CONTINUE    
      RETURN
      END
! ******************************************************         
      SUBROUTINE FLWSTS(YS,FIP,STRAN,STRRT)
!
! USER SUPPLIED SUBROUTINE TO DESCRIBE THE MATERIAL
! FLOW STRESS.
! THIS SUBROUTINE SHOWS THE EXAMPLE OF
! YS = 10. *(STRAIN RATE)**0.1
!
      IMPLICIT INTEGER*4 (I,J,K,L,M,N), REAL*8 (A-H,O-Z)
      COMMON /RIGD/RTOL,ALPH,DIAT
!
      IF(STRRT .LT. ALPH) GO TO 100
      YS = 173.2*STRRT**0
      FIP= STRRT**(-0.9)
      RETURN
!
 100  YO = 173.2*ALPH**0
      FIP= YO/ALPH
      YS = FIP*STRRT
      RETURN
      END
! ******************************************************         
      SUBROUTINE FRCBDY(RZ,URZ,NBCD,TEPS,EPS,QQ,PP,IPLNAX)
!
!  APPLY FRICTION BOUNDARY CONDITION
!
      IMPLICIT INTEGER*4 (I-N), REAL*8 (A-H, O-Z)
      COMMON /DIES/FRCFAC
      COMMON /INOT/INPT,MSSG,IUNIT,IUNI2
      DIMENSION RZ(2,1),URZ(2,1),NBCD(2,1),EPS(5),QQ(8),PP(8,1),ER(2,2),FR(2),XY(2,2),VXY(2,2)
!
      DO 100 N= 1,4
      I1=N+1
      I2=N
      IF(N .EQ. 4) I1=1
      IF(NBCD(2,I1) .NE. 3 .OR. NBCD(2,I2) .NE. 3) GO TO 100
      CALL FLWSTS(FLOW,DUM,TEPS,EPS(5))
      XY(1,1)= RZ(1,I1)
      XY(2,1)= RZ(2,I1)
      XY(1,2)=RZ(1,I2)
      XY(2,2)=RZ(2,I2)
      VXY(1,1)=URZ(1,I1)
      VXY(2,1)=URZ(2,I1)
      VXY(1,2)=URZ(1,I2)
      VXY(2,2)=URZ(2,I2)
      CALL FRCINT(XY,VXY,FLOW,FR,ER,FRCFAC,IPLNAX)
      J1=I1*2-1
      J2=I2*2-1
      QQ(J1)= QQ(J1)+FR(1)
      QQ(J2)=QQ(J2)+FR(2)
      PP(J1,J1)=PP(J1,J1)+ER(1,1)
      PP(J2,J2)=PP(J2,J2)+ER(2,2)
      PP(J1,J2)=PP(J1,J2)+ER(1,2)
      PP(J2,J1)=PP(J2,J1)+ER(2,1)
  100 CONTINUE
      RETURN
      END
! *****************************************************
      SUBROUTINE FRCINT(RZ,URZ,FLOW,FR,ER,FRCFAC,IPLNAX)
!
!  -----CHAPTER 7.3
!  -----EQUATIONS(7.14)AND(7.15)  
!  INTEGRATION METHOD : SIMPSON'S FORMULA
!  THIS ROUTINE CALCULATES THE FRICTION MATRIX
!  USED FOR BOTH TYPES OF ITERATION SCHEME
!
      IMPLICIT INTEGER*4 (I-N), REAL*8 (A-H, O-Z)
      COMMON /INOT/ INPT,MSSG,IUNIT,IUIN2
      COMMON /ITRC/ ITYP,ICONV
      DIMENSION RZ(2,1),URZ(1),QQ(8),PP(8,1)
      DIMENSION SLIV(2),ER(2,2),FR(2)
      DATA PI/3.1415926535898D0/
      DATA UA/0.0005D0/
!
!   INITIALIZE FR AND ER ARRAY
!
      DO 10 I=1,2
      FR(I)=0.
      DO 10 J =1,2
      ER(I,J) =0.
   10 CONTINUE
      NINT = 5
      FAC = DSQRT((RZ(1,2)-RZ(1,1))**2+(RZ(2,2)-RZ(2,1))**2)
      FK  = FLOW*FRCFAC/SQRT(3.)
      DH  = 2./(NINT-1)
      S   =-1.-DH
      CON = 2./PI*FK
      WD  =DH/3.*FAC*0.5*CON
!
      DO 300 N =1,NINT
      S=S+DH
      H1=0.5*(1.-S)
      H2=0.5*(1.+S)
      WDXJ = WD
      IF(IPLNAX .NE. 1) GO TO 90
      RR = H1*RZ(1,1)+H2*RZ(1,2)
      WDXJ =RR*WDXJ
   90 CONTINUE
!
      IF(N .EQ. 1 .OR. N .EQ. NINT) GO TO 100 
      NMOD=N-N/2*2
      IF(NMOD .EQ. 0) WDXJ =WDXJ*4
      IF(NMOD .EQ. 1) WDXJ =WDXJ*2
  100 CONTINUE
!
      US = H1*URZ(1)+H2*URZ(3)
      AT = DATAN(US/UA)
      IF(ITYP .EQ. 2) GO TO 200
      US2= US*US
      USA= US2+UA*UA
      CT1=AT*WDXJ
      CT2=UA/USA*WDXJ
      GO TO 250
!
!  FOR D-ITERATION CASE
!
 200  CONTINUE
      IF(DABS(US) .LE. 1.0D-5) SLOP = UA/(UA*UA+US*US)
      IF(DABS(US) .GT. 1.0D-5) SLOP = AT/US
      CT1=0.
      CT2=SLOP*WDXJ
!
!  CALCULATE CONTRIBUTION TO STIFFNESS
!
 250  CONTINUE
      FR(1) = FR(1) -H1*CT1
      FR(2) = FR(2) -H2*CT1
      ER(1,1) = ER(1,1)+H1*H1*CT2
      ER(1,2) = ER(1,2)+H1*H2*CT2
      ER(2,2) = ER(2,2)+H2*H2*CT2
      ER(2,1) = ER(1,2)
  300 CONTINUE
      RETURN
      END
! ****************************************************
      SUBROUTINE INPRED
!
!  READ INPUT FROM INPUT FILE
!
      IMPLICIT REAL*8 (A-H, O-Z), INTEGER*4 (I-N)
      CHARACTER TITLE*70
      COMMON /TITL/TITLE
      COMMON /TSTP/NINI,NCUR,NSEND,NITR,DTMAX
      COMMON /RVA1/RZ(2,100),URZ(2,100),FRZ(2,100)
      COMMON /RVA2/EPS(5,100),STS(5,100),TEPS(100)
      COMMON /INVR/ NOD(4,100),NBCD(2,100)
      COMMON /DIES/FRCFAC
      COMMON /RIGD/RTOL,ALPH,DIAT
      COMMON /MSTR/NUMNP,NUMEL,IPLNAX
      COMMON /INOT/INPT,MSSG,IUNIT,IUNI2
!
!  READ MASTER CONTROL DATA
!
      OPEN(INPT,FILE='SPID.DAT',FORM='FORMATTED',STATUS='OLD')
      READ(INPT,1000) TITLE
      READ(INPT,*)NINI,NSEND,DTMAX
      READ(INPT,*)ALPH,DIAT
      READ(INPT,*)IPLNAX
!
!  READ DIE DATA
!
      READ(INPT,*) FRCFAC
!
!       READ FEM NODE INFORMATION
!
      READ(INPT,*) NUMNP
      IF(NUMNP .GT. 300) GO TO 500
      DO 20 I=1,NUMNP
      READ (INPT,*)N,(RZ(J,N),J=1,2)
   20 CONTINUE
!      
!    READ ELEMENT INFORMATION
!
      READ (INPT,*)NUMEL
      IF(NUMEL .GT. 300) GO TO 500
      DO 40 I=1, NUMEL
      READ (INPT,*)N, (NOD(J,N), J=1,4)
   40 CONTINUE
!
!  READ BOUNDARY CONDITION DATA
!
      DO 60 N = 1, NUMNP
      DO 60 I = 1,2
      NBCD(I,N) = 0
   60 CONTINUE
!      
      DO 80 N =1, NUMNP
      READ(INPT,*)M,(NBCD(I,M), I=1,2)
      IF(M .GE. NUMNP) GO TO 100
   80 CONTINUE
  100 CONTINUE    
!  
!   READ NODE VELOCITY DATA
!
      DO 120 N=1,NUMNP
      DO 120 I=1,2
      URZ(I,N) = 0.
  120 CONTINUE
      DO 140 N =1,NUMNP
      READ (INPT,*)M,(URZ(I,M),I=1,2)
      IF(M .GE. NUMNP)GO TO 160
  140 CONTINUE
  160 CONTINUE
  220 CONTINUE
!
!   READ STRAIN DATA
!
      IF(NINI .EQ. 0) GO TO 300
      DO 240 N=1,NUMEL
      READ(INPT,*)M,TEPS(M)               
  240 CONTINUE 
  300 CLOSE(INPT)
      RETURN
!
  500 CONTINUE
      WRITE(MSSG,1010) 
      STOP
!
 1000 FORMAT(A)
 1010 FORMAT(/,'SORRY, THIS PROGRAM CANNOT HADLE MORE THAN 100','NODES OR ELEMENTS')
      END
!  *******************************************************************
      SUBROUTINE NFORCE(QQ,FRZ,LM)
!
!   ADD NODAL POINT FORCE
!
      IMPLICIT REAL*8 (A-H, O-Z),INTEGER*4 (I-N)
      DIMENSION QQ(8),FRZ(1),LM(1)
!
      DO 100 I=1,8
      N = LM(I)
      FRZ(N) = FRZ(N)- QQ(I)
  100 CONTINUE
      RETURN 
      END
!  ******************************************************************* 
      SUBROUTINE NONLIN
!
!  THIS ROUTINE CONTROLS THE ITERATIONS
!
      IMPLICIT INTEGER*4 (I-N), REAL*8 (A-H, O-Z)
      COMMON /INOT/ INPT,MSSG,IUNIT,IUNI2
      COMMON /MSTR/NUMNP,NUMEL,IPLNAX
      COMMON /TSTP/NINI,NCUR,NSEND,NITR,DTMAX
      COMMON /ITRC/ITYP,ICONV
      COMMON /CNEQ/NEQ,MBAND
      COMMON /RVA1/RZ(2,100),URZ(2,100),FRZ(2,100)
      DIMENSION UNORM(2),ENORM(2),FNORM(2)
      COMMON  A(5000),B(200)
!
      RTOL = 0.00001
!	write(zz,*)ityp
      IF(ITYP .EQ. 2) RTOL = 0.0005
      ACOEF = 1.
      NSTEL = NEQ*MBAND
      IF(NSTEL .LE. 5000 .AND. NEQ.LE.200) GO TO 10
      WRITE(MSSG,1010)
      STOP
!
   10 CONTINUE
      DO 30 N = 1,2
      UNORM(N)=0.
      ENORM(N)=0.
      FNORM(N)=0.
   30 CONTINUE
!
      ITRMAX = 20
      IF(ITYP .EQ. 2)ITRMAX =200
      DO 200 N=1,ITRMAX
      NITR = N
      CALL STIFF(B,A,NEQ,MBAND,ITYP)
      IDREC = 1
      CALL NORM(FRZ,B,FDUM,DFN,NEQ,IDREC)
      IF(ITYP .EQ. 2) DFN =0.
      CALL BANSOL(B,A,NEQ,MBAND)
      IDREC = ITYP
      CALL NORM(URZ,B,UC,EC,NEQ,IDREC)
      IF(ITYP .EQ. 1) WRITE(MSSG,1030)N
      IF(ITYP .EQ. 2) WRITE(MSSG,1050)N
      WRITE(MSSG,1070) UC,EC,DFN
      IF(N .EQ. 1) GO TO 130
      IF(EC .LT. RTOL .AND. DFN .LT. RTOL) GO TO 300
      IF(ITYP .EQ. 2) GO TO 130
      IF(EC .LT. ENORM(2)) GO TO 100
!
!  ADJUST THE ACOEF
!
      ACOEF = ACOEF*0.7
      GO TO 130
  100 CONTINUE
      IF(ENORM(1) .GT. ENORM(2) .AND. ENORM(2) .GT. EC)ACOEF= ACOEF*1.3     
      IF(ACOEF .GT. 1.) ACOEF = 1.0
!
!   VELOCITY UPDATE
!
  130 CONTINUE
      NB = 0
      DO 150 I=1,NUMNP
      DO 150 J=1,2
      NB = NB+1
      IF(ITYP .EQ. 1) URZ(J,I) =URZ(J,I)+ACOEF*B(NB)
      IF(ITYP .EQ. 2) URZ(J,I)=B(NB)
  150 CONTINUE
!
  170 CONTINUE
      UNORM(1) = UNORM(2)
      ENORM(1) = ENORM(2)
      FNORM(1) = FNORM(2)
      UNORM(2) = UC
      ENORM(2) = EC
      FNORM(2) = DFN
  200 CONTINUE
!
!  SET FLAG
!
      ICONV = 2
      RETURN
!
  300 CONTINUE
! CONVERGED CASE
! SET FLAG
!
      ICONV =1
!
      RETURN
 1010 FORMAT(/,'YOU NEED MORE  SPACE IN THE BLANK COMMON')
 1030 FORMAT(/,'N-R ITERATION NO. ',I5,/)
 1050 FORMAT(/,'DRT ITERATION NO. ',I5,/)
 1070 FORMAT('VELOCITY NORM               = ',F15.7,/,'  REL.ERROR NORM       = ',F15.7,/,'  REL.FORCE ERROR NORM = ',F15.7,/)
      END
!  ******************************************************************* 
      SUBROUTINE NORM(URZ,V,UC,EROR,NEQ,ITYP)
!
!  CALCULATE THE EROR NORM FOR LINEAR AND NONLINEAR CASE
!
      IMPLICIT INTEGER*4 (I-N), REAL*8 (A-H, O-Z)
      DIMENSION URZ(1),V(1)
!
      UC = 0.
      EROR = 0.
      DO 100 N =1,NEQ
      UC = UC+URZ(N)*URZ(N)
      IF(ITYP .EQ. 1) EROR = EROR+V(N)*V(N)
      IF(ITYP .EQ. 2) EROR = EROR+(URZ(N)-V(N))**2
  100 CONTINUE
!
      UC=DSQRT(UC)
      EROR= DSQRT(EROR)
      IF(UC .NE. 0.) EROR= EROR / UC
      RETURN
      END
!  ******************************************************************* 
      SUBROUTINE POTSOL
!
! THIS SUBROUTINE HANDLES THE POST SOLUTION PROCEDURES, IE,
!    GEOMETRY UPDATES
!    STRESS EVALUATION
!    TOTAL STRAIN EVALUATION
!
      IMPLICIT REAL*8 (A-H, O-Z), INTEGER*4 (I-N)
      COMMON /TSTP/NINI,NCUR,NSEND,NITR,DTMAX
      COMMON /MSTR/NUMNP,NUMEL,IPLNAX
      COMMON /RIGD/RTOL,ALPH,DIAT
      COMMON /RVA1/RZ(2,100),URZ(2,100),FRZ(2,100)
      COMMON /RVA2/EPS(5,100),STS(5,100),TEPS(100)
      COMMON /DIES/FRCFAC
!
!  GEOMETRY UPDATES
!
      DO 100 N = 1,NUMNP
      RZ(1,N) = RZ(1,N) +DTMAX*URZ(1,N)
      RZ(2,N) = RZ(2,N) +DTMAX*URZ(2,N)  
  100 CONTINUE    
!  
! STRESS EVALUATION
!
      DO 200 N=1,NUMEL
      AL = EPS(5,N)
      IF(AL .LT. ALPH) AL =ALPH
      CALL FLWSTS(EFSTS,STRT,TEPS(N),AL)
      EM = (EPS(1,N)+EPS(2,N)+EPS(3,N))/3.
      DO 150 I = 1,3
      STS(I,N) =2./3.*EFSTS*(EPS(I,N)-EM) / AL +DIAT * EM *3.
 150  CONTINUE
      STS(4,N)= EFSTS* EPS(4,N) /AL/3.
      STS(5,N)= EFSTS
  200 CONTINUE
!
!  UPDATE TOTAL EFFECTIVE STRAIN
!
      DO 300 N=1,NUMEL
      TEPS(N)= TEPS(N)+EPS(5,N)*DTMAX
 300  CONTINUE
      RETURN
      END
! ************************************************************
      SUBROUTINE PRTINP
!
!  THIS SUBROUTINE PRINTS THE INPUT DATA
      IMPLICIT REAL*8 (A-H, O-Z), INTEGER*4 (I-N)
      CHARACTER TITLE*70
      COMMON /TITL/TITLE
      COMMON /TSTP/NINI,NCUR,NSEND,NITR,DTMAX
      COMMON / RVA1/RZ(2,100),URZ(2,100),FRZ(2,100)
      COMMON / RVA2/EPS(5,100),STS(5,100),TEPS(100)
      COMMON /INVR/NOD(4,100),NBCD(2,100)
      COMMON /DIES/FRCFAC
      COMMON /RIGD/RTOL,ALPH,DIAT
      COMMON /MSTR/NUMNP,NUMEL,IPLNAX
      COMMON /INOT/INPT,MSSG,IUNIT,IUNI2
!
!  INPUT SUMMARY
!
      WRITE(IUNIT,1010) TITLE
      WRITE(IUNIT,1020)
      WRITE(IUNIT,1030) NINI,NSEND,DTMAX
      WRITE(IUNIT,1050) ALPH,DIAT
      WRITE(IUNIT,1070) IPLNAX
      WRITE(IUNIT,1110) FRCFAC
      WRITE(IUNIT,1130) NUMNP
      WRITE(IUNIT,1150)
      WRITE(IUNIT,1180) (N,(RZ(I,N),I=1,2),N=1,NUMNP)
!
!  PRINT NODE VELOCITY
!
      WRITE(IUNIT,1220)
      WRITE(IUNIT, 1180) (N,(URZ(I,N),I=1,2),N=1,NUMNP)
!
!  ELEMENT INFORMATION
!
      WRITE(IUNIT,1270) NUMEL
      WRITE(IUNIT,1330)
      WRITE(IUNIT,1350) (N,(NOD(I,N),I=1,4), N=1,NUMEL)
!
!  BOUNDARY CONDITION
!
      WRITE(IUNIT,1400)
      WRITE(IUNIT,1430) (N,(NBCD(I,N),I=1,2),N=1,NUMNP)
!
!  WRITE STRAIN DISTRIBUTION AT INPUT STAGE
!
      WRITE(IUNIT,1500)
      WRITE(IUNIT,1550) (N,TEPS(N),N=1,NUMEL)
      RETURN
!
 1010 FORMAT(1H1,///,5X,'OUTPUT OF S P I D ',//,5X,A,///)
 1020 FORMAT(5X,'INITIAL INPUT SUMMARY'///)
 1030 FORMAT('  INITIAL STEP NUMBER     = ',I5,/,'FINAL STEP NUMBER       = ',I5,/,'STEP SIZE IN TIME UNIT  = ',F10.5)
 1050 FORMAT('  LIMITING STRAIN RATE    = ',F15.7,/,'  PENALTY CONSTANT        = ',F15.7)
 1070 FORMAT('  DEFORMATION CODE        = ',I5,/,'      IF =1, AXISYMMETRIC ',/,'         =2, PLAIN STRAIN ')
 1110 FORMAT('  FRICTION FACTOR      =',F15.7,/)
 1130 FORMAT(///' NUMBER OF NODAL POINTS  =',I5,/)
 1150 FORMAT(//,'              NODE COORDINATES      ',//,'  NODE NO.        X-COORD        Y-COORD',/)     
 1180 FORMAT(5X,I5,5X,2F15.7)    
 1220 FORMAT(///,'                 NODE VELOCITY',//,'  NODE NO.        X-VELOCITY     Y-VELOCITY'/)  
 1270 FORMAT(//,'  NUMBER OF ELEMENTS    =',I5,/)    
 1330 FORMAT(//,'             ELEMENT CONNECTIVITY',//,'  ELE NO.    I       J       K       L   ',/)                    
 1350 FORMAT(5I7)    
 1400 FORMAT(//,'       BOUNDARY CONDITION CODE',//,'  NODE NO.   X1-CODE     X2-CODE   ',/)
 1430 FORMAT(3I7)
 1500 FORMAT(///,'   STRAIN DISTRIBUTION AT INPUT STAGE',//, '  NODE NO.          STRAIN  ',/)  
 1550 FORMAT(I5,5X,F15.7)    
      END
! ***********************************************************************
      SUBROUTINE PRTSOL
!
!  THIS SUBROUTINE PRINT THE SOLUTION RESULTS
!
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER TITLE*70
      COMMON /TITL/TITLE
      COMMON /INOT/INPT,MSSG,IUNIT,IUNI2
      COMMON /TSTP/NINI,NCUR,NSEND,NITR,DTMAX
      COMMON /MSTR/NUMNP,NUMEL,IPLNAX
      COMMON /RVA1/RZ(2,100),URZ(2,100),FRZ(2,100)
      COMMON /RVA2/EPS(5,100),STS(5,100),TEPS(100)
      COMMON /INVR/NOD(4,100),NBCD(2,100)
!
!  PRINT NODE COORDINATES
!
      WRITE(IUNIT,1010) TITLE,NCUR
      WRITE(IUNIT,1020)
      WRITE(IUNIT,1040) (N,(RZ(I,N), I=1,2), N=1, NUMNP)
!      
!   PRINT NODE VELOCITY, NODAL FORCE
!
      WRITE(IUNIT,1080)
      WRITE(IUNIT,1100) (N,(URZ(I,N),I=1,2),(FRZ(I,N), I=1,2),N=1,NUMNP) 
!
!   STRAIN RATE, STRESS, TOTAL EFFECTIVE STRAIN
!
      WRITE(IUNIT,1130)
      WRITE(IUNIT,1180) (N,(EPS(I,N),I=1,5),N=1,NUMEL)
      WRITE(IUNIT,1230)
      WRITE(IUNIT,1180) (N,(STS(I,N),I=1,5), N=1,NUMEL)
      WRITE(IUNIT,1330)
      WRITE(IUNIT,1360) (N, TEPS(N), N = 1, NUMEL)
      RETURN
!
 1010 FORMAT(1H1,///,5X,'OUTPUT OF S P I D',//,5X,A,//,10X,'SOLUTION AT STEP NUMBER = ',I5,///)
 1020 FORMAT( /,'        NODE COORDINATES',//,'    NODE NO      X-COORD      Y-COORD',/)
 1040 FORMAT(5X,I5,5X,2F15.7)
 1080 FORMAT(///,'               NODAL VELOCITY AND FORCE',//,'  NODE NO      X-VELOCITY    Y-VELOCITY','        X-FORCE       Y-FORCE',//)
 1100 FORMAT(3X,I5,3X,4F15.7)
 1130 FORMAT(///,'                  STRAIN RATE COMPONENTS',//,'  ELE. NO.       E11         E22         E33','            E12         EBAR',//)
 1180 FORMAT(I5,5F15.7)
 1230 FORMAT(//,'                   STRESS COMPONENTS',//,'  ELE. NO.       S11          S22        S33','            S12          SBAR',//)
 1330 FORMAT(///,'      TOTAL EFFECTIVE STRAIN',//,'   ELE. NO.          EFFECTIVE STRAIN',//)
 1360 FORMAT(5X,I5,5X,F15.7)
      END
! ********************************************************************
      SUBROUTINE RSTFIL
!
! GENERATE RESTART FILE
!
      IMPLICIT REAL*8 (A-H, O-Z), INTEGER*4 (I-N)
      CHARACTER TITLE*70
      COMMON /TITL/TITLE
      COMMON /TSTP/NINI,NCUR,NSEND,NITR,DTMAX
      COMMON /RVA1/RZ(2,100),URZ(2,100),FRZ(2,100)
      COMMON /RVA2/EPS(5,100),STS(5,100),TEPS(100)
      COMMON /INVR/NOD(4,100),NBCD(2,100)
      COMMON /DIES/FRCFAC
      COMMON /RIGD/RTOL,ALPH,DIAT
      COMMON /MSTR/NUMNP,NUMEL,IPLNAX
      COMMON /INOT/INPT,MSSG,IUNIT,IUNI2
!
      NN = NCUR+1
      OPEN(IUNI2,FILE='SPID.RST',STATUS='UNKNOWN',FORM='FORMATTED')
      WRITE(IUNI2,1010) TITLE
      WRITE(IUNI2,1040) NCUR,NN,DTMAX
      WRITE(IUNI2,1060) ALPH,DIAT
      WRITE(IUNI2,1080) IPLNAX
      WRITE(IUNI2,1060) FRCFAC
      WRITE(IUNI2,1080) NUMNP
      WRITE(IUNI2,1120) (N,(RZ(I,N),I=1,2),N=1,NUMNP)
      WRITE(IUNI2,1080) NUMEL
      WRITE(IUNI2,1080) (N,(NOD(I,N),I=1,4),N=1,NUMEL)
      WRITE(IUNI2,1160) (N,(NBCD(I,N),I=1,2),N=1,NUMNP)
      WRITE(IUNI2,1120) (N,(URZ(I,N),I=1,2),N=1,NUMNP)
      WRITE(IUNI2,1200) (N,TEPS(N), N=1,NUMEL)
      CLOSE(IUNI2)
      RETURN
!
 1010 FORMAT(1X,A)
 1040 FORMAT(2I10,F20.7)
 1060 FORMAT(3F20.10)
 1080 FORMAT(5I7)
 1120 FORMAT(I5,2F20.10)
 1160 FORMAT(3I7)
 1200 FORMAT(I7,F20.10)
      END
!  *********************************************************
      SUBROUTINE STIFF(B,A,NEQ,MBAND,ITYP)
!
!  STIFFNESS MATRIX GENERATION
! ITYP =1,NEWTON-RAPHSON ITERATION
! ITYP =2,DIRECT ITERATION
!
      IMPLICIT REAL*8 (A-H, O-Z), INTEGER*4 (I-N)
      COMMON /INOT/INPT,MSSG,IUNIT,IUNI2
      COMMON /RVA1/RZ(2,100),URZ(2,100),FRZ(2,100)
      COMMON /RVA2/EPS(5,100),STS(5,100), TEPS(100)
      COMMON /INVR/NOD(4,100),NBCD(2,100)
      COMMON /DIES/FRCFAC
      COMMON /MSTR/NUMNP,NUMEL,IPLNAX
      DIMENSION A(NEQ,1),B(1)
      DIMENSION RZE(2,4),URZE(2,4),NBCDE(2,4),PP(8,8),QQ(8),LM(8)
!
! INITIALIZE LOAD VECTOR, STIFFNESS MATRIX,AND
! NODAL POINT FORCE ARRAY
!
      DO 20 N = 1,NEQ
      B(N) = 0.
      DO 20 I =1,MBAND
      A(N,I)=0.
   20 CONTINUE
      DO 50 N = 1,NUMNP
      DO 50 I=1,2
   50 FRZ(I,N)=0.
!
      DO 200 N=1,NUMEL
!
! CHANGE RZ,URZ,AND NBCD FROM GLOBAL ARRANGEMENT TO ELEMENTAL
! ARRANGEMENT
!
      DO 100 I = 1,4
      I2=I*2
      I1=I2-1
      NE = NOD(I,N)
      RZE(1,I) = RZ(1,NE)
      RZE(2,I) = RZ(2,NE)
      URZE(1,I)= URZ(1,NE)
      URZE(2,I)= URZ(2,NE)
      NBCDE(1,I) = NBCD(1,NE)
      NBCDE(2,I) = NBCD(2,NE)
      LM(I2) = NOD(I,N)*2
      LM(I1) = LM(I2)-1
  100 CONTINUE
!
      CALL ELSHLF(PP,QQ,RZE,URZE,EPS(1,N),TEPS(N),IPLNAX,ITYP,N)
      IF(ITYP .EQ. 1) CALL NFORCE(QQ,FRZ,LM)
      IF(FRCFAC .NE. 0.)CALL FRCBDY(RZE,URZE,NBCDE,TEPS(N),EPS(1,N),QQ,PP,IPLNAX)
      CALL ADDBAN(B,A,NEQ,LM,QQ,PP)
  200 CONTINUE
!
!  APPLY DISPLACEMENT BOUNDARY CONDITION
!
      CALL DISBDY(URZ,NBCD,B,A,NEQ,MBAND,ITYP)
      RETURN
      END
!  ******************************************************************
      SUBROUTINE STRMTX(RZ,B,WDXJ,S,T,IPLNAX,NEL)
!
!    -----CHAPTER 6.4
!    -----EQUATIONS(6.25),(6.27) AND(6.35)
!    EVALUATE STRAIN RATE MATRIX OF QUADRILATERAL ELEMENT
!
!     B(4,8)   : STRAIN RATE MATRIX
!     RZ(2,4)  : NODE COORDINATES  
!     (S,T)    : NATURAL COORDINATE 
!
      IMPLICIT REAL*8 (A-H, O-Z), INTEGER*4 (I-N)
      COMMON /INOT/INPT,MSSG,IUNIT,IUNI2
      DIMENSION RZ(2,1), B(4,1)
!
      R12 = RZ(1,1)-RZ(1,2)
      R13 = RZ(1,1)-RZ(1,3)
      R14 = RZ(1,1)-RZ(1,4)
      R23 = RZ(1,2)-RZ(1,3)
      R24 = RZ(1,2)-RZ(1,4)
      R34 = RZ(1,3)-RZ(1,4)
!
      Z12 = RZ(2,1) - RZ(2,2)
      Z13 = RZ(2,1) - RZ(2,3)
      Z14 = RZ(2,1) - RZ(2,4)
      Z23 = RZ(2,2) - RZ(2,3)
      Z24 = RZ(2,2) - RZ(2,4)
      Z34 = RZ(2,3) - RZ(2,4)
!
!  REFER EQUATION (6.34)
!
      DXJ8 = ( (R13*Z24-R24*Z13)+(R34*Z12-R12*Z34)*S+(R23*Z14-R14*Z23)*T)
      DXJ = DXJ8/8.
      IF(DXJ .GT. 0.) GO TO 10
      WRITE(MSSG,1010) NEL
      WRITE(MSSG,1030) DXJ,S,T
      STOP
   10 CONTINUE
!
      X1 = ( Z24-Z34*S-Z23*T) / DXJ8
      X2 = (-Z13+Z34*S+Z14*T) / DXJ8 
      X3 = (-Z24+Z12*S-Z14*T) / DXJ8 
      X4 = ( Z13-Z12*S+Z23*T) / DXJ8 
!
      Y1 = (-R24+R34*S+R23*T) / DXJ8
      Y2 = ( R13-R34*S-R14*T) / DXJ8 
      Y3 = ( R24-R12*S+R14*T) / DXJ8 
      Y4 = (-R13+R12*S-R23*T) / DXJ8 
!
      DO 20 I=1,4
      DO 20 J=1,8
      B(I,J) = 0.
   20 CONTINUE
!
      B(1,1) = X1
      B(1,3) = X2
      B(1,5) = X3
      B(1,7) = X4
      B(2,2) = Y1
      B(2,4) = Y2
      B(2,6) = Y3
      B(2,8) = Y4
!
      WDXJ = DXJ
      IF(IPLNAX .NE. 1) GO TO 40
      Q1 = (1.-S)*(1.-T)*0.25
      Q2 = (1.+S)*(1.-T)*0.25
      Q3 = (1.+S)*(1.+T)*0.25
      Q4 = (1.-S)*(1.+T)*0.25
!
      R = Q1*RZ(1,1)+Q2*RZ(1,2)+Q3*RZ(1,3)+Q4*RZ(1,4)
      B(3,1) = Q1 / R
      B(3,3) = Q2 / R
      B(3,5) = Q3 / R
      B(3,7) = Q4 / R
      WDXJ = WDXJ * R
!
   40 CONTINUE
      B(4,1) = Y1
      B(4,3) = Y2
      B(4,5) = Y3
      B(4,7) = Y4
      B(4,2) = X1
      B(4,4) = X2
      B(4,6) = X3
      B(4,8) = X4
      RETURN
!
 1010 FORMAT(/,'SORRY,NEGATIVE JACOBIAN DETECTED AT ELEMENT NO.',I5)
 1030 FORMAT(' DXJ,S,T = ',3F15.7)
      END 
! *****************************************************************      
      SUBROUTINE VSPLON(QQ,PP,B,URZ,EPS,WDXJ,IDREC)
!
!    ---- CHAPTER 6.5 AND 7.1
!    ---- EQUATIONS (6.43),(6.44)
!    REDUCED INTEGRATION OF VOLUME STRAIN RATE
!
!    PP = ELEMENTAL STIFFNESS MATRIX
!    QQ = ELEMENTAL LOAD VECTOR
!    B  = STRAIN RATE MATRIX 
!
      IMPLICIT REAL*8 (A-H, O-Z), INTEGER*4 (I-N)
      COMMON /RIGD/RTOL,ALPH,DIAT
      DIMENSION PP(8,8), QQ(8),B(4,8), URZ(1),EPS(1)
      DIMENSION D(6),XX(8)
      DATA D/3*0.6666666666666667D0, 3*0.3333333333333333D0/
!
!   GENERATE DILATATIONAL STRAIN RATE MATRIX
!
      DO 20 I = 1,8
      XX(I) = B(1,I)+B(2,I)+B(3,I)
   20 CONTINUE
!
!  CALCULATE STRAIN RATE COMPONENTS
!
      DO 40 I = 1,5
      EPS(I) = 0.
   40 CONTINUE
      XVOL = 0.
      DO 60 J =1,8
      XVOL = XVOL +XX(J)*URZ(J)
      DO 60 I= 1,4
      EPS(I) = EPS(I) +B(I,J)*URZ(J)
   60 CONTINUE
      EB2 = (EPS(1)**2+EPS(2)**2+EPS(3)**2)*D(1)+EPS(4)**2*D(4)
      EPS(5) = DSQRT(EB2)
!
!  EVALUATE VOLUMETRIC CONTRIBUTION OF STIFFNESS MATRIX
!
      DO 80 I = 1,8
      IF(IDREC .EQ. 1)QQ(I) =QQ(I) -DIAT*WDXJ*XVOL*XX(I)
      TEM = DIAT*WDXJ*XX(I) 
      DO 80 J = I,8
      PP(I,J) = PP(I,J)+TEM*XX(J)
      PP(J,I) = PP(I,J)
   80 CONTINUE
      RETURN
      END
!  **************************************************************
      SUBROUTINE VSPLST(QQ,PP,B,URZ,TEPS,WDXJ,IDREC)
!
!   ---- CHAPTER 6.5
!   ---- EQUATIONS (6.43),(6.44),(6.46)
! FOUR POINTS INTEGRATION OF VOLUME STRAIN RATE
!
! PP  = ELEMENTAL STIFFNESS MATRIX
! QQ  = ELEMENTAL LOAD VECTOR
! B   = STRAIN RATE MATRIX
!
      IMPLICIT REAL*8 (A-H, O-Z), INTEGER*4 (I-N)
      COMMON /TSTP/NINI,NCUR,NSEND,NITR,DTMAX
      COMMON /RIGD/RTOL,ALPH,DIAT
      DIMENSION PP(8,8),QQ(8),B(4,8),URZ(1)
      DIMENSION D(6),FDV(8),E(4),XX(8)
      DATA D/3*0.6666666666666667D0,3*0.3333333333333333D0/
!
! ELIMINATE DIALATATIONAL COMPONENT FROM STRAIN RATE MATRIX
!
      DO 20 I= 1,8
      XX(I) = (B(1,I)+B(2,I)+B(3,I))/3.
   20 CONTINUE
      DO 40 I = 1,8
      DO 40 J=1,3
      B(J,I) = B(J,I)-XX(I)
   40 CONTINUE
!
! CALCULATE STRAIN RATE
!
      DO 60  J = 1,4
      E(J) = 0.
      DO 60  I = 1,8
      E(J) = E(J) +B(J,I)*URZ(I)
   60 CONTINUE
      EFSR2 = D(1)*E(1)*E(1)+ D(2)*E(2)*E(2)+ D(3)*E(3)*E(3)+D(4)*E(4)*E(4)
      IF(NITR .EQ. 1 .AND. NCUR .EQ. NINI .AND. IDREC .EQ. 2)EFSR2 = (ALPH*100.)**2
      ALPH2 = ALPH**2
      IF(EFSR2 .LT. ALPH2) EFSR2 =ALPH2
      EFSR = DSQRT(EFSR2)
      CALL FLWSTS(EFSTS,STRAT,TEPS,EFSR)
!
! CALCULATE FIRST DERIVATE OF EFSR **2
!
      DO 80  I = 1, 8
      FDV(I) = 0.
      DO 80  J = 1, 4
      FDV(I) = FDV(I)+D(J)*E(J)*B(J,I)
   80 CONTINUE
!
!  ADD POINT CONTRIBUTION TO STIFFNESS MATRIX
!
      F1 = EFSTS / EFSR *WDXJ
      IF(IDREC .EQ. 2) GO TO 200
      F2 = STRAT /EFSR2*WDXJ-F1 / EFSR2
      DO 120  I = 1, 8
      QQ(I) = QQ(I)- FDV(I)*F1
      DO 110  J = I, 8
      TEM = 0.
      DO 100 K=1,4
      TEM = TEM +D(K)*B(K,I)*B(K,J)
  100 CONTINUE
      PP(I,J) = PP(I,J) + TEM*F1
      IF(EFSR2 .LT. ALPH2) GO TO 105
      PP(I,J) = PP(I,J)+FDV(I)*FDV(J)*F2
  105 PP(J,I) = PP(I,J)
  110 CONTINUE
  120 CONTINUE
      RETURN
!
  200 CONTINUE
      DO 300  I = 1, 8
      DO 280  J = I, 8
      TEM = 0.
      DO 250 K = 1,4
      TEM = TEM +D(K)*B(K,I)*B(K,J)
  250 CONTINUE
      PP(I,J) = PP(I,J) +TEM*F1
      PP(J,I)=PP(I,J)
  280 CONTINUE
  300 CONTINUE
      RETURN
      END
