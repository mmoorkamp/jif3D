      SUBROUTINE EPOL(T, NX, NZ, IONOS, IATM, DXJ, 
     @DZI,RO, RIONOS, HX_REAL, HX_IMAG, EY_REAL, EY_IMAG)
C
c               MODELLING MT 2-D (from Tartis, 1987)
C
C                       E-POLARISATION ( MODE TE )
C
C
C       INPUT :       
C                     - T      Period in 1/s
C                     - DXJ    Grid cell size in x direction in km
C                     - DZI    Grid cell size in z direction in km
C                     - RO     Resistivities
C                     - NZ     Total number of layers ( <--> DZI )
C                     - NX     Total number of horizontal cells
C                     - IONOS    Number of layers in the ionosphere
C                     - IATM     Number of layers in the athmosphere
C                     - RIONOS   Resistivity of the ionosphere
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCC 
      COMPLEX*16 :: E0_B,V,V1,E0Z,CHE
      REAL*8  ::    T, RIONOS, EPS_V
      INTEGER*4 ::  NZ, NX, IONOS, IATM
      REAL*8, dimension(NX+1)  ::   DXJ
      REAL*8, dimension(NZ+1)  ::   DZI
      REAL*8, dimension(NX*(NZ-IONOS-IATM)) :: HX_IMAG,
     @HX_REAL,EY_REAL, EY_IMAG, RO
      REAL*8, allocatable, dimension (:,:) :: RO2_B
      REAL*8, allocatable, dimension (:) :: BUFF,DXM,DZM,DXMH,DZMH 
      COMPLEX*16, allocatable, dimension (:,:) :: HX2_B,EY2_B,HZ2_B,E,
     @UT,A 
      COMPLEX*16, allocatable, dimension (:) :: EPLUS, B, AIT, BI, CHE1
       INTEGER*4 alloc_status, dealloc_status
CCCCCCCCCCCCCCCCCCCCCCCCCCC 
C...... initialising parameters
C
      PI=DATAN(1.D0)*4.D0
      C=3.D0*1.D10
      AMU=4.D0*PI*1.D-7
      EPS_V = 1.D-15
      E0_B =(1.d0,0.d0)

C...... frequency
      F=1.D0/T
C
C
      RION=RIONOS
      ISED=NZ-IONOS-IATM
      N0=IONOS
      N1=IATM
      N2=ISED
      N=NZ
      
CCCCCCCCCCCCCCCCCCCCCCCCCCC      
C     Make the two dimensional arrays and allocate memory

      allocate(RO2_B(ISED+1,NX+1),stat=alloc_status)
      allocate(HX2_B(ISED+1,NX+1),stat=alloc_status)
      allocate(HZ2_B(ISED+1,NX+1),stat=alloc_status)
      allocate(EY2_B(ISED+1,NX+1),stat=alloc_status)
      allocate(E(NZ+7,2),stat=alloc_status)
      
      allocate(A(NZ+1,NZ+1),stat=alloc_status)
      allocate(AIT(NZ+1),stat=alloc_status)
      allocate(B(NZ+1),stat=alloc_status)
      allocate(CHE1(NZ+1),stat=alloc_status)
      allocate(BI(NZ+1),stat=alloc_status)
      allocate(UT(NZ+1,NX+1),stat=alloc_status)
      
      allocate(DXM(NX),stat=alloc_status)
      allocate(DZM(NZ),stat=alloc_status)
      allocate(DXMH(NX),stat=alloc_status)
      allocate(DZMH(NZ),stat=alloc_status)
      
      allocate(BUFF(2*NZ*(NZ-1)*NX),stat=alloc_status)
      allocate(EPLUS(NZ+1),stat=alloc_status)
      
      DO 500 I=1,(NZ-IONOS-IATM)
       DO 500 J=1,NX
          RO2_B(I,J) =RO((J-1)*(NZ-IONOS-IATM) + I)
500   CONTINUE 

CCCCCCCCCCCCCCCCCCCCCCCCCCC           
C       Transformation to cgs units and determination of the 
C       structure
C
      DO 25 I=1,NZ
25       DZI(I)=DZI(I)*1.D2  
      DO 701 J=1,NX
      DXJ(J)=DXJ(J)*1.D2
      DO 21 I=1,N2
       RO2_B(I,J)=RO2_B(I,J)*1.D11
21       CONTINUE
701       CONTINUE
C
      RION=RION*1.D11
C
      NM1=N-1
      NM2=2*NM1
      DO 24 I=1,NM1
      DZMH(I)=1./DZI(I)+1./DZI(I+1)
24       DZM(I)=.5D0*(DZI(I)+DZI(I+1))
      NXM1=NX-1
      DO 22 J=1,NXM1
      DXMH(J)=1.D0/DXJ(J)+1.D0/DXJ(J+1)
22       DXM(J)=.5D0*(DXJ(J)+DXJ(J+1))
C
C       CALCUL DE UT
C
      IF(N0.EQ.0) GOTO 102
      U=8.D0*PI*PI*F
      DO 101 J=1,NX
      IF(RION.EQ.0.0) THEN
       RION = EPS_V
      ENDIF
       UT(N0,J)=DCMPLX(0.D0,U/RION)
101   CONTINUE
      U1=2.D0*PI*F
      U=-U1*U1/(C*C)
      NOM1=N0-1
      DO 702 I=1,NOM1
      DO 103 J=1,NX
103       UT(I,J)=U
702       CONTINUE
102       IF(N1.EQ.0) GOTO 104
      U1=2.D0*PI*F
      U=-U1*U1/(C*C)
      I1=N0+1
      I2=N0+N1
      DO 703 I=I1,I2
      DO 105 J=1,NX
105       UT(I,J)=U
703       CONTINUE
104       U=8.D0*PI*PI*F
      I1=N0+N1+1
      DO 704 I=I1,NZ
      IV=I-N0-N1
      DO 106 J=1,NX
       IF(RO2_B(IV,J).EQ.0.0) THEN
         RO2_B(IV,J) = EPS_V
       ENDIF
       UT(I,J)=DCMPLX(0.D0,U/RO2_B(IV,J))
106   CONTINUE
704       CONTINUE
C
      E0Z=E0_B/DZI(1)
C
      KECS=1
      DO 706 I=1,N
      B(I)=0.D0
      DO 110 J=1,N
110       A(I,J)=0.D0
706       CONTINUE
C
C
C              DEFINITION DE LA MATRICE A1 DE DIMENSION N
C
C
      DO 115 I=1,NM1
      V1=UT(I,1)*DZI(I)*DXJ(1)+UT(I+1,1)*DZI(I+1)*
     #DXJ(1)
      U2=DZM(I)/DXJ(1)+.5D0*DXJ(1)*DZMH(I)
      A(I,I)=.25D0*V1+U2
115       A(I,I+1)=-.5D0*DXJ(1)/DZI(I+1)
      A(1,N)=-DZM(1)/DXJ(1)
      V1=UT(1,1)*DZI(1)*DXJ(1)+UT(1,2)*DZI(1)*DXJ(2)+
     #UT(2,1)*DZI(2)*DXJ(1)+UT(2,2)*DZI(2)*DXJ(2)
      V2=DZM(1)*DXMH(1)+DXM(1)*DZMH(1)
      A(N,N)=.25D0*V1+V2
      A(NM1,N)=0.D0
C
C      DETERMINATION DE LA MATRICE B1
C
      B(1)=.5D0*E0Z*DXJ(1)
      B(N)=E0Z*DXM(1)
C
C
C       TRIANGULISATION DE LA MATRICE A
C
      L=1
      IL=1
10       V=1.D0/A(1,1)
      B(1)=B(1)*V
      BI(IL)=B(1)
      DO 20 J=2,N
      AIT(J)=A(1,J)
      A(1,J)=A(1,J)*V
20       B(J)=B(J)-B(1)*AIT(J)
      DO 707 I=2,N
      DO 30 J=I,N
30       A(I,J)=A(I,J)-A(1,J)*AIT(I)
707       CONTINUE
C
      IL=IL+1
      DO 90 J=2,N
C      CALL WRITEC(A(1,J),KECS,2)
      BUFF(KECS)=DREAL(A(1,J))
      BUFF(KECS+1)=DIMAG(A(1,J))
90       KECS=KECS+2
      IF(IL.NE.N) GOTO 91
      KEC=KECS
C       CALL WRITEC(BI(1),KECS,NM2)
      DO 9998 KLI=1,NM1
      BUFF(KEC+KLI-1)=DREAL(BI(KLI))
      BUFF(KEC+KLI)=DIMAG(BI(KLI))
9998  KEC=KEC+1
      KECS=KECS+NM2
91       CONTINUE
C
C
C       PASSAGE DE LA MATRICE AK A LA MATRICE AK+1
C
C
      DO 708 I=1,NM1
      B(I)=B(I+1)
      DO 40 J=1,NM1
40       A(I,J)=A(I+1,J+1)
708       CONTINUE
C
      IF(IL.EQ.N.AND.L.EQ.NX) GOTO 60
C
      DO 45 I=1,N
45       A(I,N)=0.D0
      B(N)=0.D0
      A(NM1,N)=-DXM(L)/DZI(IL)
      IF(IL.NE.N) GOTO 50
      IL=1
      L=L+1
      A(NM1,N)=0.D0
      B(N)=E0Z*DXM(L)
50       A(1,N)=-DZM(IL)/DXJ(L)
      IF(L.EQ.NX) GOTO 55
      V1=UT(IL,L)*DZI(IL)*DXJ(L)+UT(IL+1,L)*DZI(IL+1)*DXJ(L)+
     #UT(IL,L+1)*DZI(IL)*DXJ(L+1)+
     #UT(IL+1,L+1)*DZI(IL+1)*DXJ(L+1)
      V2=DZM(IL)*DXMH(L)+DXM(L)*DZMH(IL)
      A(N,N)=.25D0*V1+V2
      GOTO 10
55       CONTINUE
      V1=UT(IL,L)*DZI(IL)*DXJ(L)+UT(IL+1,L)*DZI(IL+1)*DXJ(L)
      V2=DZM(IL)/DXJ(L)+.5D0*DXJ(L)*DZMH(IL)
      A(N,N)=.25*V1+V2
      A(NM1,N)=-.5D0*DXJ(L)/DZI(IL)
      IF(IL.NE.1) GOTO 10
      A(NM1,N)=0.D0
      B(N)=.5D0*E0Z*DXJ(L)
      GOTO 10
C
C       TRIANGULATION DE LA DERNIERE MATRICE
C
60        IL=0
70       IL=IL+1
      I2=IL+1
C
      V=1.D0/A(IL,IL)
      B(IL)=B(IL)*V
      DO 75 J=I2,NM1
      V1=A(IL,J)
      AIT(J)=V1
      A(IL,J)=V1*V
75       B(J)=B(J)-V1*B(IL)
      DO 709 I=I2,NM1
      DO 80 J=I,NM1
80       A(I,J)=A(I,J)-A(IL,J)*AIT(I)
709       CONTINUE
      IF(I2.NE.NM1) GOTO 70
C
C       EXTRACTION DE RO DE UT
C
      I3=N0+N1
      J=NX+1
85       I=ISED+1
      J=J-1
87       I=I-1
      II3=I+I3
      RO2_B(I,J)=U/DIMAG(UT(II3,J))
      IF(I.GT.1) GOTO 87
      IF(J.GT.1) GOTO 85
C
C
C       CALCUL DU CHAMP  EY
C
C
      I4=N0+N1+1
      NM1=N-1
      DZ1=DZI(I4)
      U1=T/(DZ1*2.D0*PI)
      DO 711 I=1,NM1
      DO 1 J=1,2
1       E(I,J)=0.D0
711       CONTINUE
      I=NM1
      J=NX+1
      E(I,1)=B(NM1)/A(I,I)
2       CONTINUE
      I=I-1
      IF(I.EQ.0) GOTO 3
      I1=I+1
      E(I,1)=B(I)-A(I,I1)*E(I1,1)
4       CONTINUE
      I1=I1+1
      IF(I1.GT.NM1) GOTO 2
      E(I,1)=E(I,1)-A(I,I1)*E(I1,1)
      GOTO 4
3       CONTINUE
      DO 1005 L=1,NM1
1005       E(L,2)=E(L,1)
c       IF(J.EQ.NX+1)
c    @  WRITE(3,'(2E13.7)') (E(L,2),L=1,NM1)
      DO 1006 L=1,NM1
1006       E(L,1)=0.D0


        U=2.D0*PI/T*AMU
        DO 3333 IZZ=N0+N1,NM1
        U1=DZI(IZZ+1)/1.D2*2.D0*PI/T*AMU
        CHE=(E(IZZ+1,2)-E(IZZ,2))/DCMPLX(0.D0,U1)
        IF(J.LE.NX) THEN
        EY2_B(IZZ-N0-N1+1,J)=.5D0*(E(IZZ,2)+EPLUS(IZZ))
        HZ2_B(IZZ-N0-N1+1,J)=(EPLUS(IZZ)-E(IZZ,2))/DXJ(J)*1.D2
     @/DCMPLX(0.D0,-U)
        HX2_B(IZZ-N1-N0+1,J)=.5D0*(CHE+CHE1(IZZ))
        ENDIF
        EPLUS(IZZ)=E(IZZ,2)
        CHE1(IZZ)=CHE
 3333       CONTINUE
C
C
      J=J-1
      IF(J.EQ.0) GOTO 4000
      IF(J.NE.NX) GOTO 8
       KECS=KECS-N*NM2
      GOTO 11
8       IF(J.NE.1) GOTO 1022
      KECS=1
      GOTO  11
1022       CONTINUE
      NF2=N*2
      KECS=KECS-NF2*NM2
11       CONTINUE
C
      DO 710 I=1,NM1
      DO 1012 JV=1,NM1
C       READEC(A(I,JV),KECS,2)
      AR=BUFF(KECS)
      AI=BUFF(KECS+1)
      A(I,JV)=DCMPLX(AR,AI)
1012       KECS=KECS+2
710       CONTINUE
      KEC=KECS
C       CALL READEC(B(1),KECS,NM2)
      DO 9996 KIJ=1,NM1
      AR=BUFF(KEC+KIJ-1)
      AI=BUFF(KEC+KIJ)
      B(KIJ)=DCMPLX(AR,AI)
9996  KEC=KEC+1
      KECS=KECS+NM2
C
C
      I=N
1013       I=I-1
      IF(I.EQ.0) GOTO 3
C
      IF(I.EQ.NM1) GOTO 1014
C
      I1=I+1
      I2=1
      E(I,1)=B(I)-A(I,I2)*E(I1,1)
1016       IF(I1.EQ.NM1) GOTO 1015
      I1=I1+1
      I2=I2+1
      E(I,1)=E(I,1)-A(I,I2)*E(I1,1)
      GOTO 1016
1015       I1=0
1017       I2=I2+1
      I1=I1+1
      IF(I1.GT.I) GOTO 1013
      E(I,1)=E(I,1)-A(I,I2)*E(I1,2)
      GOTO 1017
1014       I1=1
      I2=I1
      E(I,1)=B(I)-A(I,I1)*E(I1,2)
      GOTO 1017
4000       CONTINUE

C
C   NORMALISATION PAR RAPPORT A CHAMPX(1,1)
C
      CHE=HX2_B(1,1)
C
      DO 1201 J=1,NX
      DO 1201 I=1,NZ-N0-N1
c  en nT
      HX2_B(I,J)=-HX2_B(I,J)/CHE            
      HZ2_B(I,J)=-HZ2_B(I,J)/CHE
c en microV/m/nT
      EY2_B(I,J)=-EY2_B(I,J)/CHE/4.D-4/PI 
 1201       CONTINUE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 3000 I=1,NZ
 3000       DZI(I)=DZI(I)*1.D-2
      DO 3001 J=1,NX
 3001       DXJ(J)=DXJ(J)*1.D-2
 
      DO 3002 I=1,NZ-N0-N1
        DO 3002 J=1,NX
        
C           Transfer the E-Field to V/m
            EY2_B(I,J) = EY2_B(I,J)/1.D6
C           Transfer the B to a H field and nT to T
            HX2_B(I,J) = (HX2_B(I,J))/(4*PI*1.D2)
        
            HX_REAL((J-1)*(NZ -N0-N1) + I) = DREAL(HX2_B(I,J))
            HX_IMAG((J-1)*(NZ -N0-N1) + I) = DIMAG(HX2_B(I,J))
            
            EY_REAL((J-1)*(NZ -N0-N1) + I) = DREAL(EY2_B(I,J))
            EY_IMAG((J-1)*(NZ -N0-N1) + I) = DIMAG(EY2_B(I,J))
            
3002  CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCC

      deallocate(RO2_B,stat=dealloc_status)
      deallocate(EY2_B,stat=dealloc_status)
      deallocate(HX2_B,stat=dealloc_status)
      deallocate(HZ2_B,stat=dealloc_status)
      deallocate(E,stat=dealloc_status)
      
      deallocate(A,stat=dealloc_status)
      deallocate(AIT,stat=dealloc_status)
      deallocate(B,stat=dealloc_status)
      deallocate(CHE1,stat=dealloc_status)
      deallocate(BI,stat=dealloc_status)
      deallocate(UT,stat=dealloc_status)
      
      deallocate(DXM,stat=dealloc_status)
      deallocate(DZM,stat=dealloc_status)
      deallocate(DXMH,stat=dealloc_status)
      deallocate(DZMH,stat=dealloc_status)
      
      deallocate(BUFF,stat=dealloc_status)
      deallocate(EPLUS,stat=dealloc_status)
CCCCCCCCCCCCCCCCCCCCCCCCCCC
      RETURN
      END
