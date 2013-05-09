      SUBROUTINE HPOL(T, NX, NZ, DXJ, DZI, RO, HY_REAL, HY_IMAG, 
     @EX_REAL, EX_IMAG, EZ_REAL, EZ_IMAG)
C
C
C                     MODELING MT 2-D (from Tartis, 1987, modified Heincke 2006)
C
C                        H-POLARISATION ( MODE TM )
C       
C       In this case, the horizontal component of the magnetic field is constant at the surface
C       (The first row corresponds to the air) 
C
C       INPUT :       
C                     - T      Period in 1/s
C                     - DXJ    Grid cell size in x direction in km
C                     - DZI    Grid cell size in z direction in km
C                     - RO     Resistivities
C                     - NZ     Total number of layers ( <--> DZI )
C                     - NX     Total number of horizontal cells
C                     - HY_REAL,HY_IMAG   HY-component of the field
c                     - EX_REAL,EX_IMAG   EX-component of the field
c                     - EZ_REAL,EZ_IMAG   EZ-component of the field
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCC      
 
      COMPLEX*16 :: H0   
      REAL*8  ::    T
      INTEGER*4 ::  NZ, NX
      REAL*8, dimension(NX+1)  ::   DXJ
      REAL*8, dimension(NZ+1)  ::   DZI
      REAL*8, dimension(NX*NZ)  :: RO, HY_REAL,HY_IMAG,EX_REAL,EX_IMAG 
      REAL*8, dimension(NX*NZ)  :: EZ_REAL,EZ_IMAG 
      REAL*8, allocatable, dimension (:,:) :: RO2
      REAL*8, allocatable, dimension (:) :: BUFF, DXM, DZM 
      COMPLEX*16, allocatable, dimension (:,:) :: HY_B2,EX_B2,A,H,EZ_B2 
      COMPLEX*16, allocatable, dimension (:) ::  B,AIT,BI,EXX,EXX1,
     @HPLUS
      COMPLEX*16 V,V1,HYI1,HYI2,AA,BB
      INTEGER*4 alloc_status, dealloc_status
CCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      
      PI=ATAN(1.D0)*4.D0
      h0=(1.d0,0.d0)
      
C.... FREQUENCY
      F=1./T
C
C
      N=NZ
      NM1=N-1
 
CCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Make the two dimensional arrays and allocate memory

      allocate(RO2(NZ+1,NX+1),stat=alloc_status)
      allocate(HY_B2(NZ+1,NX+1),stat=alloc_status)
      allocate(EX_B2(NZ+1,NX+1),stat=alloc_status)
      allocate(EZ_B2(NZ+1,NX+1),stat=alloc_status)
      allocate(H(NZ+1,2),stat=alloc_status)
      
      allocate(A(NZ+1,NZ+1),stat=alloc_status)
      allocate(B(NZ+1),stat=alloc_status)
      allocate(AIT(NZ+1),stat=alloc_status)
      allocate(BI(NZ+1),stat=alloc_status)
      allocate(EXX(NZ+1),stat=alloc_status)
      allocate(EXX1(NZ+1),stat=alloc_status)
        
      allocate(DXM(NX),stat=alloc_status)
      allocate(DZM(NZ),stat=alloc_status)
      
      allocate(BUFF(2*NZ*(NZ-1)*NX),stat=alloc_status)
      allocate(HPLUS(NZ+1),stat=alloc_status)
       
      DO 500 I=1,NZ
        DO 500 J=1,NX
          RO2(I,J) =RO((J-1)*NZ + I)
500   CONTINUE    

CCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Transformation in cgs units and determination of the structure
C
      DO 20 I=1,N
20       DZI(I)=DZI(I)*1.D2
      DO 21 J=1,NX
21       DXJ(J)=DXJ(J)*1.D2
C
      DO 301 J=1,NX
      DXM(J)=.5D0*(DXJ(J)+DXJ(J+1))
      DO 200 I=1,N
      DZM(I)=.5*(DZI(I)+DZI(I+1))
200       RO2(I,J)=RO2(I,J)*1.D11
301       CONTINUE
C
C
C
C
      KECS=1
      DO 302 I=1,N
      B(I)=0.D0
      DO 110 J=1,N
110       A(I,J)=0.D0
302       CONTINUE
C
C
C              DEFINITION DE LA MATRICE A1 DE DIMENSION N
C
C
      DO 111 I=1,N
111       RO2(I,NX+1)=RO2(I,NX)
      DXJ(NX+1)=DXJ(NX)
      DZI(N+1)=DZI(N)
      DXM(NX)=DXJ(NX)
      U=8.D0*PI*PI*F
      DO 112 I=1,NM1
      U1=DXJ(1)/DZI(I)
      U2=DXJ(1)/DZI(I+1)
      U3=DZI(I)/DXJ(1)
      U4=DZI(I+1)/DXJ(1)
      UI=.5D0*U*DXJ(1)*DZM(I)
      UR=.5D0*(RO2(I,1)*(U1+U3)+RO2(I+1,1)*(U2+U4))
      A(I,I)=DCMPLX(UR,UI)
112       A(I,I+1)=-.5D0*RO2(I+1,1)*U2
C
      J=1
      I=1
      A(NM1,N)=0.D0
      U1=RO2(I,J)*DZI(I)/DXJ(J)
      U2=RO2(I,J)*DXJ(J)/DZI(I)
      U3=RO2(I+1,J)*DZI(I+1)/DXJ(J)
      U4=RO2(I+1,J)*DXJ(J)/DZI(I+1)
      U5=RO2(I,J+1)*DZI(I)/DXJ(J+1)
      U6=RO2(I,J+1)*DXJ(J+1)/DZI(I)
      U7=RO2(I+1,J+1)*DZI(I+1)/DXJ(J+1)
      U8=RO2(I+1,J+1)*DXJ(J+1)/DZI(I+1)
      A(1,N)=-.5D0*(U1+U3)
      UR=.5D0*(U1+U2+U3+U4+U5+U6+U7+U8)
      UI=U*DZM(I)*DXM(J)
      A(N,N)=DCMPLX(UR,UI)
C
C       DEFINITION DE LA MATRICE B1
C
C
      B(1)=.5D0*H0*U2
      B(N)=.5D0*H0*(U2+U6)
C
C       TRIANGULISATION DE LA MATRICE A
C
      L=1
      IL=1
115       V=1.D0/A(1,1)
      B(1)=B(1)*V
      BI(IL)=B(1)
      DO 120 J=2,N
      AIT(J)=A(1,J)
      A(1,J)=A(1,J)*V
120       B(J)=B(J)-B(1)*AIT(J)
      DO 303 I=2,N
      DO 121 J=I,N
121       A(I,J)=A(I,J)-A(1,J)*AIT(I)
303       CONTINUE
C
      IL=IL+1
      DO 122 J=2,N
C      CALL WRITEC(A(1,J),KECS,2)
      BUFF(KECS)=DREAL(A(1,J))
      BUFF(KECS+1)=DIMAG(A(1,J))
122       KECS=KECS+2
      IF(IL.NE.N) GOTO 1220
      M=2*NM1
      KEC=KECS
C       CALL WRITEC(BI(1),KECS,M)
      DO 9998 KLI=1,NM1
      BUFF(KEC+KLI-1)=DREAL(BI(KLI))
      BUFF(KEC+KLI)=DIMAG(BI(KLI))
9998  KEC=KEC+1
      KECS=KECS+M
1220       CONTINUE
C
C
C       PASSAGE DE LA MATRICE AK A LA MATRICE AK+1
C
C
      DO 304 I=1,NM1
      B(I)=B(I+1)
      DO 123 J=1,NM1
123       A(I,J)=A(I+1,J+1)
304       CONTINUE
C
      IF(IL.EQ.N.AND.L.EQ.NX) GOTO 130
C
      DO 124 I=1,N
124       A(I,N)=0.D0
      B(N)=0.D0
C
      U1=RO2(IL,L)*DZI(IL)/DXJ(L)
      U2=RO2(IL,L)*DXJ(L)/DZI(IL)
      U3=RO2(IL+1,L)*DZI(IL+1)/DXJ(L)
      U4=RO2(IL+1,L)*DXJ(L)/DZI(IL+1)
      U5=RO2(IL,L+1)*DZI(IL)/DXJ(L+1)
      U6=RO2(IL,L+1)*DXJ(L+1)/DZI(IL)
      U7=RO2(IL+1,L+1)*DZI(IL+1)/DXJ(L+1)
      U8=RO2(IL+1,L+1)*DXJ(L+1)/DZI(IL+1)
      A(NM1,N)=-.5D0*(U2+U6)
      IF(IL.NE.N) GOTO 125
      IL=1
      L=L+1
      U1=RO2(IL,L)*DZI(IL)/DXJ(L)
      U2=RO2(IL,L)*DXJ(L)/DZI(IL)
      U3=RO2(IL+1,L)*DZI(IL+1)/DXJ(L)
      U4=RO2(IL+1,L)*DXJ(L)/DZI(IL+1)
      U5=RO2(IL,L+1)*DZI(IL)/DXJ(L+1)
      U6=RO2(IL,L+1)*DXJ(L+1)/DZI(IL)
      U7=RO2(IL+1,L+1)*DZI(IL+1)/DXJ(L+1)
      U8=RO2(IL+1,L+1)*DXJ(L+1)/DZI(IL+1)
      A(NM1,N)=0.D0
      B(N)=.5D0*H0*(U2+U6)
125       A(1,N)=-.5D0*(U1+U3)
      UR=.5D0*(U1+U2+U3+U4+U5+U6+U7+U8)
      UI=U*DZM(IL)*DXM(L)
      A(N,N)=DCMPLX(UR,UI)
      IF(L.NE.NX) GOTO 115
      A(N,N)=.5D0*A(N,N)
      A(NM1,N)=.5D0*A(NM1,N)
      B(N)=.5D0*B(N)
      GOTO 115
C
C       TRIANGULATION DE LA DERNIERE MATRICE
C
130        IL=0
131       IL=IL+1
      I2=IL+1
C
      V=1.D0/A(IL,IL)
      B(IL)=B(IL)*V
      DO 133 J=I2,NM1
      V1=A(IL,J)
      AIT(J)=V1
      A(IL,J)=V1*V
133       B(J)=B(J)-V1*B(IL)
      DO 305 I=I2,NM1
      DO 135 J=I,NM1
135       A(I,J)=A(I,J)-A(IL,J)*AIT(I)
305       CONTINUE
      IF(I2.NE.NM1) GOTO 131
C
C
C       CALCUL DU CHAMP MAGNETIQUE HY
C
C
C
C       CALCUL DES CONSTANTES DEN,CH0,CH1,CH2,CH3
C
C      U01=1.D0/(4.D0*PI)
C
      DO 306 I=1,NM1
      DO 137 J=1,2
137       H(I,J)=0.D0
306       CONTINUE
C
       I=NM1
      J=NX+1
      H(I,1)=B(NM1)/A(I,I)
138       I=I-1
      IF(I.EQ.0) GOTO 140
      I1=I+1
      H(I,1)=B(I)-A(I,I1)*H(I1,1)
139       I1=I1+1
      IF(I1.GT.NM1) GOTO 138
      H(I,1)=H(I,1)-A(I,I1)*H(I1,1)
      GOTO 139
140       CONTINUE
      DO 141 L=1,NM1
      H(L,2)=H(L,1)
141       H(L,1)=0.D0
C       WRITE(3,'(I3,4X,8E13.7)') J,(H(L,2),L=1,4)
C
        IF(J.LE.NX) HY_B2(1,J)=H0
        DO 3333 IZZ=2,NM1+1
        IF(J.LE.NX) THEN
c   en nT
        HY_B2(IZZ,J)=.5D0*(H(IZZ-1,2)+HPLUS(IZZ))
c en microV/m/nT
        EX_B2(IZZ-1,J)=-1.D-11*RO2(IZZ-1,J)/DZI(IZZ-1)*1.D2*
     @               (HY_B2(IZZ,J)-HY_B2(IZZ-1,J))/4.D-4/PI 
        ENDIF
        HPLUS(IZZ)=H(IZZ-1,2)
 3333       CONTINUE
C
      J=J-1
      IF(J.EQ.0) GOTO 400
      IF(J.NE.NX) GOTO 152
      KECS=KECS-2*N*NM1
      GOTO 160
152       IF(J.NE.1) GOTO 154
      KECS=1
      GOTO  160
154       CONTINUE
      KECS=KECS-4*N*NM1
160       CONTINUE
C
      DO 161 I=1,NM1
      DO 159 JV=1,NM1
C       READEC(A(I,JV),KECS,2)
      AR=BUFF(KECS)
      AI=BUFF(KECS+1)
      A(I,JV)=DCMPLX(AR,AI)
159       KECS=KECS+2
161       CONTINUE
      M=2*NM1
      KEC=KECS
C       CALL READEC(B(1),KECS,M)
      DO 9996 KIJ=1,NM1
      AR=BUFF(KEC+KIJ-1)
      AI=BUFF(KEC+KIJ)
      B(KIJ)=DCMPLX(AR,AI)
9996  KEC=KEC+1
      KECS=KECS+M
C
C
      I=N
162       I=I-1
      IF(I.EQ.0) GOTO 140
C
      IF(I.EQ.NM1) GOTO 170
C
      I1=I+1
      I2=1
      H(I,1)=B(I)-A(I,I2)*H(I1,1)
163       IF(I1.EQ.NM1) GOTO 164
      I1=I1+1
      I2=I2+1
      H(I,1)=H(I,1)-A(I,I2)*H(I1,1)
      GOTO 163
164       I1=0
165       I2=I2+1
      I1=I1+1
      IF(I1.GT.I) GOTO 162
      H(I,1)=H(I,1)-A(I,I2)*H(I1,2)
      GOTO 165
170       I1=1
      I2=I1
      H(I,1)=B(I)-A(I,I1)*H(I1,2)
      GOTO 165
400       CONTINUE
C CALCUL DE EZ
C
      COEF=1.D-5/4./PI
      DO 1777 I=1,NZ
      DO 1777 J=2,NX
      IF(I.EQ.1) THEN
      EZ_B2(1,J)=0.D0
      ELSE
      X=DXJ(J-1)/2.
      AA=2.*(HY_B2(I,J)-HY_B2(I,J-1))/(DXJ(J)+DXJ(J-1))
      BB=HY_B2(I,J-1)
      HYI2=AA*X+BB
            IF(J.EQ.2) THEN
            X=-DXJ(1)/2.
            HYI1=AA*X+BB
            ENDIF
      EZ_B2(I,J-1)=RO2(I,J-1)*(HYI2-HYI1)/DXJ(J-1)*COEF
      HYI1=HYI2
            IF(J.EQ.NX) THEN
            X=DXJ(NX-1)/2.+DXJ(NX)
            HYI2=AA*X+BB
            EZ_B2(I,NX)=RO2(I,NX)*(HYI2-HYI1)/DXJ(NX)*COEF
            ENDIF
      ENDIF
 1777       CONTINUE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 3000 I=1,N
 3000       DZI(I)=DZI(I)*1.D-2
      DO 3001 J=1,NX
 3001       DXJ(J)=DXJ(J)*1.D-2
 
      DO 3002 I=1,NZ
        DO 3002 J=1,NX
C           Transfer the E-Field to V/m
            EX_B2(I,J) = EX_B2(I,J)/1.D6
            EZ_B2(I,J) = EZ_B2(I,J)/1.D6
C           Transfer the B to a H field and nT to T
            HY_B2(I,J) = (HY_B2(I,J))/(4*PI*1.D2)

            HY_REAL((J-1)*NZ + I) = DREAL(HY_B2(I,J))
            HY_IMAG((J-1)*NZ + I) = DIMAG(HY_B2(I,J))
            
            EX_REAL((J-1)*NZ + I) = DREAL(EX_B2(I,J))
            EX_IMAG((J-1)*NZ + I) = DIMAG(EX_B2(I,J))
            
            EZ_REAL((J-1)*NZ + I) = DREAL(EZ_B2(I,J))
            EZ_IMAG((J-1)*NZ + I) = DIMAG(EZ_B2(I,J))
            
3002  CONTINUE 
 
      deallocate(RO2,stat=dealloc_status)
      deallocate(HY_B2,stat=dealloc_status)
      deallocate(EX_B2,stat=dealloc_status)
      deallocate(EZ_B2,stat=dealloc_status)
      deallocate(H,stat=dealloc_status)
      
      deallocate(A,stat=dealloc_status)
      deallocate(B,stat=dealloc_status)
      deallocate(AIT,stat=dealloc_status)
      deallocate(BI,stat=dealloc_status)
      deallocate(EXX,stat=dealloc_status)
      deallocate(EXX1,stat=dealloc_status)
      
      deallocate(DXM,stat=dealloc_status)
      deallocate(DZM,stat=dealloc_status)
      
      deallocate(BUFF,stat=dealloc_status)
      deallocate(HPLUS,stat=dealloc_status)
CCCCCCCCCCCCCCCCCCCCCCCCCCC

      RETURN
      END
C
C
