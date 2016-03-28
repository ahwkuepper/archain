
        PROGRAM RECOIL


        REAL*8 MBH1, MBH2, Q12, ETAQ
        REAL*8 VBH1(3), VBH2(3), E1(3), E2(3), NORM(3)
        REAL*8 VBHABS1, VBHABS2, NORMABS
        REAL*8 AMC, BMC, HC, XIC, VC11, VAC, VBC, VCC
        REAL*8 VMASS, VPERP, VPAR
        REAL*8 ALPHA1(3), ALPHA2(3), SBH1(3), SBH2(3)
        REAL*8 SABSBH1, SABSBH2, ARAND1, ARAND2
        REAL*8 INVEC(3), INVECABS, DELTAVEC(3), DELTAVECABS
        REAL*8 PHIDELTA, STILDE(3), ALPHAPERP
        REAL*8 VKICK(3)

        REAL*8 Clight, GAUSS, G

        Clight = 300000.0
        G = 0.0043

C       NUMERICAL COEFFICIENTS
        AMC = 1.2e4        !km/s
        BMC = -0.93
        HC = 6.9e3         !km/s
        XIC = 2.5307274154 !145.0 deg in rad
        VC11 = 3677.76     !km/s
        VAC = 2481.21      !km/s
        VBC = 1792.45      !km/s
        VCC = 1506.52      !km/s

        DO I=1,1000000

C       INPUT PARAMETERS
        MBH1 = RAND(0)*1.e8
        MBH2 = RAND(0)*1.e8
        CALL GETGAUSS(GAUSS)
        VBH1(1) = GAUSS*1000.0
        CALL GETGAUSS(GAUSS)
        VBH1(2) = GAUSS*1000.0
        CALL GETGAUSS(GAUSS)
        VBH1(3) = GAUSS*1000.0
        CALL GETGAUSS(GAUSS)
        VBH2(1) = GAUSS*1000.0
        CALL GETGAUSS(GAUSS)
        VBH2(2) = GAUSS*1000.0
        CALL GETGAUSS(GAUSS)
        VBH2(3) = GAUSS*1000.0


        CALL GETGAUSS(GAUSS)
        SBH1(1) = GAUSS !E1 spin component -> S = a*G*MBH**2/c with a E {0,1}
        CALL GETGAUSS(GAUSS)
        SBH1(2) = GAUSS !E2 spin component
        CALL GETGAUSS(GAUSS)
        SBH1(3) = GAUSS !NORM spin component
        SABSBH1 = sqrt(SBH1(1)**2+SBH1(2)**2+SBH1(3)**2)
        ARAND1 = RAND(0)

        SBH1(1) = ARAND1*G*MBH1**2/Clight*SBH1(1)/SABSBH1
        SBH1(2) = ARAND1*G*MBH1**2/Clight*SBH1(2)/SABSBH1
        SBH1(3) = ARAND1*G*MBH1**2/Clight*SBH1(3)/SABSBH1

        CALL GETGAUSS(GAUSS)
        SBH2(1) = GAUSS !E1 spin component
        CALL GETGAUSS(GAUSS)
        SBH2(2) = GAUSS !E2 spin component
        CALL GETGAUSS(GAUSS)
        SBH2(3) = GAUSS !NORM spin component
        SABSBH2 = sqrt(SBH2(1)**2+SBH2(2)**2+SBH2(3)**2)
        ARAND2 = RAND(0)

        SBH2(1) = ARAND2*G*MBH2**2/Clight*SBH2(1)/SABSBH2
        SBH2(2) = ARAND2*G*MBH2**2/Clight*SBH2(2)/SABSBH2
        SBH2(3) = ARAND2*G*MBH2**2/Clight*SBH2(3)/SABSBH2

        CALL GETGAUSS(GAUSS)
        INVEC(1) = GAUSS  !arbitrary infall vector
        CALL GETGAUSS(GAUSS)
        INVEC(2) = GAUSS
        CALL GETGAUSS(GAUSS)
        INVEC(3) = GAUSS
        INVECABS = sqrt(INVEC(1)**2+INVEC(2)**2+INVEC(3)**2)

C       MASS RATIO
        Q12 = MBH1/MBH2
        ETAQ = Q12/(1.0+Q12)**2

C       CONSTRUCT UNIT VECTORS
        VBHABS1 = sqrt(VBH1(1)**2+VBH1(2)**2+VBH1(3)**2)
        VBHABS2 = sqrt(VBH2(1)**2+VBH2(2)**2+VBH2(3)**2)

C       USE VBH1 AS FIRST UNIT VECTOR IN ORBITAL PLANE
        E1(1) = VBH1(1)/VBHABS1
        E1(2) = VBH1(2)/VBHABS1
        E1(3) = VBH1(3)/VBHABS1

C       CONSTRUCT NORMAL VECTOR
        NORM(1) = VBH1(2)*VBH2(3)-VBH1(3)*VBH2(2)
        NORM(2) = VBH1(3)*VBH2(1)-VBH1(1)*VBH2(3)
        NORM(3) = VBH1(1)*VBH2(2)-VBH1(2)*VBH2(1)
        NORMABS = sqrt(NORM(1)**2+NORM(2)**2+NORM(3)**2)
        NORM(1) = NORM(1)/NORMABS
        NORM(2) = NORM(2)/NORMABS
        NORM(3) = NORM(3)/NORMABS

C       SECOND UNIT VECTOR IN ORBITAL PLANE
        E2(1) = NORM(2)*E1(3)-NORM(3)*E1(2)
        E2(2) = NORM(3)*E1(1)-NORM(1)*E1(3)
        E2(3) = NORM(1)*E1(2)-NORM(2)*E1(1)

C        WRITE(*,*) "  E1: ", E1(1), E1(2), E1(3)
C        WRITE(*,*) "  E2: ", E2(1), E2(2), E2(3)
C        WRITE(*,*) "NORM: ", NORM(1), NORM(2), NORM(3)


C       FIRST KICK COMPONENT
        VMASS = AMC*ETAQ**2*(1.0-Q12)/(1.0+Q12)*(1.0+BMC*ETAQ)
C        WRITE(*,*) "VMASS: ", VMASS


C       DIMENSIONLESS SPINS
        ALPHA1(1) = Clight/G*SBH1(1)/MBH1**2
        ALPHA1(2) = Clight/G*SBH1(2)/MBH1**2
        ALPHA1(3) = Clight/G*SBH1(3)/MBH1**2
        ALPHA2(1) = Clight/G*SBH2(1)/MBH2**2
        ALPHA2(2) = Clight/G*SBH2(2)/MBH2**2
        ALPHA2(3) = Clight/G*SBH2(3)/MBH2**2


C       SECOND KICK COMPONENT
        VPERP = HC*ETAQ**2/(1.0+Q12)*(ALPHA2(3)-Q12*ALPHA1(3))
C        WRITE(*,*) "VPERP: ", VPERP


C       CONSTRUCT IN-PLANE COMPONENT
        DELTAVEC(1) = (MBH1+MBH2)*(SBH2(1)/MBH2-SBH1(1)/MBH1)
        DELTAVEC(2) = (MBH1+MBH2)*(SBH2(2)/MBH2-SBH1(2)/MBH1)
        DELTAVEC(3) = 0.0
        DELTAVECABS = sqrt(DELTAVEC(1)**2+DELTAVEC(2)**2)

C       ANGLE BETWEEN INFALL VECTOR AND IN-PLANE COMPONENT
        PHIDELTA = (DELTAVEC(1)*INVEC(1)+DELTAVEC(2)*INVEC(2)+
     &              DELTAVEC(3)*INVEC(3))/(DELTAVECABS*INVECABS)
        PHIDELTA = ACOS(PHIDELTA)

        STILDE(1) = 2.0*(ALPHA2(1)+Q12*ALPHA1(1))/(1.0+Q12)**2
        STILDE(2) = 2.0*(ALPHA2(2)+Q12*ALPHA1(2))/(1.0+Q12)**2
        STILDE(3) = 2.0*(ALPHA2(3)+Q12*ALPHA1(3))/(1.0+Q12)**2

        ALPHAPERP = sqrt((ALPHA2(1)-Q12*ALPHA1(1))**2+
     &                   (ALPHA2(2)-Q12*ALPHA1(2))**2)


C       THIRD KICK COMPONENT
        VPAR = 16.0*ETAQ**2/(1.0+Q12)*(VC11+VAC*STILDE(3)+
     &         VBC*STILDE(3)**2+VCC*STILDE(3)**3)*
     &         ALPHAPERP*cos(PHIDELTA)
C        WRITE(*,*) "VPAR: ", VPAR


        VKICK(1) = VMASS*E1(1)+VPERP*(cos(XIC)*E1(1)+sin(XIC)*E2(1))
     &             +VPAR*NORM(1)
        VKICK(2) = VMASS*E1(2)+VPERP*(cos(XIC)*E1(2)+sin(XIC)*E2(2))
     &             +VPAR*NORM(2)
        VKICK(3) = VMASS*E1(3)+VPERP*(cos(XIC)*E1(3)+sin(XIC)*E2(3))
     &             +VPAR*NORM(3)


        WRITE(*,*) MBH1, VBH1(1), VBH1(2), VBH1(3)
     &           , SBH1(1), SBH1(2), SBH1(3)
     &           , MBH2, VBH2(1), VBH2(2), VBH2(3)
     &           , SBH2(1), SBH2(2), SBH2(3)
     &           , VKICK(1), VKICK(2), VKICK(3)
     &           , VMASS, VPERP, VPAR
     &           , SQRT(VKICK(1)**2+VKICK(2)**2+VKICK(3)**2)

        END DO


        END



************************************************************
************************************************************



        SUBROUTINE GETGAUSS(GAUSS)

        REAL*8 GAUSS
        REAL*8 XX, YY, QQ

        QQ = 2.0

        DO WHILE (QQ.GT.1.0)
            XX = 2.0*RAND(0)-1.0
            YY = 2.0*RAND(0)-1.0
            QQ = XX*XX + YY*YY
        END DO

        GAUSS = SQRT(-2.0*LOG(QQ)/QQ)*XX

        RETURN

        END

