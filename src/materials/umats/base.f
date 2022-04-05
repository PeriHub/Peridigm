      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
C    INCLUDE 'ABA_PARAM.INC'
      implicit real(8) (a-h,o-z)
C
      CHARACTER*80 CMNAME
      DOUBLE PRECISION, DIMENSION(NTENS,NTENS) :: C
      DOUBLE PRECISION, DIMENSION(NTENS) :: STRANNP1
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),
     2 DDSDDT(NTENS),DRPLDE(NTENS),
     3 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     5 JSTEP(4)
C
C  EVALUATE NEW STRESS TENSOR
C    
      DO K1=1,NTENS
            STRANNP1(K1) = STRAN(K1) + DSTRAN(K1)
            DO K2=1,NTENS
                  C(K1,K2) = 0.0
            END DO
      END DO
      IF (NTENS == 6) THEN
            C(1,1) = PROPS(1)-PROPS(2)*PROPS(2)/PROPS(1)
            C(2,2) = PROPS(1)-PROPS(2)*PROPS(2)/PROPS(1)
            C(3,3) = PROPS(1)-PROPS(2)*PROPS(2)/PROPS(1)
            C(1,2) = PROPS(2)-PROPS(2)*PROPS(2)/PROPS(1)
            C(2,1) = PROPS(2)-PROPS(2)*PROPS(2)/PROPS(1)
            C(3,1) = PROPS(2)-PROPS(2)*PROPS(2)/PROPS(1)
            C(1,3) = PROPS(2)-PROPS(2)*PROPS(2)/PROPS(1)
            C(2,3) = PROPS(2)-PROPS(2)*PROPS(2)/PROPS(1)
            C(3,2) = PROPS(2)-PROPS(2)*PROPS(2)/PROPS(1)
            C(4,4) = PROPS(3)
            C(5,5) = PROPS(3)
            C(6,6) = PROPS(3)
      ELSE IF (NTENS == 4) THEN
            C(1,1) = PROPS(1)-PROPS(2)*PROPS(2)/PROPS(1)
            C(2,2) = PROPS(1)-PROPS(2)*PROPS(2)/PROPS(1)
            C(3,3) = PROPS(1)-PROPS(2)*PROPS(2)/PROPS(1)
            C(1,2) = PROPS(2)-PROPS(2)*PROPS(2)/PROPS(1)
            C(2,1) = PROPS(2)-PROPS(2)*PROPS(2)/PROPS(1)
            C(3,1) = PROPS(2)-PROPS(2)*PROPS(2)/PROPS(1)
            C(1,3) = PROPS(2)-PROPS(2)*PROPS(2)/PROPS(1)
            C(2,3) = PROPS(2)-PROPS(2)*PROPS(2)/PROPS(1)
            C(3,2) = PROPS(2)-PROPS(2)*PROPS(2)/PROPS(1)
            C(4,4) = PROPS(3)
      ELSE
            C(1,1) = PROPS(1)-PROPS(2)*PROPS(2)/PROPS(1)
            C(2,2) = PROPS(1)-PROPS(2)*PROPS(2)/PROPS(1)
            C(1,2) = PROPS(2)-PROPS(2)*PROPS(2)/PROPS(1)
            C(2,1) = PROPS(2)-PROPS(2)*PROPS(2)/PROPS(1)
            C(3,3) = PROPS(3)

      END IF      

      STRESS = MATMUL(C, STRANNP1)
      RETURN
      END