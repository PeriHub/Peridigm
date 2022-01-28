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
      DOUBLE PRECISION, DIMENSION(6,6) :: C
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),
     2 DDSDDT(NTENS),DRPLDE(NTENS),
     3 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     5 JSTEP(4)
C
C  EVALUATE NEW STRESS TENSOR
C     
      C(1,1) = PROPS(1)-PROPS(3)*PROPS(3)/PROPS(7)
      C(1,2) = PROPS(2)-PROPS(3)*PROPS(8)/PROPS(7)
      C(1,3) = 0.0
      C(1,4) = 0.0
      C(1,5) = 0.0
      C(1,6) = 0.0
      C(2,1) = PROPS(2)-PROPS(3)*PROPS(8)/PROPS(7)
      C(2,2) = PROPS(7)-PROPS(3)*PROPS(8)/PROPS(7)
      C(2,3) = 0.0
      C(2,4) = 0.0
      C(2,5) = 0.0
      C(2,6) = 0.0
      C(3,1) = 0.0
      C(3,2) = 0.0
      C(3,3) = 0.0
      C(3,4) = 0.0
      C(3,5) = 0.0
      C(3,6) = 0.0
      C(4,1) = 0.0
      C(4,2) = 0.0
      C(4,3) = 0.0
      C(4,4) = 0.0
      C(4,5) = 0.0
      C(4,6) = 0.0
      C(5,1) = 0.0
      C(5,2) = 0.0
      C(5,3) = 0.0
      C(5,4) = 0.0
      C(5,5) = 0.0
      C(5,6) = 0.0
      C(6,1) = 0.0
      C(6,2) = 0.0
      C(6,3) = 0.0
      C(6,4) = 0.0
      C(6,5) = 0.0
      C(6,6) = PROPS(21)
      STRESS(1)= (STRAN(2)+STRAN(4))*C(1,6) + (STRAN(3)+STRAN(5))*C(1,5) + (STRAN(6)+STRAN(4))*C(1,4) + C(1,1)*STRAN(1) + C(1,2)*STRAN(5) + C(1,3)*STRAN(3)
      STRESS(2)= (STRAN(2)+STRAN(4))*C(6,6) + (STRAN(3)+STRAN(5))*C(6,5) + (STRAN(6)+STRAN(4))*C(6,4) + C(6,1)*STRAN(1) + C(6,2)*STRAN(5) + C(6,3)*STRAN(3)
      STRESS(3)= (STRAN(2)+STRAN(4))*C(5,6) + (STRAN(3)+STRAN(5))*C(5,5) + (STRAN(6)+STRAN(4))*C(5,4) + C(5,1)*STRAN(1) + C(5,2)*STRAN(5) + C(5,3)*STRAN(3)
      STRESS(4)= (STRAN(2)+STRAN(4))*C(2,6) + (STRAN(3)+STRAN(5))*C(2,5) + (STRAN(6)+STRAN(4))*C(2,4) + C(2,1)*STRAN(1) + C(2,2)*STRAN(5) + C(2,3)*STRAN(3)
      STRESS(5)= (STRAN(2)+STRAN(4))*C(4,6) + (STRAN(3)+STRAN(5))*C(4,5) + (STRAN(6)+STRAN(4))*C(4,4) + C(4,1)*STRAN(1) + C(4,2)*STRAN(5) + C(4,3)*STRAN(3)
      STRESS(6)= (STRAN(2)+STRAN(4))*C(3,6) + (STRAN(3)+STRAN(5))*C(3,5) + (STRAN(6)+STRAN(4))*C(3,4) + C(3,1)*STRAN(1) + C(3,2)*STRAN(5) + C(3,3)*STRAN(3)
      RETURN
      END