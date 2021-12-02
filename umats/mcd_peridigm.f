C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %                                                   PERIDIGM VERSION                                                    %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     MCODAC - Modular COmposite Damage Analysis Code
C     Peridigm Interface
C     Created on 04.07.2019	  	  	
C     
C     Version: 1.0 
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     Requirements and dependencies:
C        - 
C
C     Changelog:
C        - 		   							    	   
C
C     @author: garb_ma            [DLR-FA,STM Braunschweig]
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  	

C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
C     UPERMAT
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		        
SUBROUTINE UPERMAT(STRESSNEW, STATENEW,
    &           STRESSOLD, STATEOLD, STRAININC, PROPS, 
    &           STEPTIME, TOTALTIME, DT, NodeID, 
    &	         MATID, NPOINT, NDIR, NSHR, NSTATEV, NPROPS)
    &           BIND(C, NAME="UPERMAT")
     USE MCODAC 
     IMPLICIT NONE
C
C      Externals
     PROCEDURE(PERIDIGM_UPERMAT_INTERFACE), POINTER :: MCD_PERI_MAT 	
C
C      Input variables (Read) 
     INTEGER(i4), INTENT(IN) :: MATID, NodeID
     INTEGER(i4), INTENT(IN) :: NPOINT, NDIR, NSHR, NSTATEV, NPROPS	  
C	  
     REAL(DP), INTENT(IN) :: STEPTIME, TOTALTIME, DT 
C	  
     REAL(DP), INTENT(IN) :: PROPS(NPROPS)
C	  
     REAL(DP), INTENT(IN) :: STRESSOLD(NDIR+NSHR)
     REAL(DP), INTENT(IN) :: STATEOLD(NSTATEV)
     REAL(DP), INTENT(IN) :: STRAININC(NDIR+NSHR) 
C	  	  
C      Return variables (Write)
     REAL(DP), INTENT(OUT) :: STRESSNEW(NDIR+NSHR) 
     REAL(DP), INTENT(OUT) :: STATENEW(NSTATEV)	
C     ======================================================================================
C     UPERMAT directory:
C              (1) Sensor: 					Key = 1001
C              (2) Lamina:  				Key = 1002 
C              (3) Shokrieh: 				Key = 1003
C              (4) Fiber: 					Key = 1004
C              (5) Resin: 					Key = 1005
C              (6) NeoHookean (from Muesli):	Key = 1006
C              (7) Hail:					Key = 1007
C     ======================================================================================
C      Procedure 
     IF (ASSOCIATED(MCD_PERI_MAT)) NULLIFY(MCD_PERI_MAT)
C	  
C	  Use procedure pointer to access material sharing the same standard interface.
     IF (MATID .EQ. 1002) MCD_PERI_MAT => UPERMAT_Lamina  
     IF (MATID .EQ. 1006) MCD_PERI_MAT => UPERMAT_NeoHookean
C	  
C	  ! Call selected material
     IF (ASSOCIATED(MCD_PERI_MAT)) THEN 
        CALL MCD_PERI_MAT(STRESSNEW, STATENEW,
    &        STRESSOLD, STATEOLD, STRAININC, PROPS, 
    &        STEPTIME, TOTALTIME, DT, 	 
    &	      MATID, NPOINT, NDIR, NSHR, NSTATEV, NPROPS)		  
     END IF	  
C
C      Dealloacte local function pointer with no further use	 
     IF (ASSOCIATED(MCD_PERI_MAT)) NULLIFY(MCD_PERI_MAT)	 	
C	 
     RETURN
     END SUBROUTINE UPERMAT
