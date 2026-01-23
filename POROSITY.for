C=========================================================================
C 孔隙率计算文件，主要包含下列子程序：
C POROSLPINT               功能为计算有效孔隙率并据此对于滑移系强度进行修正
C POROHARD                 功能为计算有效孔隙率并据此对于硬化刚度矩阵进行修正
C CALC_DPORDG              功能为计算孔隙率对滑移的导数 
C CALC_DPLMDG              功能为计算平均体积塑性应变对滑移导数
C CALC_WEAK_FUNCTION       功能为计算软化函数
C CALC_EFFECTIVE_POROSITY  功能为计算有效孔隙率
C=========================================================================


C-POROSLPINT 功能为计算有效孔隙率并据此对于滑移系强度进行修正
      SUBROUTINE PORO_SLIP_INIT(GSLIP0)
        USE PARABANK
        IMPLICIT NONE
        INTEGER I
        REAL(KIND=8)::GSLIP0(ND), W
        CALL CALC_WEAK_FUNCTION(W, POROS0, PORO_A)
        DO I=1,NSLPTL
          GSLIP0(I) = GSLIP0(I) * W
        END DO
      RETURN  
      END SUBROUTINE PORO_SLIP_INIT


C-POROHARD 功能为计算有效孔隙率并据此对于硬化刚度矩阵进行修正
      SUBROUTINE PORO_HARD(H, POROS, EQVSTR, HYSTR, EQVPL, TAUSLP, 
     2              GSLIP, DPORDG)
        USE PARABANK
        IMPLICIT NONE
        INTEGER I, J
        REAL(KIND=8)::H(ND, ND), POROS, DPORDG(ND)
        REAL(KIND=8)::HYSTR, EQVSTR, EQVPL, TAUSLP(ND), GSLIP(ND)
        real(KIND=8)::PI=3.1415926535D0
        real(KIND=8)::W, L
        
        CALL CALC_WEAK_FUNCTION(W, POROS, PORO_A)
        CALL CALC_DPORDG(DPORDG, POROS, EQVSTR, HYSTR, EQVPL, TAUSLP)
        DO I=1,NSLPTL
          DO J=1,NSLPTL
            L=GSLIP(I)/W*PORO_A*DCOSH(PORO_A*POROS)**(-2.0D0)*DPORDG(J)
            H(I,J)=H(I,J)*W-L
          END DO
        END DO
      RETURN  
      END SUBROUTINE PORO_HARD


C-CALC_DPORDG 功能为计算孔隙率对滑移的导数       
      SUBROUTINE CALC_DPORDG(DFDG, POROS, EQVSTR, HYSTR, EQVPL, TAUSLP)
        USE PARABANK
        IMPLICIT NONE
        INTEGER I, J
        REAL(KIND=8)::DFDG(ND), POROS
        REAL(KIND=8)::X, Y !(滑移应力与等效塑性应变的比值)
        REAL(KIND=8)::DFDGN !(孔隙率形核部分对等效塑性应变的导数)
        REAL(KIND=8)::DFDGE !(孔隙率生长部分对等效塑性应变的导数)
        REAL(KIND=8)::HYSTR, EQVSTR, EQVPL, TAUSLP(ND)
        real(KIND=8)::PI = 3.1415926535D0               

        DO I = 1, ND
          X = TAUSLP(I)/EQVSTR
          Y = HYSTR/EQVSTR
          IF (X.GE.-1.0D3.AND.X.LE.1.0D3) THEN !(2025/05/26 lbt:如果X为NaN，则将其置为0)
            X = X
          ELSE
            X = 0.0D0
          END IF
          IF (Y.GE.-1.0D3.AND.Y.LE.1.0D3) THEN !(2025/05/26 lbt:如果Y为NaN，则将其置为0)
            Y = Y
          ELSE
            Y = 0.0D0
          END IF
          DFDGN = INCLUSF/(DSQRT(2.0D0*PI)*SEVLEPN) *
     &    DEXP(-0.5D0*((EQVPL-EVALEPN)/SEVLEPN)**2.0D0) * X
          DFDGE = (PORO_B * DSINH((PORO_N-0.5D0)/(PORO_N+0.5D0)*Y) *
     &    ((1.0D0-POROS)**(-PORO_N)-(1.0D0-POROS))) * DSIGN(1.0D0, X)
          DFDG(I) = DFDGN + DFDGE
        END DO
      RETURN  
      END SUBROUTINE CALC_DPORDG


C-POROEVOL 功能为对孔隙率与等效塑性应变进行迭代        
      SUBROUTINE PORO_EVOL(DGAMMA, POROS, EQVPL, EQVSTR, TAUSLP, DPORDG)
        USE PARABANK
        IMPLICIT NONE
        INTEGER I, J
        REAL(KIND=8)::DGAMMA(ND), DPOROS, POROS, EQVPL, DPORDG(ND)
        REAL(KIND=8)::EQVSTR, DEQVPL, TAUSLP(ND)
        REAL(KIND=8)::W, X
        DPOROS = 0.0D0
        DEQVPL = 0.0D0
        DO I = 1, ND
          X = TAUSLP(I)/EQVSTR
          IF (X.GE.-1.0D3.AND.X.LE.1.0D3) THEN !(2025/05/26 lbt:如果X为NaN，则将其置为0)
            X = X
          ELSE
            X = 0.0D0
          END IF    
          DPOROS = DPOROS + DGAMMA(I) * DPORDG(I)
          DEQVPL = DEQVPL + DGAMMA(I) * X
        END DO
        EQVPL = EQVPL + DEQVPL
        POROS = POROS + DPOROS
      RETURN
      END SUBROUTINE PORO_EVOL


C-CALC_WEAK_FUNCTION 功能为计算软化函数
      SUBROUTINE CALC_WEAK_FUNCTION(W, POROS, PORO_A)
        IMPLICIT NONE
        REAL(KIND=8)::POROS, W, PORO_A
        W = (1-DTANH(PORO_A*POROS))
      RETURN  
      END SUBROUTINE CALC_WEAK_FUNCTION




