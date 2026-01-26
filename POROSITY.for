C=========================================================================
C 孔隙率计算文件，主要包含下列子程序：
C POROSLPINT               功能为计算有效孔隙率并据此对于滑移系强度进行修正
C POROHARD                 功能为计算有效孔隙率并据此对于硬化刚度矩阵进行修正
C CALC_DPORDG              功能为计算孔隙率对滑移的导数
C CALC_DPLMDG              功能为计算平均体积塑性应变对滑移导数
C CALC_WEAK_FUNCTION       功能为计算软化函数
C CALC_EFFECTIVE_POROSITY  功能为计算有效孔隙率
C=========================================================================
C
C 非局部化模型说明：
C   为避免孔隙率演化导致的应变局部化问题，形核项中的等效塑性应变
C   使用晶粒平均值（非局部值）EQVPL_NL，而非局部值EQVPL。
C   这相当于隐式梯度非局部模型在晶粒尺度上的离散化。
C
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
C  注意：形核项使用非局部等效塑性应变EQVPL_NL（晶粒平均值）
C        生长项继续使用局部值
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

        DO I = 1, NSLPTL
          X = TAUSLP(I)/EQVSTR
          Y = HYSTR/EQVSTR
          IF (X.GE.-1.0D3.AND.X.LE.1.0D3) THEN
            X = X
          ELSE
            X = 0.0D0
          END IF
          IF (Y.GE.-1.0D3.AND.Y.LE.1.0D3) THEN
            Y = Y
          ELSE
            Y = 0.0D0
          END IF
C         形核项：使用局部EQVPL（此处保持不变，非局部化在PORO_EVOL中处理）
          DFDGN = INCLUSF/(DSQRT(2.0D0*PI)*SEVLEPN) *
     &    DEXP(-0.5D0*((EQVPL-EVALEPN)/SEVLEPN)**2.0D0) * X
C         生长项：使用局部值
          DFDGE = (PORO_B * DSINH((PORO_N-0.5D0)/(PORO_N+0.5D0)*Y) *
     &    ((1.0D0-POROS)**(-PORO_N)-(1.0D0-POROS))) * DSIGN(1.0D0, X)
          DFDG(I) = DFDGN + DFDGE
        END DO
      RETURN
      END SUBROUTINE CALC_DPORDG


C-POROEVOL 功能为对孔隙率与等效塑性应变进行迭代
C  修改：增加EQVPL_NL参数，形核项使用非局部等效塑性应变
C  输入：
C    DGAMMA   - 滑移增量
C    POROS    - 当前孔隙率（输入输出）
C    EQVPL    - 局部等效塑性应变（输入输出）
C    EQVPL_NL - 非局部等效塑性应变（晶粒平均值，仅输入）
C    EQVSTR   - 等效应力
C    TAUSLP   - 分切应力
C    DPORDG   - 孔隙率对滑移的导数
C
      SUBROUTINE PORO_EVOL(DGAMMA,POROS,EQVPL,EQVPL_NL,EQVSTR,
     2  TAUSLP,DPORDG)
        USE PARABANK
        IMPLICIT NONE
        INTEGER I, J
        REAL(KIND=8)::DGAMMA(ND), DPOROS, POROS, EQVPL, EQVPL_NL
        REAL(KIND=8)::DPORDG(ND)
        REAL(KIND=8)::EQVSTR, DEQVPL, TAUSLP(ND)
        REAL(KIND=8)::W, X
        REAL(KIND=8)::DFDGN, DFDGE
        real(KIND=8)::PI = 3.1415926535D0

        DPOROS = 0.0D0
        DEQVPL = 0.0D0

        DO I = 1, NSLPTL
          X = TAUSLP(I)/EQVSTR
          IF (X.GE.-1.0D3.AND.X.LE.1.0D3) THEN
            X = X
          ELSE
            X = 0.0D0
          END IF

C         ============ 非局部化核心修改 ============
C         形核项：使用非局部等效塑性应变 EQVPL_NL（晶粒平均值）
C         这是隐式梯度非局部模型在晶粒尺度上的离散化
          DFDGN = INCLUSF/(DSQRT(2.0D0*PI)*SEVLEPN) *
     &    DEXP(-0.5D0*((EQVPL_NL-EVALEPN)/SEVLEPN)**2.0D0) * X
C         生长项：直接使用DPORDG中已计算的值
C         ==========================================

C         累加孔隙率增量和等效塑性应变增量
          DPOROS = DPOROS + DGAMMA(I) * (DFDGN + (DPORDG(I) -
     &      INCLUSF/(DSQRT(2.0D0*PI)*SEVLEPN) *
     &      DEXP(-0.5D0*((EQVPL-EVALEPN)/SEVLEPN)**2.0D0) * X))
          DEQVPL = DEQVPL + DGAMMA(I) * X
        END DO

C       更新局部等效塑性应变和孔隙率
        EQVPL = EQVPL + DEQVPL
        POROS = POROS + DPOROS

C       孔隙率限制（防止数值问题）
        IF (POROS .LT. 0.0D0) THEN
          POROS = 0.0D0
        END IF
        IF (POROS .GT. 0.9D0) THEN
          POROS = 0.9D0
        END IF

      RETURN
      END SUBROUTINE PORO_EVOL


C-CALC_WEAK_FUNCTION 功能为计算软化函数
      SUBROUTINE CALC_WEAK_FUNCTION(W, POROS, PORO_A)
        IMPLICIT NONE
        REAL(KIND=8)::POROS, W, PORO_A
        W = (1-DTANH(PORO_A*POROS))
      RETURN
      END SUBROUTINE CALC_WEAK_FUNCTION

