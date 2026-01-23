C-----------------------------------------------------------------------------------------
      SUBROUTINE DIS_IRRA_INIT(NOEL,NMD,NVOID,
     2 NLOOP,NPRE,EMFP,SVOID,SLOOP,SPRE,
     3 TOTMD,TOTVOID,TOTLOOP,TOTPRE,TOTEMFP,SLPDIR,SLPNOR,
     4 HVOID,HLOOP,HPRE)
      USE PARABANK
      IMPLICIT NONE
      INTEGER::NOEL
      DOUBLE PRECISION::NMD(ND), NVOID(ND), NLOOP(ND), NPRE(ND)
      DOUBLE PRECISION::EMFP(ND), FDD(ND)
      DOUBLE PRECISION::TOTMD, TOTVOID, TOTLOOP, TOTPRE, TOTEMFP
      DOUBLE PRECISION::SVOID, SLOOP, SPRE            !每个缺陷的
      DOUBLE PRECISION::HVOID, HLOOP, HPRE      !每个缺陷的硬化系数
      DOUBLE PRECISION::SLPDIR(3,ND), SLPNOR(3,ND)  
      REAL*8 TEMP    
      INTEGER I,J,ID

!     位错信息的赋值
      ID = 0
      TOTMD=0.0D0
      DO I=1, NSET
        DO J=1, NSLIP(I)
          ID = ID+1
          NMD(ID)=DISLOCATION_DENSITY(I)
          TOTMD=TOTMD+NMD(ID)
        END DO
      END DO

!     赋值空隙密度，环状缺陷密度，预沉淀密度
      IF (IRADON.EQ.1) THEN
        DO I=1,NSLPTL
          NVOID(I)=0.0D0
          NLOOP(I)=0.0D0
          NPRE(I)=0.0D0
          TOTVOID=TOTVOID+NVOID(I)
          TOTLOOP=TOTLOOP+NLOOP(I)
          TOTPRE=TOTPRE+NPRE(I)
        END DO
      END IF

!     计算林位错密度
!     初始化
      DO I=1,ND
        FDD(I)=0.0D0
      END DO
!     计算
      DO I=1,NSLPTL
        FDD(I)=0.0D0
        DO J=1,NSLPTL
          FDD(I)=FDD(I)+NMD(J)*dABS(SLPNOR(1,I)*SLPNOR(1,J)+
     2    SLPNOR(2,I)*SLPNOR(2,J)+SLPNOR(3,I)*SLPNOR(3,J))
        END DO
      END DO

!     计算平均自由程
!     初始化
      TOTEMFP=0.0D0
      DO I=1,ND
        EMFP(I)=0.0D0
      END DO
      IF (IRADON.EQ.1) THEN
        DO I=1,NSLPTL
            TEMP=1.0D0*DSQRT(FDD(I))+DSQRT(NVOID(I)*SVOID)+DSQRT(NPRE(I)*SPRE)+DGRAIN !单位mm
            EMFP(I)=1.0D0/TEMP !单位mm
            TOTEMFP=TOTEMFP+EMFP(I)
        END DO
      ELSE
        DO I=1,NSLPTL
            TEMP=1.0D0*DSQRT(FDD(I))+DGRAIN !单位mm
            EMFP(I)=1.0D0/TEMP !单位mm
            TOTEMFP=TOTEMFP+EMFP(I)
        END DO
      END IF
      RETURN
      END SUBROUTINE DIS_IRRA_INIT
C-----------------------------------------------------------------------------------------      
      
C-----------------------------------------------------------------------------------------      
      SUBROUTINE DIS_GSLP_INIT(GSLIP0,NMD,NVOID,NLOOP,
     2  NPRE,SVOID,SLOOP,SPRE)
      USE PARABANK
      IMPLICIT NONE
      DOUBLE PRECISION::GSLIP0(ND), ARHO(ND)
      DOUBLE PRECISION::NMD(ND), NVOID(ND), NLOOP(ND), NPRE(ND)
      DOUBLE PRECISION::SVOID, SLOOP, SPRE     
      DOUBLE PRECISION::EMFP(ND)
      INTEGER I,J
!     首先统计每个滑移系上的位错密度导致的交互硬化总和
      DO I=1,NSLPTL
        ARHO(I) = 0.0D0 !ARHO即为硬化相和密度的乘积开根号
        DO J=1,NSLPTL
          ARHO(I)=ARHO(I)+NMD(J)*HARDM(I,J)
        END DO
        ARHO(I)=DSQRT(ARHO(I))
      END DO

      DO I=NSLPTL+1,ND
        ARHO(I)=0.0D0
      END DO

!     计算各个滑移系的临界剪应力
      IF (IRADON.EQ.1) THEN
!!    子程序使用DBH模型计算滑移系的滑移阻力, 之后详细书写
!!    子程序使用DBH模型计算滑移系的滑移阻力, 之后详细书写
!!    子程序使用DBH模型计算滑移系的滑移阻力, 之后详细书写
      ELSE
        DO I=1,NSLPTL
          GSLIP0(I)=MPNS*SIGMA0(I)+GMOD*BURGER(I)*HRHO*ARHO(I)
        END DO
      END IF

!     通过将某个滑移系的临界剪切应力设置为一个非常大的数值，来关闭该滑移系(用于层状组织的BOR取向关系)
      IF (LAYERSTRUCTURE.EQ.1) THEN
        DO I=1,NSLPTL
          IF (I.EQ.1 .OR. I.EQ.2 .OR. I.EQ.4 .OR. I.EQ.5) THEN
            GSLIP0(I)=GSLIP0(I)
          ELSE
            GSLIP0(I)=GSLIP0(I)+50.0D0
          END IF
        END DO
      END IF

      RETURN      
      END SUBROUTINE DIS_GSLP_INIT
C-----------------------------------------------------------------------------------------

C-----------------------------------------------------------------------------------------
      SUBROUTINE TWIN_TAU_INIT(TAU0)
      USE PARABANK
      IMPLICIT NONE
      DOUBLE PRECISION::TAU0(NTWTL)
      INTEGER I,J,N,TWIN_ID
C    给每个孪晶系赋予初始切变阻力
      TWIN_ID = 0
      DO I=1, TWSET
        DO J=1, NTWIN(I)
          TWIN_ID = TWIN_ID+1
          TAU0(TWIN_ID) = TAU0_TWIN(I)
        END DO
      END DO
      RETURN
      END SUBROUTINE TWIN_TAU_INIT
C-----------------------------------------------------------------------------------------    

C-----------------------------------------------------------------------------------------
      SUBROUTINE DIS_IRRA_HARD(H,NOEL,NMD,NVOID,NLOOP,
     2  NPRE,SVOID,SLOOP,SPRE,DMDDG,DLOOPDG,EMFP)
      USE PARABANK
      IMPLICIT NONE
      INTEGER::NOEL
      DOUBLE PRECISION::ARHO(ND),H(ND,ND)
      DOUBLE PRECISION::NMD(ND), NVOID(ND), NLOOP(ND), NPRE(ND)
      DOUBLE PRECISION::SVOID, SLOOP, SPRE, EMFP(ND)
      DOUBLE PRECISION::DMDDG(ND), DLOOPDG(ND), P, TEMP
      INTEGER I,J,K
!     首先统计每个滑移系上的位错密度导致的交互硬化总和
!     初始化
      DO I=1,ND
        ARHO(I)=0.0D0
      END DO
      DO I=1,NSLPTL
        DO J=1,NSLPTL
          ARHO(I)=ARHO(I)+NMD(J)*HARDM(I,J) !ARHO即为硬化阵和密度的乘积开根号
        END DO
        ARHO(I)=DSQRT(ARHO(I))
      END DO

!     计算总硬化矩阵
      IF (IRADON.EQ.1) THEN
!!    子程序使用DBH模型计算滑移系的滑移阻力, 之后详细书写
!!    子程序使用DBH模型计算滑移系的滑移阻力, 之后详细书写
!!    子程序使用DBH模型计算滑移系的滑移阻力, 之后详细书写
      ELSE
!     计算中间变量P
        K=1 !K用于控制HBCKUP IS NAN不用输出太多遍
!     初始化
        DO I=1,ND
          DO J=1,ND
            H(I,J)=0.0D0
          END DO
        END DO
!     计算硬化矩阵H
        DO I=1,NSLPTL
          P=ARHO(I)
          DO J=1,NSLPTL
            TEMP=0.5D0*P**(-1.0D0)*GMOD*BURGER(I)*HRHO*HARDM(I,J)*DMDDG(J)                           
            IF (ISNAN(TEMP)) THEN
              IF(K.EQ.1) THEN
                write(*,*) 'HBCKUP IS NAN!'
                PROBLEM=1
                K=0
              END IF
            ELSE
              H(I,J)=TEMP
            END IF
          END DO
        END DO
      END IF

      IF (PROBLEM.EQ.1) THEN
        write(114514, *) 'HRHO',HRHO
        write(114514, *) 'ARHO',ARHO
        write(114514, *) 'DMDDG',DMDDG
      END IF

      RETURN 
      END SUBROUTINE DIS_IRRA_HARD
C-----------------------------------------------------------------------------------------

C-----------------------------------------------------------------------------------------
      SUBROUTINE TWIN_HARD(H,SLIP_SUM,TWIN_SUM)
      USE PARABANK
      IMPLICIT NONE
      INTEGER::I,J
      DOUBLE PRECISION::H(ND,ND),SLIP_SUM,TWIN_SUM
      ! 滑移系对于滑移系的硬化矩阵H(1:NS,1:NS)(由DIS_IRRA_HARD给出，此处不需要修改)
      ! 滑移系对孪晶系的硬化矩阵H(NS+1:ND,1:NS)(认为滑移系对于孪晶系硬化为0)
      DO I=NSLPTL+1,ND
        DO J=1,NSLPTL
          H(I,J)=0.0D0
        END DO
      END DO
      ! 孪晶系对滑移系的硬化矩阵H(1:NS,NS+1:ND)
      DO I=1,NSLPTL
        DO J=NSLPTL+1,ND
          H(I,J)=H_SL_TW*SLIP_SUM**TWIN_D  
        END DO
      END DO
      ! 孪晶系对孪晶系的硬化矩阵H(NS+1:ND,NS+1:ND)
      DO I=NSLPTL+1,ND
        DO J=NSLPTL+1,ND
          H(I,J)=H_TW_TW*TWIN_SUM**TWIN_B
        END DO
      END DO
      IF (PROBLEM.EQ.1) THEN
        write(11451, *) 'TWIN HARDENING MATRIX H:'
        DO I=1,ND
          write(114516, *) H(I,1:ND)
        END DO
      END IF
      RETURN
      END SUBROUTINE TWIN_HARD

C-----------------------------------------------------------------------------------------
      SUBROUTINE DIS_IRRA_EVOLUTION(DMDDG,DLOOPDG,NOEL,EMFP,TOTEMFP,TOTMD,
     2  NMD,DGAMMA,SLPDIR,SLPNOR)
      USE PARABANK
      IMPLICIT NONE
      INTEGER::NOEL
      DOUBLE PRECISION::TEMP, TEMP2, EMFP(ND), TOTEMFP
      DOUBLE PRECISION::DMDDG(ND),DLOOPDG(ND)
      DOUBLE PRECISION::NMD(ND),DNMD(ND),DGAMMA(ND),TOTMD
      DOUBLE PRECISION::SLPDIR(3,ND), SLPNOR(3,ND), FDD(ND)
      INTEGER::I, J
      IF (IRADON.EQ.1) THEN
!!    子程序使用DBH模型计算滑移系的滑移阻力, 之后详细书写
!!    子程序使用DBH模型计算滑移系的滑移阻力, 之后详细书写
!!    子程序使用DBH模型计算滑移系的滑移阻力, 之后详细书写      
      ELSE
        
!     计算林位错密度
        DO I=1,NSLPTL
          FDD(I)=0.0D0
          DO J=1,NSLPTL
            FDD(I)=FDD(I)+NMD(J)*dABS(SLPNOR(1,I)*SLPNOR(1,J)+
     2      SLPNOR(2,I)*SLPNOR(2,J)+SLPNOR(3,I)*SLPNOR(3,J))
          END DO
        END DO

!     计算平均自由程
!     初始化
        TOTEMFP=0.0D0
        DO I=1,ND
          EMFP(I)=0.0D0
        END DO
!     计算 
        DO I=1,NSLPTL
            TEMP=1.0D0*DSQRT(FDD(I))+DGRAIN !单位mm
            EMFP(I)=1.0D0/TEMP !单位mm
            TOTEMFP=TOTEMFP+EMFP(I)
        END DO        
!     计算缺陷对于滑移的偏导DMDDG 
!     初始化
        DO I =1, ND
          DMDDG(I)=0.0D0
          DLOOPDG(I)=0.0D0
          DNMD(I)=0.0D0
        END DO
!     计算
        TOTMD=0.0D0
        DO I=1, NSLPTL
          TEMP2=PMUL/BURGER(I)/EMFP(I)-PDYN*NMD(I)
          IF (ISNAN(TEMP2)) THEN
            write(*,*) 'DMDDG NAN'
            DMDDG(I)=DMDDG(I)
            PROBLEM=1
          ELSE
            IF (TEMP2.LT.0.0D0) THEN
              write(*,*) 'DMDDG LOW THAN 0'
              TEMP2=0.0D0
            END IF
            DMDDG(I)=TEMP2*DSIGN(1.0D0, DGAMMA(I))
          END IF
!     更新位错密度
          DNMD(I) = DMDDG(I) * (DGAMMA(I)) 
          NMD(I) = NMD(I) + DNMD(I)
          TOTMD = TOTMD + NMD(I)
        END DO
        DLOOPDG=0.0D0
!     如果出错，输出有效平均自由程，用于调试
        IF (PROBLEM.EQ.1) THEN
          write(114515,'(A)') 'DMDDG:'
          write(114515,'(100(1X,E8.1))') (DMDDG(I), I=1,ND)
          write(114515,'(A)') 'DGAMMA:'
          write(114515,'(100(1X,E8.1))') (DGAMMA(I), I=1,ND)
          write(114515,'(A)') 'PMUL:'
          write(114515,'(1X,E8.1)') PMUL
          write(114515,'(A)') 'DNMD:'
          write(114515,'(100(1X,E8.1))') (DNMD(I), I=1,ND)
          write(114515,'(A)') 'NMD:'
          write(114515,'(100(1X,E8.1))') (NMD(I), I=1,ND)
          write(114515,'(A)') 'FDD:'
          write(114515,'(100(1X,E8.1))') (FDD(I), I=1,ND)
          write(114515,'(A)') 'EMFP:'
          write(114515,'(100(1X,E8.1))') (EMFP(I), I=1,ND)
        END IF
        IF (PROBLEM.EQ.1) THEN
          stop
        END IF
      END IF
      RETURN
      END SUBROUTINE DIS_IRRA_EVOLUTION
C-----------------------------------------------------------------------------------------
      