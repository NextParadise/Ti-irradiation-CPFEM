!***********************************************************************
! 文件名: TI_IHCPFEM.for
! 功能描述: UMAT子程序 - 钛合金辐照硬化模型
!     
!   作者: 刘博韬
!   单位: SJTU
!   创建日期: 2025-11-27
!   最后修改: 2025-12-2
!   版本: v0.5
!     
!   子程序类型: UMAT (User Material)
!   适用版本: Abaqus 2020 及以上
!   编译器: Intel Fortran Compiler 2019
!     
!   功能说明:
!       本子程序实现了辐照硬化下钛合金的晶体塑性增量理论本构模型
!       - 特性1: 针对Ti合金的晶格进行编写
!     
C***********************************************************************
C     材料参数定义 (PROPS数组):- constant = 50
C***********************************************************************
C     +-----------------------------------------------------------------+
C     | 索引范围                  | 变量说明                            |
C     +-----------------------------------------------------------------+
C     | 弹性常数变量 (1-21)                                             |
C     +-----------------------------------------------------------------+
C     | 1      ~ 2                | 各向同性                            |
C     | 1      ~ 3                | 立方系 C11,C12,C44                  |
C     | 1      ~ 21               | 各向异性                            |
C     +-----------------------------------------------------------------+
C     | 滑移系几何信息 (25-40)                                          |
C     +-----------------------------------------------------------------+
C     | 25     ~ 27               | 局域坐标系第一方向                  |
C     | 28     ~ 30               | 全局坐标系第一方向                  |
C     | 33     ~ 35               | 局域坐标系第二方向                  |
C     | 36     ~ 38               | 全局坐标系第二方向                  |
C     +-----------------------------------------------------------------+
!     
C***********************************************************************
C     状态变量定义 (STATEV数组) - 通用化索引
C***********************************************************************
C     定义: ND = NSLPTL + NTWTL (滑移系总数 + 孪晶系总数)
C     总大小: 42 * ND (当 ND=30 时为 900)
C
C     +-----------------------------------------------------------------+
C     | 索引范围 (通用公式)         | 变量说明                          |
C     +-----------------------------------------------------------------+
C     | 基本滑移系+孪晶系变量 (0 ~ 3*ND)                                |
C     +-----------------------------------------------------------------+
C     | 0*ND + 1   ~  1*ND          | 各滑移系+孪晶系强度(Strength)     |
C     | 1*ND + 1   ~  2*ND          | 各滑移系+孪晶系累积滑移量         |
C     | 2*ND + 1   ~  3*ND          | 各滑移系+孪晶系分切应力(RSS)      |
C     +-----------------------------------------------------------------+
C     | 滑移系+孪晶系几何信息 (3 ~ 9*ND)                                |
C     +-----------------------------------------------------------------+
C     | 3*ND + 1   ~  4*ND          | 滑移面法向 l 分量 (Global)        |
C     | 4*ND + 1   ~  5*ND          | 滑移面法向 m 分量 (Global)        |
C     | 5*ND + 1   ~  6*ND          | 滑移面法向 n 分量 (Global)        |
C     | 6*ND + 1   ~  7*ND          | 滑移方向 i 分量 (Global)          |
C     | 7*ND + 1   ~  8*ND          | 滑移方向 j 分量 (Global)          |
C     | 8*ND + 1   ~  9*ND          | 滑移方向 k 分量 (Global)          |
C     +-----------------------------------------------------------------+
C     | 滑移系演化变量 (9 ~ 12*ND)                                      |
C     +-----------------------------------------------------------------+
C     | 9*ND + 1   ~ 10*ND          | 各滑移系 |DGAMMA| 累积和          |
C     | 10*ND + 1  ~ 11*ND          | 各滑移系角度 ANG                  |
C     | 11*ND + 1  ~ 12*ND          | 各滑移系 SMD, SMDMAX              |
C     +-----------------------------------------------------------------+
C     | 变形梯度与张量 (12 ~ 14*ND) - 固定占用 2*ND 空间                |
C     +-----------------------------------------------------------------+
C     | 12*ND + 1  ~ 12*ND + 36     | 刚度矩阵张量 D(I,J) (6x6)         |
C     | 12*ND + 37 ~ 12*ND + 51     | (保留空间 / Padding)              |
C     | 14*ND - 9 ~ 14*ND           | 旋转张量 RROT(3,3) (结束于 420)   |
C     +-----------------------------------------------------------------+
C     | 位错密度变量 (14 ~ 18*ND)                                       |
C     +-----------------------------------------------------------------+
C     | 14*ND + 1  ~ 15*ND          | 各滑移系可动位错密度 (Mobile)     |
C     | 15*ND + 1  ~ 16*ND          | 备用                              |
C     | 16*ND + 1  ~ 17*ND          | 备用                              |
C     | 17*ND + 1  ~ 18*ND          | 备用                              |
C     +-----------------------------------------------------------------+
C     | 辐照缺陷密度 (18 ~ 21*ND)                                       |
C     +-----------------------------------------------------------------+
C     | 18*ND + 1  ~ 19*ND          | 各滑移系空隙密度 (Voids)          |
C     | 19*ND + 1  ~ 20*ND          | 各滑移系位错环密度 (Loops)        |
C     | 20*ND + 1  ~ 21*ND          | 各滑移系沉淀相密度 (Precipitates) |
C     +-----------------------------------------------------------------+
C     | 微观结构参数 (21 ~ 22*ND)                                       |
C     +-----------------------------------------------------------------+
C     | 21*ND + 1  ~ 22*ND          | 各滑移系平均自由程 (MFP)          |
C     +-----------------------------------------------------------------+
C     | 密度演化导数 (22 ~ 29*ND)                                       |
C     +-----------------------------------------------------------------+
C     | 22*ND + 1  ~ 23*ND          | d(可动位错)/d(滑移)               |
C     | 23*ND + 1  ~ 24*ND          | 备用!!可用于孪晶                  |
C     | 24*ND + 1  ~ 25*ND          | d(空隙)/d(滑移)                   |
C     | 25*ND + 1  ~ 26*ND          | d(位错环)/d(滑移)                 |
C     | 26*ND + 1  ~ 27*ND          | d(沉淀相)/d(滑移)                 |
C     | 27*ND + 1  ~ 28*ND          | 备用                              |
C     | 28*ND + 1  ~ 29*ND          | 备用                              |
C     +-----------------------------------------------------------------+
C     | 全局统计量 (位于数组末尾区域)                                   |
C     +-----------------------------------------------------------------+
C     | 29*ND + 1                   | 总可动位错密度 TOTMD              |
C     | 29*ND + 2                   | 备用                              |
C     | 29*ND + 3                   | 备用                              |
C     | 29*ND + 4                   | 备用                              |
C     | 29*ND + 5                   | 总空位密度 TOTVOID                |
C     | 29*ND + 6                   | 总位错环密度 TOTLOOP              |
C     | 29*ND + 7                   | 总沉淀相密度 TOTPRE               |
C     | 29*ND + 8                   | 备用                              |
C     | 29*ND + 9                   | 总自由程 TOTEMFP                  |
C     +-----------------------------------------------------------------+
C     | 滑移系分组信息 (位于 30*ND 之前)                                |
C     +-----------------------------------------------------------------+
C     | NSTATV - 30                 | 整体的塑性等效变形                |
C     | NSTATV - 29                 | 整体的孔隙率                      |
C     | NSTATV - 28                 | 各个滑移系切应变的加和            |
C     | NSTATV - 27                 | 各个孪晶系体积分数的加和          |
C     | NSTATV - 9                  | 第1组滑移系数量                   |
C     | NSTATV - 8                  | 第2组滑移系数量                   |
C     | NSTATV - 7                  | 第3组滑移系数量                   |
C     | NSTATV - 6                  | 第4组滑移系数量                   |
C     | NSTATV - 5                  | 第5组滑移系数量                   |
C     |                             |                                   |
C     | NSTATV - 3                  | 滑移组总数 NSET                   |
C     | NSTATV - 2                  | 孪晶系总数 NTWTL                  |
C     | NSTATV - 1                  | 滑移系总数 NSLPTL                 |
C     | NSTATV                      | 滑移系总数 + 孪晶系总数 ND        |
C     +-----------------------------------------------------------------+
!     
!     参考文献:
!       [1] 作者, 论文名, 期刊, 年份
!     
!     修改历史:
!       2025-11-27  v0.1  初始版本
!       2025-12-21  v0.2  引入孪晶
!       2025-12-25  v0.5  引入模块处理全局变量
!     
!     注意事项:
!       1. 需要定义NSTATV状态变量数量
!       2. 材料参数通过INP文件的*USER MATERIAL卡片输入
!***********************************************************************
!=======================================================================
!                   晶体塑性有限元常用常量参数模块
      MODULE PARABANK
  !---------------------------------------------------------------------
  ! 1. 晶体学本构(大部分在HCP_slip_twin.for中赋值)
  !---------------------------------------------------------------------
      INTEGER NSLPTL, NTWTL             ! 滑移系总数,孪晶系总数
      INTEGER NSLIP(1:5), NTWIN(1:2)    ! 每个滑移系类型有多少个滑移系
      INTEGER NSET, TWSET               ! 滑移组数量,孪晶组数量
      REAL*8 SIGMA0(1:42)               ! 各个滑移系的派纳力
      REAL*8 MPNS                       ! 滑移系的派纳力拟合参数
      REAL*8 BURGER(1:42)               ! 各个滑移系的Burgers矢量大小
      REAL*8 GMOD                       ! 剪切模量
      REAL*8 DGRAIN                     ! 晶粒尺寸
      REAL*8 HARDM(1:42,1:42)           ! 位错交互硬化矩阵
      REAL*8 HRHO, PMUL, PDYN           ! 位错硬化系数与密度演化参数
      REAL*8 DISLOCATION_DENSITY(1:5) ! 各滑移系的初始位错密度
      REAL*8 TAU0_TWIN(1:2)             ! 各个孪晶系的初始切应力
      REAL*8 H_SL_TW                    ! 孪晶对滑移交互硬化系数
      REAL*8 TWIN_D                     ! 孪晶对滑移交互硬化系数
      REAL*8 H_TW_TW                    ! 孪晶对孪晶交互硬化系数
      REAL*8 TWIN_B                     ! 孪晶对孪晶交互硬化系数
  !---------------------------------------------------------------------
  ! 2. 网格与晶粒拓扑
  !---------------------------------------------------------------------        
      INTEGER ND                        ! HCP滑移系为30, 30滑移系+12孪晶系为42
      PARAMETER(ND=42)
      INTEGER ELENUM                    ! 几何模型中的单元数
      PARAMETER(ELENUM=50000)
      INTEGER GRNUM                     ! 几何模型中的晶粒数
      PARAMETER(GRNUM=13)
  !---------------------------------------------------------------------
  ! 3. 求解器与流动法则控制
  !---------------------------------------------------------------------
      REAL*8 THETA                      ! 隐式求解参数theta
      PARAMETER(THETA=0.5D0)
      REAL*8 GAMDOT0                    ! 指数流动法则中的本构参数GAMDOT0
      PARAMETER(GAMDOT0=0.001D0)
      REAL*8 qEXP                       ! 指数流动法则中的本构参数qEXP
      PARAMETER(qEXP=20.0D0)
  !---------------------------------------------------------------------
  ! 4. 本构开关
  !---------------------------------------------------------------------
      INTEGER NLEGOM                    ! 小变形还是大变形
      PARAMETER (NLEGOM=1)            
      INTEGER IRADON                    ! 是否考虑辐照影响，1-考虑，0-不考虑
      PARAMETER (IRADON=0)
      INTEGER POROSITY                  ! 是否考虑孔隙率影响，1-考虑，0-不考虑
      PARAMETER (POROSITY=0)
      INTEGER LAYERSTRUCTURE            ! 是否考虑层状组织影响，1-层状组织，0-等轴组织
  !---------------------------------------------------------------------
  ! 4. 调试参数
  !---------------------------------------------------------------------      
      INTEGER PROBLEM                   ! 调试参数，1-出现问题，0-没有问题
      END MODULE PARABANK
!=======================================================================

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      USE PARABANK

      CHARACTER(len=80)::CMNAME
      character (len=200) site

      integer lenoutdir
!     ABAQUS变量声明
      REAL(kind=8)::STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

      INTEGER NOEL, KSTEP, KINC, IGRAIN, I, J, K, ID, IDBACK,
     2       J1, J2, I1, I2, L, ITRATN, ITRMAX,
     3       NITRTN, NDI, IDNOR, IDDIR, IE, IE_SW, NSHR,
     4       IDG, IJ2, NTENS, NSTATV, NPROPS, NPT, LAYER, KSPT
      INTEGER INDX(ND), ITRM(3)
      
      !     初始全局局部变化变量声明
      REAL(kind=8)::XRD(3),XND(3),XTD(3),PI,A(3),
     1 CROT(3,3),PHI1,PHI,PHI2,A11,A12,A13,A21,A22,A23,A31,
     2 A32, A33, PROPC(16)

!     CPFE变量
      REAL(kind=8)::SLPDIR(3,ND), SLPNOR(3,ND), SLPDEF(6,ND),
     3          SLPSPN(3,ND), DSPDIR(3,ND), DSPNOR(3,ND),
     4          DLOCAL(6,6), D(6,6), ROTD(6,6), ROTATE(3,3),
     5          FSLIP(ND), DFDXSP(ND), DDEMSD(6,ND),
     6          H(ND,ND), DDGDDE(ND,6),
     7          DSTRES(6), DELATS(6), DSPIN(3), DVGRAD(3,3),
     8          DGAMMA(ND), DTAUSP(ND), DGSLIP(ND),
     9          WORKST(ND,ND), TERM(3,3), TRM0(3,3), RROT(3,3)
      REAL*8 CHECK, GSHEAR, E11, E12, GAMERR, SLIP_SUM, TWIN_SUM
     3       DDCMP, DEV, TERM1, DTIME, TAUSLP, GSLIP,
     4       X, TERM2, TERM3, TERM4, TEMPWK, PNEWDT,
     5       ROMNOR, ROMDIR, DIRC, ANG, NORC, SMD, RESIDU, XMIS,
     6       SSE, SPD, SCD, RPL, DRPLDT, TEMP, DTEMP, CELENT

!     力学变量声明      
      REAL*8 HYSTR, EQVSTR, DEVSTR(6), EQVSTN, DEQVPLDT, EQVPL
!     孔隙率变量声明
      REAL(kind=8)::EQVSTR0, POROS, PORO_A, W

!     位错与缺陷变量声明（辐照硬化）
      REAL(kind=8)::NMD(ND), NVOID(ND), NLOOP(ND), NPRE(ND), EMFP(ND)
      REAL*8 TOTMD, TOTVOID, TOTLOOP, TOTPRE, TOTEMFP
      REAL*8 SVOID, SLOOP, SPRE
      REAL*8 HVOID, HLOOP, HPRE
      REAL(kind=8) :: E_DSPIN(3,3)
      REAL(kind=8) :: RRROT(3,3)
      REAL(kind=8) :: DMDDG(ND)
      REAL(kind=8) :: DLOOPDG(ND)

      DATA PI /3.14159265D0/
      DATA XRD/1.0D0, 0.0D0, 0.0D0/
      DATA XND/0.0D0, 0.0D0, 1.0D0/
      DATA XTD/0.0D0, 1.0D0, 0.0D0/
      REAL*8 RDX, RDY, RDZ, RD, TDX, TDY, TDZ, TD
      COMMON/CRYSPAR/ RDX, RDY, RDZ, RD, TDX, TDY, TDZ, TDEAL

!     晶粒方向声明与定义
      INTEGER::GELES(ELENUM)
      REAL(kind=8)::FGRAIN(3,GRNUM)

      INCLUDE 'GEULER.for'
      INCLUDE 'ELEGRN.for'
      INCLUDE 'FIT.for'
!------------------ 定义结束，可以放程序了 ------------------！
!     用于断点调试的代码，平时注释掉
c      write(*,*) "Please input an interger:"
c      read(*,*)

      PROBLEM = 0

!     输出文件变量声明与定义
      call getoutdir(site,lenoutdir)
      site=trim(site)//'/' 
      OPEN(UNIT=114514, FILE =trim(site)//'0test.txt')
      OPEN(UNIT=114515, FILE =trim(site)//'1DGAMMA.txt')
      OPEN(UNIT=114516, FILE =trim(site)//'2TAU_GSLIP.txt')
      OPEN(UNIT=114517, FILE =trim(site)//'3GSLIP0.txt')

!     初始化位错交互矩阵(这里设置了1-ND，是因为此时只有ND被定义了，实际有效大小是1-NSLPTL)
      DO I=1,ND
        DO J=1,ND
          IF (I.EQ.J) THEN
            HARDM(I,J)=1.0D0 !!这里可以按照PPT中的HCP硬化来修改
          ELSE
            HARDM(I,J)=1.4D0 !!交叉硬化系数，可以调整
          END IF
        END DO
      END DO

!     根据PROP25到30对坐标进行转换
      IF (KSTEP.EQ.1.AND.KINC.LE.1) THEN
        RD=DSQRT(PROPS(25)**2.0D0+PROPS(26)**2.0D0+PROPS(27)**2.0D0)
        IF (RD.EQ.0.0D0) THEN
          WRITE (6,*) '***ERROR*** - rolling dir = 0'
          STOP
        END IF
        RDX=PROPS(25)/RD
        RDY=PROPS(26)/RD
        RDZ=PROPS(27)/RD
        TD=DSQRT(PROPS(28)**2.0D0+PROPS(29)**2.0D0+PROPS(30)**2.0D0)
        IF (TD.EQ.0.0D0) THEN
          WRITE (6,*) '***ERROR*** - transverse dir = 0'
          STOP
        END IF
        TDX=PROPS(28)/TD
        TDY=PROPS(29)/TD
        TDZ=PROPS(30)/TD
        IF (DABS(RDX*TDX+RDY*TDY+RDZ*TDZ).GE.1.0D-3) THEN
          WRITE (6,*) '***ERROR*** - rolling dir is not normal'
          STOP
        END IF
      END IF

!     计算余弦矩阵
      A(1)=FGRAIN(1,IGRAIN)
      A(2)=FGRAIN(2,IGRAIN)
      A(3)=FGRAIN(3,IGRAIN)

!     计算欧拉角
      PHI1=PI/180.0D0*A(1)
      PHI=PI/180.0D0*A(2)
      PHI2=PI/180.0D0*A(3)

!     从轧制坐标到单晶基坐标的变化矩阵A[ij]
      A11=DCOS(PHI1)*DCOS(PHI2)-DSIN(PHI1)*DSIN(PHI2)*DCOS(PHI)
      A12=DSIN(PHI1)*DCOS(PHI2)+DCOS(PHI1)*DSIN(PHI2)*DCOS(PHI)
      A13=DSIN(PHI2)*DSIN(PHI)
      A21=-DCOS(PHI1)*DSIN(PHI2)-DSIN(PHI1)*DCOS(PHI2)*DCOS(PHI)
      A22=-DSIN(PHI1)*DSIN(PHI2)+DCOS(PHI1)*DCOS(PHI2)*DCOS(PHI)
      A23=DCOS(PHI2)*DSIN(PHI)
      A31=DSIN(PHI1)*DSIN(PHI)
      A32=-DCOS(PHI1)*DSIN(PHI)
      A33=DCOS(PHI)

!     全局坐标系下两个坐标的方向余弦      
      PROPC(1)=1.0D0
      PROPC(2)=0.0D0
      PROPC(3)=0.0D0
      PROPC(4)=A11*RDX+A12*TDX+A13*(RDY*TDZ-RDZ*TDY)
      PROPC(5)=A11*RDY+A12*TDY+A13*(RDZ*TDX-RDX*TDZ)
      PROPC(6)=A11*RDZ+A12*TDZ+A13*(RDX*TDY-RDY*TDX)
      PROPC(7)=0.0D0
      PROPC(8)=0.0D0      
      PROPC(9)=0.0D0
      PROPC(10)=1.0D0
      PROPC(11)=0.0D0
      PROPC(12)=A21*RDX+A22*TDX+A23*(RDY*TDZ-RDZ*TDY)
      PROPC(13)=A21*RDY+A22*TDY+A23*(RDZ*TDX-RDX*TDZ)
      PROPC(14)=A21*RDZ+A22*TDZ+A23*(RDX*TDY-RDY*TDX)
      PROPC(15)=0.0D0
      PROPC(16)=0.0D0

!     初始化局部坐标系下的弹性矩阵
      DO J=1,6
         DO I=1,6
            DLOCAL(I,J)=0.0D0
         END DO
      END DO

!     设置局部坐标系下的弹性矩阵(此处为HCP结构)
      DLOCAL(1,1)=PROPS(1)
      DLOCAL(2,2)=PROPS(1)
      DLOCAL(1,2)=PROPS(2)
      DLOCAL(2,1)=PROPS(2)
      DLOCAL(1,3)=PROPS(3)
      DLOCAL(3,1)=PROPS(3)
      DLOCAL(2,3)=PROPS(3)
      DLOCAL(3,2)=PROPS(3)
      DLOCAL(3,3)=PROPS(4)
      DLOCAL(5,5)=PROPS(5)
      DLOCAL(6,6)=PROPS(5)
      DLOCAL(4,4)=(DLOCAL(1,1)-DLOCAL(1,2))/2.

!     旋转矩阵
      DO I=1,3
       DO J=1,3
         CROT(I,J)=0.0D0
       END DO
      END DO
      CALL ROTATION (PROPC, ROTATE)

!     作用到全局坐标系下
      DO J=1,3
         J1=1+J/3
         J2=2+J/2
         DO I=1,3
            I1=1+I/3
            I2=2+I/2
            ROTD(I,J)=ROTATE(I,J)**2.0D0
            ROTD(I,J+3)=2.0D0*ROTATE(I,J1)*ROTATE(I,J2)
            ROTD(I+3,J)=ROTATE(I1,J)*ROTATE(I2,J)
            ROTD(I+3,J+3)=ROTATE(I1,J1)*ROTATE(I2,J2)+
     2                    ROTATE(I1,J2)*ROTATE(I2,J1)
         END DO
      END DO

!     最终得到全局坐标系下的切线刚度矩阵
      DO J=1,6
         DO I=1,6
            D(I,J)=0.0D0
         END DO
      END DO
      DO J=1,6
         DO I=1,J
            DO K=1,6
               DO L=1,6
                  D(I,J)=D(I,J)+DLOCAL(K,L)*ROTD(I,K)*ROTD(J,L)
               END DO
            END DO
            D(J,I)=D(I,J)
         END DO
      END DO

!     得到DSPIN旋转速度张量增量
      IF (NLGEOM.NE.0) THEN
         DO J=1,3
            DO I=1,3
               TERM(I,J)=DROT(J,I)
               TRM0(I,J)=DROT(J,I)
            END DO
            TERM(J,J)=TERM(J,J)+1.0D0
            TRM0(J,J)=TRM0(J,J)-1.0D0
         END DO
         CALL LUDCMP (TERM, 3, 3, ITRM, DDCMP)
         DO J=1,3
            CALL LUBKSB (TERM, 3, 3, ITRM, TRM0(1,J))
         END DO
         DSPIN(1)=TRM0(2,1)-TRM0(1,2)
         DSPIN(2)=TRM0(1,3)-TRM0(3,1)
         DSPIN(3)=TRM0(3,2)-TRM0(2,3)
      END IF  

!     计算因为弹性部分导致的体积变化
      DEV=0.0D0
      DO I=1,NDI
         DEV=DEV+DSTRAN(I)
      END DO
      NITRTN = -1

!======================= 开始迭代 =========================!
1000  CONTINUE
      NITRTN = NITRTN+1
!     判断是否是第一次状态，如果是的话进行初始化
! ---------------------- 初始模块 -------------------------!
      IF (STATEV(1).EQ.0.0D0) THEN
        IF (NOEL.eq.1) THEN
            write(*,*) 'Initial state, Time=', Time
        END IF
        DO I=1,ND
          NMD(I)=0.0D0 !定义可动位错密度
          NVOID(I)=0.0D0 !定义空位密度
          NLOOP(I)=0.0D0 !定义位错环密度
          NPRE(I)=0.0D0 !定义析出相密度
          EMFP(I)=0.0D0 !定义有效平均自由程
          DMDDG(I)=0.0D0
          DLOOPDG(I)=0.0D0
        END DO
        TOTMD=0.0D0 !所有滑移系的可动位错
        TOTVOID=0.0D0 !所有滑移系的空位总和
        TOTLOOP=0.0D0 !所有滑移系的位错环总和
        TOTPRE=0.0D0 !所有滑移系的析出相密度总和
        TOTEMFP=0.0D0 !所有滑移系的平均有效自由程总和

!     得到HCP中的所有初始滑移系信息
        CALL HCP_slip_twin(SLPDIR,SLPNOR,ROTATE) !!目前使用的子程序是HCP_slip_twin
        IF (ND.LT.NSLPTL) THEN
            WRITE (6,*)
     2 '***ERROR - parameter ND chosen by the present user is less than
     3             the total number of slip systems NSLPTL'
            STOP
        END IF
          DO J=1,ND !!计算每个滑移系+孪晶系的SLPDEF
            SLPDEF(1,J)=SLPDIR(1,J)*SLPNOR(1,J)
            SLPDEF(2,J)=SLPDIR(2,J)*SLPNOR(2,J)
            SLPDEF(3,J)=SLPDIR(3,J)*SLPNOR(3,J)
            SLPDEF(4,J)=SLPDIR(1,J)*SLPNOR(2,J)+SLPDIR(2,J)*SLPNOR(1,J)
            SLPDEF(5,J)=SLPDIR(1,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(1,J)
            SLPDEF(6,J)=SLPDIR(2,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(2,J)      
          END DO
!     初始化状态变量：滑移系方向和法向
        STATEV(NSTATV-3)=DBLE(NSET)   !倒数第三个变量储存滑移组数量
        STATEV(NSTATV-2)=DBLE(NTWTL)  !倒数第二个变量储存孪晶系数量
        STATEV(NSTATV-1)=DBLE(NSLPTL) !倒数第一个变量储存滑移系数量
        DO I=1, NSET
          STATEV(NSTATV-10+I)=dble(NSLIP(I)) !!储存每个滑移组的滑移系数，这里的NSET包括孪晶的伪滑移系
        END DO
        IDNOR=3*ND
        IDDIR=6*ND
        DO J=1, ND !将滑移系方向和法相储存在STATEV里
          DO I=1,3
            IDNOR=IDNOR+1
            STATEV(IDNOR)=SLPNOR(I,J)
            IDDIR=IDDIR+1
            STATEV(IDDIR)=SLPDIR(I,J)
          END DO
        END DO

!     对滑移系的缺陷数据和位错数据进行初始化
        CALL DIS_IRRA_INIT(NOEL,NMD,NVOID,NLOOP,NPRE,EMFP,
     2  SVOID,SLOOP,SPRE,TOTMD,TOTVOID,TOTLOOP,TOTPRE,TOTEMFP,SLPDIR,SLPNOR,
     3  HVOID, HLOOP, HPRE) !!需要修改
        DO I=1,NSLPTL
          STATEV(14*ND+I)=NMD(I)
          STATEV(18*ND+I)=NVOID(I)
          STATEV(19*ND+I)=NLOOP(I)
          STATEV(20*ND+I)=NPRE(I)
          STATEV(21*ND+I)=EMFP(I)
        END DO
        STATEV(29*ND+1)=TOTMD
        STATEV(29*ND+5)=TOTVOID
        STATEV(29*ND+6)=TOTLOOP
        STATEV(29*ND+7)=TOTPRE
        STATEV(29*ND+9)=TOTEMFP
        
        CALL DIS_GSLP_INIT(STATEV(1),NMD,
     2  NVOID,NLOOP,NPRE,SVOID,SLOOP,SPRE)!!需要修改

!     初始化孪晶系的切变阻力
        CALL TWIN_TAU_INIT(STATEV(NSLPTL+1)) !!需要修改
c        DO I=1,ND
c          write(114517,*) 'TAU0' , I, STATEV(I)
c        END DO

!     初始化每个滑移系+孪晶系的剪切应变
        DO I=1,ND
          STATEV(ND+I)=0.0D0 !各个滑移系产生的滑移量
          STATEV(9*ND+I)=0.0D0 !各个滑移系的累积滑移量的绝对值
        END DO
        SLIP_SUM=0.0D0
        TWIN_SUM=0.0D0
        STATEV(NSTATV-28)=SLIP_SUM !各个滑移系切应变的加和
        STATEV(NSTATV-27)=TWIN_SUM !各个孪晶系体积分数的加和

!     初始化每个滑移系+孪晶系的分切应力(存在STATEV(2*ND+I)中)
        DO I=1,ND
          TERM1=0.0D0
          DO J=1,NTENS
            IF (J.LE.NDI) THEN
              TERM1=TERM1+SLPDEF(J,I)*STRESS(J)
            ELSE
              TERM1=TERM1+SLPDEF(J-NDI+3,I)*STRESS(J)
            END IF
          END DO
          STATEV(2*ND+I)=TERM1
        END DO
!     初始化旋转矩阵RROT，用于计算后续
        DO I=1, 3
          DO J=1, 3
            RROT(I,J)=0.0D0
          IF(I.EQ.J)RROT(I,J)=1.0D0
          END DO
        END DO
! -------------------- 结束初始模块 -----------------------!
!     如果不是初始状态，从状态变量中读取信息      
      ELSE 
        IF(NOEL.eq.1) THEN
            write(*,*) 'Current state, Time=', Time
        END IF
! -------------------- 读取变量模块 -----------------------!
!     提取滑移系信息
        NSET=INT(STATEV(NSTATV-3))   !提取滑移组总数
        NTWTL=INT(STATEV(NSTATV-2))  !提取孪晶系总数
        NSLPTL=INT(STATEV(NSTATV-1)) !提取滑移系总数
!     提取每个滑移组的滑移系数量
        DO I=1, NSET
          NSLIP(I)=INT(STATEV(NSTATV-10+I))
        END DO
!     从状态变量中提取滑移系+孪晶系方向和法向
        IDNOR=3*ND
        IDDIR=6*ND
        DO J=1, ND
          DO I=1,3
            IDNOR=IDNOR+1
            SLPNOR(I,J)=STATEV(IDNOR)
            IDDIR=IDDIR+1
            SLPDIR(I,J)=STATEV(IDDIR) 
          END DO
        END DO
!     从状态变量中提取位错密度和缺陷密度
        DO I=1, NSLPTL
          NMD(I)=STATEV(14*ND+I)
          NVOID(I)=STATEV(18*ND+I)
          NLOOP(I)=STATEV(19*ND+I)
          NPRE(I)=STATEV(20*ND+I)
          EMFP(I)=STATEV(21*ND+I)
          DMDDG(I)=STATEV(22*ND+I)
          DLOOPDG(I)=STATEV(25*ND+I)
        END DO
        TOTMD=STATEV(29*ND+1)
        TOTVOID=STATEV(29*ND+5)
        TOTLOOP=STATEV(29*ND+6)
        TOTPRE=STATEV(29*ND+7)
        TOTEMFP=STATEV(29*ND+9)
!     计算此时的施密特因子(滑移系+孪晶系)
        DO J=1,ND
           SLPDEF(1,J)=SLPDIR(1,J)*SLPNOR(1,J)
           SLPDEF(2,J)=SLPDIR(2,J)*SLPNOR(2,J)
           SLPDEF(3,J)=SLPDIR(3,J)*SLPNOR(3,J)
           SLPDEF(4,J)=SLPDIR(1,J)*SLPNOR(2,J)+SLPDIR(2,J)*SLPNOR(1,J)
           SLPDEF(5,J)=SLPDIR(1,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(1,J)
           SLPDEF(6,J)=SLPDIR(2,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(2,J)
        END DO        
!    读取此时的旋转张量RROT
        K=0
        DO I=1,3
          DO J=1,3
            K=K+1
            RROT(I,J)=STATEV(14*ND-10+K)
          END DO
        END DO
!    从状态变量中读取全局刚度矩阵D
        IE_SW=12*ND
        IE=0
        DO I=1,6
          DO J=1,6
            IE=1+IE
            D(I,J)=STATEV(IE_SW+IE)
          ENDDO
        ENDDO
! ------------------- 结束读取变量模块 ---------------------!
      END IF

!======================= 计算主体 =========================!
!     计算滑移系+孪晶系旋转张量
      IF (NLGEOM.NE.0) THEN
         DO J=1,ND
            SLPSPN(1,J)=0.5D0*(SLPDIR(1,J)*SLPNOR(2,J)-SLPDIR(2,J)*SLPNOR(1,J))
            SLPSPN(2,J)=0.5D0*(SLPDIR(3,J)*SLPNOR(1,J)-SLPDIR(1,J)*SLPNOR(3,J))
            SLPSPN(3,J)=0.5D0*(SLPDIR(2,J)*SLPNOR(3,J)-SLPDIR(3,J)*SLPNOR(2,J))
         END DO
      END IF
!     双点乘弹性模量和施密特因子DDEMSD
      DO J=1,ND
        DO I=1,6
          DDEMSD(I,J)=0.0D0
          DO K=1,6
            DDEMSD(I,J)=DDEMSD(I,J)+D(K,I)*SLPDEF(K,J)
          END DO
        END DO
      END DO
      IF (NLGEOM.NE.0) THEN
        DO J=1,ND
          DDEMSD(4,J)=DDEMSD(4,J)-SLPSPN(1,J)*STRESS(1)
          DDEMSD(5,J)=DDEMSD(5,J)+SLPSPN(2,J)*STRESS(1)
          IF (NDI.GT.1) THEN
            DDEMSD(4,J)=DDEMSD(4,J)+SLPSPN(1,J)*STRESS(2)
            DDEMSD(6,J)=DDEMSD(6,J)-SLPSPN(3,J)*STRESS(2)
          END IF
          IF (NDI.GT.2) THEN
            DDEMSD(5,J)=DDEMSD(5,J)-SLPSPN(2,J)*STRESS(3)
            DDEMSD(6,J)=DDEMSD(6,J)+SLPSPN(3,J)*STRESS(3)
          END IF
          IF (NSHR.GE.1) THEN
            DDEMSD(1,J)=DDEMSD(1,J)+SLPSPN(1,J)*STRESS(NDI+1)
            DDEMSD(2,J)=DDEMSD(2,J)-SLPSPN(1,J)*STRESS(NDI+1)
            DDEMSD(5,J)=DDEMSD(5,J)-SLPSPN(3,J)*STRESS(NDI+1)
            DDEMSD(6,J)=DDEMSD(6,J)+SLPSPN(2,J)*STRESS(NDI+1)
          END IF
          IF (NSHR.GE.2) THEN
            DDEMSD(1,J)=DDEMSD(1,J)-SLPSPN(2,J)*STRESS(NDI+2)
            DDEMSD(3,J)=DDEMSD(3,J)+SLPSPN(2,J)*STRESS(NDI+2)
            DDEMSD(4,J)=DDEMSD(4,J)+SLPSPN(3,J)*STRESS(NDI+2)
            DDEMSD(6,J)=DDEMSD(6,J)-SLPSPN(1,J)*STRESS(NDI+2)
          END IF
          IF (NSHR.EQ.3) THEN
            DDEMSD(2,J)=DDEMSD(2,J)+SLPSPN(3,J)*STRESS(NDI+3)
            DDEMSD(3,J)=DDEMSD(3,J)-SLPSPN(3,J)*STRESS(NDI+3)
            DDEMSD(4,J)=DDEMSD(4,J)-SLPSPN(2,J)*STRESS(NDI+3)
            DDEMSD(5,J)=DDEMSD(5,J)+SLPSPN(1,J)*STRESS(NDI+3)
          END IF
        END DO
      END IF
!     调用子程序计算滑移率FSLIP和剪切应力DFDXSP          
      CALL STRAINRATE (STATEV(2*ND+1),STATEV(1),FSLIP(1),DFDXSP(1)) !!需要修改
! ------------------- 读取计算力学性质 ---------------------!     
!     计算正应力
      HYSTR=0.0D0
      DO I=1,NDI
        HYSTR=HYSTR+STRESS(I)
      END DO
      HYSTR=HYSTR/3.0D0
!     计算偏应力
      DO I=1,6
        DEVSTR(I)=0.0D0
      END DO 
      DO I=1, NDI
        DEVSTR(I)=STRESS(I)-HYSTR
      END DO
      DO I=1, NSHR
        DEVSTR(3+I)=STRESS(NDI+I)
      END DO
      IF (NDI.LT.3)THEN
        DO I=NDI+1, 3
          DEVSTR(I)=0.0D0-HYSTR
        END DO
      END IF
!     计算等效应力
      EQVSTR=0.0D0
      DO I=1,6
        EQVSTR=EQVSTR+DEVSTR(I)*DEVSTR(I)
      END DO
      EQVSTR=DSQRT(1.5D0*EQVSTR)
!     计算等效塑性应变率
      DEQVPLDT=0.0D0
      DO I =1,ND
        DEQVPLDT=DEQVPLDT+STATEV(2*ND+I)*FSLIP(I)
        !FSLIP(I)在STRAINRATE子程序中计算得到,为每个滑移系+孪晶系的切变率
      END DO
      IF(EQVSTR.NE.0.0D0) DEQVPLDT=W/EQVSTR
!     从状态变量中读取总的等效塑性应变
      EQVPL = STATEV(NSTATV-30)

      IF (POROSITY.EQ.1) THEN
!     从状态变量中读取孔隙率
        POROS = STATEV(NSTATV-29)
!     计算基体的流动应力
        W = 0.0D0
!!        CALL CALC_WEAK_FUNCTION(W, POROS, PORO_A)
        EQVSTR0=EQVSTR/W
      END IF
! ------------------ 结束读取计算力学性质 -------------------!

! ---------------------- 进入硬化模块 ----------------------!
      CALL DIS_IRRA_HARD(H,NOEL,NMD,NVOID,NLOOP,NPRE,
     2  SVOID,SLOOP,SPRE,DMDDG,DLOOPDG,EMFP) !!需要修改
      IF (POROSITY.EQ.1) THEN
!!        CALL POROHARD()!!需要修改
      END IF
c      IF (NOEL.EQ.1) write(114514, *) 'H', H 
      CALL TWIN_HARD(H,STATEV(NSTATV-28),STATEV(NSTATV-27)) !!需要修改
! --------------------- 结束硬化模块 -----------------------!

! ---------------------- 计算WORKST -----------------------!
      TERM1=THETA*DTIME
!     读取滑移系+孪晶系的剪切应力和强度
      DO I=1,ND
        TAUSLP=STATEV(2*ND+I)
        GSLIP=STATEV(I)
        TERM2=TERM1*DFDXSP(I)/GSLIP
        TERM3=TERM1*X*DFDXSP(I)/GSLIP
        DO J=1,ND
          TERM4=0.0D0
          DO K=1,6
            TERM4=TERM4+DDEMSD(K,I)*SLPDEF(K,J)
          END DO 
          TEMPWK=TERM2*TERM4+H(I,J)*TERM3*DSIGN(1.D0,FSLIP(J))
!     PNEWDT为Abaqus的时间步长控制量，如果TEMPWK出现NaN则减小时间步长
          IF (ISNAN(TEMPWK)) THEN
            PNEWDT=0.1D0 
          ELSE
            WORKST(I,J)=TEMPWK
          END IF
          IF (NITRTN.GT.0) THEN
            IF (ISNAN(TEMPWK)) THEN
              PNEWDT=0.1D0
            ELSE
              WORKST(I,J)=TEMPWK
            END IF
          END IF
        END DO
        WORKST(I,I)=WORKST(I,I)+1.0D0
!     计算滑移率增量DGAMMA 
        X=TAUSLP/GSLIP
        DGAMMA(I)=0.0D0
        DO J=1,NDI
            DGAMMA(I)=DGAMMA(I)+DDEMSD(J,I)*DSTRAN(J)
          END DO
          IF (NSHR.GT.0) THEN
            DO J=1,NSHR
              DGAMMA(I)=DGAMMA(I)+DDEMSD(J+3,I)*DSTRAN(J+NDI)
            END DO
          END IF
        DGAMMA(I)=DGAMMA(I)*TERM2+FSLIP(I)*DTIME
      END DO
        
      CALL LUDCMP (WORKST, ND, ND, INDX, DDCMP)
      CALL LUBKSB (WORKST, ND, ND, INDX, DGAMMA)

      DO I=1,ND
        IF (ISNAN(DGAMMA(I))) THEN
          DGAMMA=0.0D0
          PNEWDT=0.1D0
        END IF
      END DO
! -------------------- 结束计算WORKST --------------------!
!====================== 计算主体结束 ======================!

!====================== 更新状态变量 ======================!      
!     更新滑移应变和累积滑移应变
      DO I=1,ND
        STATEV(I+ND)=STATEV(I+ND)+DGAMMA(I)
        STATEV(9*ND+I)=STATEV(9*ND+I)+DABS(DGAMMA(I))
      END DO
      DO I=1,NSLPTL
        SLIP_SUM = SLIP_SUM+STATEV(9*ND+I)
      END DO
      !!这一部分指的是拉伸孪晶部分的参考滑移系数是0.175，即切应变/孪晶体积分数=0.175
      DO I=NSLPTL+1,NSLPTL+NTWIN(1) 
        TWIN_SUM = TWIN_SUM+STATEV(9*ND+I)/0.175D0
      END DO
      !!这一部分指的是压缩孪晶部分的参考滑移系数是0.175，即切应变/孪晶体积分数=0.218
      DO I=NSLPTL+NTWIN(1)+1,ND
        TWIN_SUM = TWIN_SUM+STATEV(9*ND+I)/0.218D0
      END DO
      STATEV(NSTATV-28)=SLIP_SUM !更新所有滑移系切应变的加和
      STATEV(NSTATV-27)=TWIN_SUM !更新所有孪晶系切应变的加和

!     计算滑移系+孪晶系强度增量
      DO I=1,ND
        DGSLIP(I)=0.0D0
        DO J=1,ND
          DGSLIP(I)=DGSLIP(I)+H(I,J)*DGAMMA(J)
        END DO
      END DO

!     更新滑移系+孪晶系强度
      DO I=1,ND
        STATEV(I)=STATEV(I)+DGSLIP(I)
      END DO

!     计算晶格拉伸应变
      DO J=1,6
         DELATS(J)=0.0D0
      END DO
      DO J=1,3
         IF (J.LE.NDI) DELATS(J)=DSTRAN(J)
         DO I=1,ND
            DELATS(J)=DELATS(J)-SLPDEF(J,I)*DGAMMA(I)
         END DO
      END DO
      DO J=1,3
         IF (J.LE.NSHR) DELATS(J+3)=DSTRAN(J+NDI)
         DO I=1,ND
            DELATS(J+3)=DELATS(J+3)-SLPDEF(J+3,I)*DGAMMA(I)
         END DO
      END DO

!     计算晶格拉伸所对应的变形梯度增量DVGRAD(仅仅在大应变时需要计算)
      IF (NLGEOM.NE.0) THEN
        DO J=1,3
          DO I=1,3
            IF (I.EQ.J) THEN
              DVGRAD(I,J)=DELATS(I)
            ELSE
              DVGRAD(I,J)=DELATS(I+J+1)
            END IF
          END DO
        END DO
          DO J=1,3
            DO I=1,J
              IF (J.GT.I) THEN
                IJ2=I+J-2
                IF (MOD(IJ2,2).EQ.1) THEN
                  TERM1=1.0D0
                ELSE
                  TERM1=-1.0D0
                END IF
                  DVGRAD(I,J)=DVGRAD(I,J)+TERM1*DSPIN(IJ2)
                  DVGRAD(J,I)=DVGRAD(J,I)-TERM1*DSPIN(IJ2)
                DO K=1,ND
                  DVGRAD(I,J)=DVGRAD(I,J)-TERM1*DGAMMA(K)*SLPSPN(IJ2,K)
                  DVGRAD(J,I)=DVGRAD(J,I)+TERM1*DGAMMA(K)*SLPSPN(IJ2,K)
                END DO
              END IF
          END DO
        END DO
      END IF

!     计算滑移系的分切应力
      DO I=1,ND
        DTAUSP(I)=0.0D0
        DO J=1,6
          DTAUSP(I)=DTAUSP(I)+DDEMSD(J,I)*DELATS(J)
        END DO
      END DO

!888888888888888888888888888888888!
!!    DEBUG: 检查DTAUSP
c      WRITE(114516,*) 'DTAUSP'
c      WRITE(114516,*)  DTAUSP
c      WRITE(114515,*) 'DGAMMA'
c      WRITE(114515,*)  DGAMMA
!888888888888888888888888888888888!

!     更新滑移系的分切应力
      DO I=1,ND
        STATEV(2*ND+I)=STATEV(2*ND+I)+DTAUSP(I)
      END DO

!     计算材料应力增量
      IF (NLGEOM.EQ.0) THEN
        DO I=1,NTENS
          DSTRES(I)=0.0D0
        END DO
      ELSE
        DO I=1,NTENS
          DSTRES(I)=-STRESS(I)*DEV
        END DO
      END IF
      DO I=1,NDI
        DO J=1,NDI
          DSTRES(I)=DSTRES(I)+D(I,J)*DSTRAN(J)
        END DO
        IF (NSHR.GT.0) THEN
          DO J=1,NSHR
            DSTRES(I)=DSTRES(I)+D(I,J+3)*DSTRAN(J+NDI)
          END DO
        END IF
        DO J=1,ND   
          DSTRES(I)=DSTRES(I)-DDEMSD(I,J)*DGAMMA(J)   
        END DO
      END DO
      IF (NSHR.GT.0) THEN
        DO I=1,NSHR
          DO J=1,NDI
            DSTRES(I+NDI)=DSTRES(I+NDI)+D(I+3,J)*DSTRAN(J)
          END DO
          DO J=1,NSHR
            DSTRES(I+NDI)=DSTRES(I+NDI)+D(I+3,J+3)*DSTRAN(J+NDI)
          END DO
          DO J=1,ND
            DSTRES(I+NDI)=DSTRES(I+NDI)-DDEMSD(I+3,J)*DGAMMA(J)
          END DO
        END DO
      END IF
!     将更新应力存入状态变量中
      DO I=1,NTENS
        STRESS(I)=STRESS(I)+DSTRES(I)
      END DO
!     更新滑移系方向和法向
C      DO I=1, NSLPTL
C        ANOR(I)=0.0D0
C        ADIR(I)=0.0D0
C      END DO
!     计算滑移系法向增量DSPNOR和滑移系方向增量DSPDIR
      IF (NLGEOM.NE.0) THEN
        DO J=1,ND
          DO I=1,3
            DSPNOR(I,J)=0.0D0
            DSPDIR(I,J)=0.0D0
            DO K=1,3
              DSPNOR(I,J)=DSPNOR(I,J)-SLPNOR(K,J)*DVGRAD(K,I)
              DSPDIR(I,J)=DSPDIR(I,J)+SLPDIR(K,J)*DVGRAD(I,K)
            END DO
          END DO
        END DO
!     更新滑移系法向和方向
        IDNOR=3*ND
        IDDIR=6*ND
        DO J=1,ND
          DO I=1,3
            IDNOR=IDNOR+1
            STATEV(IDNOR)=STATEV(IDNOR)+DSPNOR(I,J)
            IDDIR=IDDIR+1
            STATEV(IDDIR)=STATEV(IDDIR)+DSPDIR(I,J)
          END DO
!     计算滑移系法向和方向的模长ROMNOR和ROMDIR，用于后续归一化
          ROMNOR=DSQRT(STATEV(IDNOR-2)*STATEV(IDNOR-2)+STATEV(IDNOR-1)*
     2                  STATEV(IDNOR-1)+STATEV(IDNOR)*STATEV(IDNOR))
          ROMDIR=DSQRT(STATEV(IDDIR-2)*STATEV(IDDIR-2)+STATEV(IDDIR-1)*
     2                  STATEV(IDDIR-1)+STATEV(IDDIR)*STATEV(IDDIR))
!     计算滑移系方向和Z轴的夹角余弦DIRC和夹角ANG，并储存
         DIRC=STATEV(IDDIR)/ROMDIR
         ANG=DACOS(DIRC)*180.0D0/PI
         STATEV(10*ND+J)=ANG
!     计算滑移系法向和Z轴的夹角余弦NORC和与对Z轴的施密特因子SMD，并储存      
         NORC=STATEV(IDNOR)/ROMNOR
         SMD=NORC*DIRC
         STATEV(11*ND+J)=SMD
        END DO
      END IF
!     位错与辐照缺陷的演化
      CALL DIS_IRRA_EVOLUTION(DMDDG(1),DLOOPDG(1),NOEL,EMFP(1),
     2  TOTEMFP,TOTMD,NMD,DGAMMA,SLPDIR,SLPNOR)
!     更新位错密度和缺陷密度到状态变量中
      DO I=1, NSLPTL
        STATEV(14*ND+I)=NMD(I)
        STATEV(18*ND+I)=NVOID(I)
        STATEV(19*ND+I)=NLOOP(I)
        STATEV(20*ND+I)=NPRE(I)
        STATEV(21*ND+I)=EMFP(I)
        STATEV(22*ND+I)=DMDDG(I)
        STATEV(25*ND+I)=DLOOPDG(I)
        
      END DO
      STATEV(29*ND+1)=TOTMD
      STATEV(29*ND+5)=TOTVOID
      STATEV(29*ND+6)=TOTLOOP
      STATEV(29*ND+7)=TOTPRE
      STATEV(29*ND+0)=TOTEMFP
!     修改和存储与孔隙率相关的变量
      IF (POROSITY.EQ.1) THEN
c        CALL POROEVOL()!!需要修改
        STATEV(NSTATV-29)=POROS
        STATEV(NSTATV-30)=EQVPL
        IF (NOEL.eq.1) THEN
          write(6,*) 'Current porosity=', POROS
        END IF
      END IF
!     计算DDSDDE的前置矩阵DDGDDE(滑移率增量对于变形增量的偏导)
      TERM1=THETA*DTIME
      DO I=1,NTENS
         DO J=1,ND
            TAUSLP=STATEV(2*ND+J)
            GSLIP=STATEV(J)    
            X=TAUSLP/GSLIP
            TERM2=TERM1*DFDXSP(J)/GSLIP
            IF (I.LE.NDI) THEN
               DDGDDE(J,I)=TERM2*DDEMSD(I,J)
            ELSE
               DDGDDE(J,I)=TERM2*DDEMSD(I-NDI+3,J)
            END IF
         END DO
         CALL LUBKSB (WORKST, ND, ND, INDX, DDGDDE(1,I))
      END DO
!     计算切线刚度矩阵DDSDDE
!     切线刚度矩阵的弹性部分
      DO J=1,NTENS
        DO I=1,NTENS
          DDSDDE(I,J)=0.0D0
        END DO
      END DO
      DO J=1,NDI
        DO I=1,NDI
          DDSDDE(I,J)=D(I,J)
          IF (NLGEOM.NE.0) DDSDDE(I,J)=DDSDDE(I,J)-STRESS(I)
        END DO
      END DO
      IF (NSHR.GT.0) THEN
        DO J=1,NSHR
          DO I=1,NSHR
            DDSDDE(I+NDI,J+NDI)=D(I+3,J+3)
          END DO
          DO I=1,NDI
          DDSDDE(I,J+NDI)=D(I,J+3)
          DDSDDE(J+NDI,I)=D(J+3,I)
          IF (NLGEOM.NE.0) DDSDDE(J+NDI,I)=DDSDDE(J+NDI,I)-STRESS(J+NDI)
          END DO
        END DO
      END IF
!     切线刚度矩阵的塑性部分
      DO J=1,NDI
        DO I=1,NDI
          DO K=1,ND             
            DDSDDE(I,J)=DDSDDE(I,J)-DDEMSD(I,K)*DDGDDE(K,J)
          END DO
        END DO
      END DO
      IF (NSHR.GT.0) THEN
        DO J=1,NSHR
          DO I=1,NSHR
            DO K=1,ND
              DDSDDE(I+NDI,J+NDI)=DDSDDE(I+NDI,J+NDI)-DDEMSD(I+3,K)*DDGDDE(K,J+NDI)
            END DO
          END DO
          DO I=1,NDI
            DO K=1,ND
              DDSDDE(I,J+NDI)=DDSDDE(I,J+NDI)-DDEMSD(I,K)*DDGDDE(K,J+NDI)
              DDSDDE(J+NDI,I)=DDSDDE(J+NDI,I)-DDEMSD(J+3,K)*DDGDDE(K,I)
            END DO
          END DO
        END DO
      END IF
!888888888888888888888888888888888!
!!    DEBUG: 检查DDSDDE
c      WRITE(114517,*) 'DDSDDE'
c      WRITE(114517,*)  DDSDDE
c      WRITE(114517,*) 'D'
c      WRITE(114517,*)  D
!888888888888888888888888888888888!
!     计算晶格旋转张量(来自弹性应变的旋转)
!     初始化晶格旋转增量E_DSPIN
      DO I=1,3
        DO J=1,3
          E_DSPIN(I,J)=0.0D0
        ENDDO
      ENDDO
!     计算晶格旋转增量E_DSPIN=P_DSTRAN-SUM(DGAMMA*SLPSPN)
      DO J=1,3
        DO I=1,J
          IF (J.GT.I) THEN
            IJ2=I+J-2
            IF (MOD(IJ2,2).EQ.1) THEN
              TERM1=1.0D0
            ELSE
              TERM1=-1.0D0
            END IF
            E_DSPIN(I,J)=E_DSPIN(I,J)+TERM1*DSPIN(IJ2)
            E_DSPIN(J,I)=E_DSPIN(J,I)-TERM1*DSPIN(IJ2)
            DO K=1,ND
              E_DSPIN(I,J)=E_DSPIN(I,J)-TERM1*DGAMMA(K)*SLPSPN(IJ2,K)
              E_DSPIN(J,I)=E_DSPIN(J,I)+TERM1*DGAMMA(K)*SLPSPN(IJ2,K)
            END DO
          END IF
        END DO
      END DO
      DO I = 1, 3
        E_DSPIN(I,I) = E_DSPIN(I,I)+1.0D0
      END DO
!     更新旋转矩阵RROT=E_DSPIN*RROT
      CALL MULTIP(E_DSPIN,3,3,RROT,3, RRROT)        
      DO I=1,3
        DO J=1,3
          RROT(I,J)=RRROT(I,J)
        END DO
      END DO
      K=0
      DO I=1,3
        DO J=1,3
          K=K+1
          STATEV(14*ND-10+K)=RROT(I,J)
        END DO
      END DO
!     根据新的旋转矩阵计算当前情况下的取向矩阵CROT(晶体局部坐标系到全局坐标系)
      CALL MULTIP(RROT,3,3,ROTATE,3,CROT)    
!     根据新的取向矩阵CROT计算材料的全局弹性刚度矩阵
      CALL ELASTO(DLOCAL, CROT, D)
!     将更新后的刚度矩阵D存入状态变量中
      IE_SW=12*ND
      IE=0
      DO I=1,6
        DO J=1,6
          IE=1+IE
          STATEV(IE_SW+IE)=D(I,J)
        END DO
      END DO
!888888888888888888888888888888888!
!!    DEBUG: 检查D
c      WRITE(114517,*) 'D_AFTER'
c      WRITE(114517,*)  D
!888888888888888888888888888888888!
!==================== 更新状态变量结束 =====================!
!     主程序结束
      RETURN
      END
!***********************************************************************

!***************************  基础子程序  ******************************
!-------------------- 1. 矩阵乘法子程序MULTIP -------------------------!
      SUBROUTINE MULTIP(A,N,M,B,L,C)
      IMPLICIT NONE
C-----  MATRIX MULTIPLY
      INTEGER I,J,K,N,M,L
      real(kind=8)::A(N,M),B(M,L),C(N,L)
      DO 10 I=1,N
         DO 10 J=1,L
            C(I,J)=0.0D0
            DO 10 K=1,M
               C(I,J)=C(I,J)+A(I,K)*B(K,J)
10    CONTINUE
      RETURN
      END
!----------------------------------------------------------------------!

!----------------- 2. 全局弹性张量计算子程序ELASTO --------------------!
      SUBROUTINE ELASTO(EMATX0,T,EMATX)
      IMPLICIT NONE
C-----  Calculate the transformed E
      real(kind=8)::T1(6,6),ET(6,6),T1T(6,6),EMATX0(6,6),T(3,3),
     2              EMATX(6,6)
      REAL*8 TWO
      DATA TWO /2.0D0/
      REAL*8 T11,T12,T13,T21,T22,T23,T31,T32,T33
      INTEGER I,J
      T11=T(1,1)
      T12=T(1,2)
      T13=T(1,3)
      T21=T(2,1)
      T22=T(2,2)
      T23=T(2,3)
      T31=T(3,1)
      T32=T(3,2)
      T33=T(3,3)
      T1(1,1)=T11*T11
      T1(1,2)=T12*T12
      T1(1,3)=T13*T13
      T1(1,4)=T11*T12
      T1(1,5)=T11*T13
      T1(1,6)=T12*T13
      T1(2,1)=T21*T21
      T1(2,2)=T22*T22
      T1(2,3)=T23*T23
      T1(2,4)=T21*T22
      T1(2,5)=T21*T23
      T1(2,6)=T22*T23
      T1(3,1)=T31*T31
      T1(3,2)=T32*T32
      T1(3,3)=T33*T33
      T1(3,4)=T31*T32
      T1(3,5)=T31*T33
      T1(3,6)=T32*T33
      T1(4,1)=TWO*T11*T21
      T1(4,2)=TWO*T12*T22
      T1(4,3)=TWO*T13*T23
      T1(4,4)=T11*T22+T21*T12
      T1(4,5)=T13*T21+T11*T23
      T1(4,6)=T13*T22+T12*T23
      T1(5,1)=TWO*T11*T31
      T1(5,2)=TWO*T12*T32
      T1(5,3)=TWO*T13*T33
      T1(5,4)=T12*T31+T11*T32
      T1(5,5)=T13*T31+T11*T33
      T1(5,6)=T13*T32+T12*T33
      T1(6,1)=TWO*T21*T31
      T1(6,2)=TWO*T22*T32
      T1(6,3)=TWO*T23*T33
      T1(6,4)=T22*T31+T21*T32
      T1(6,5)=T23*T31+T21*T33
      T1(6,6)=T23*T32+T22*T33
      CALL MULTIP(EMATX0,6,6,T1,6,ET)
      DO 10 I=1,6
         DO 10 J=1,6
            T1T(J,I)=T1(I,J)
10    CONTINUE
      CALL MULTIP(T1T,6,6,ET,6,EMATX)
      RETURN
      END
!----------------------------------------------------------------------!

!-------------------- 3. 旋转矩阵计算子程序ELASTO ---------------------!
      SUBROUTINE ROTATION (PROP, ROTATE)
C-----  This subroutine calculates the rotation matrix, i.e. the
C     direction cosines of cubic crystal [100], [010] and [001]
C     directions in global system
C-----  The rotation matrix is stored in the array ROTATE.
C-----  Use single precision on cray
      IMPLICIT NONE
      REAL(KIND=8)::PROP(16), ROTATE(3,3), TERM1(3,3), TERM2(3,3)
      INTEGER::INDX(3)
      INTEGER  I,J,K
      REAL*8 DCMP, ANGLE1, ANGLE2
    
C-----  Subroutines:
C       CROSS  -- cross product of two vectors
C       LUDCMP -- LU decomposition
C       LUBKSB -- linear equation solver based on LU decomposition
C                 method (must call LUDCMP first)
C-----  PROP -- constants characterizing the crystal orientation
C            PROP(1) - PROP(3) -- direction of the first vector in
C                                 local cubic crystal system
C            PROP(4) - PROP(6) -- direction of the first vector in
C                                 global system
C            PROP(9) - PROP(11)-- direction of the second vector in
C                                 local cubic crystal system
C            PROP(12)- PROP(14)-- direction of the second vector in
C                                 global system
C-----  ROTATE -- rotation matrix (OUTPUT):
C            ROTATE(i,1) -- direction cosines of direction [1 0 0] in
C                           local cubic crystal system
C            ROTATE(i,2) -- direction cosines of direction [0 1 0] in
C                           local cubic crystal system
C            ROTATE(i,3) -- direction cosines of direction [0 0 1] in
C                           local cubic crystal system
C-----  local matrix: TERM1
      CALL CROSS (PROP(1), PROP(9), TERM1, ANGLE1)
C-----  LU decomposition of TERM1
      CALL LUDCMP (TERM1, 3, 3, INDX, DCMP)
C-----  inverse matrix of TERM1: TERM2
      DO J=1,3
         DO I=1,3
            IF (I.EQ.J) THEN
               TERM2(I,J)=1.0D0
            ELSE
               TERM2(I,J)=0.0D0
            END IF
         END DO
      END DO
      DO J=1,3
         CALL LUBKSB (TERM1, 3, 3, INDX, TERM2(1,J))
      END DO
C-----  global matrix: TERM1
      CALL CROSS (PROP(4), PROP(12), TERM1, ANGLE2)
C-----  Check: the angle between first and second vector in local and
C     global systems must be the same.  The relative difference must be
C     less than 0.1%.
      IF (DABS(ANGLE1/ANGLE2-1.0D0).GT.0.001D0) THEN
         WRITE (6,*)
     2      '***ERROR - angles between two vectors are not the same'
         STOP
      END IF
C-----  rotation matrix: ROTATE
      DO J=1,3
         DO I=1,3
            ROTATE(I,J)=0.0D0
            DO K=1,3
               ROTATE(I,J)=ROTATE(I,J)+TERM1(I,K)*TERM2(K,J)
            END DO
         END DO
      END DO
      RETURN
      END
!----------------------------------------------------------------------!
      
!---------------------- 4. 向量叉乘子程序CROSS ------------------------!
      SUBROUTINE CROSS (A, B, C, ANGLE)
C-----  (1) normalize vectors A and B to unit vectors
C       (2) store A, B and A*B (cross product) in C
           IMPLICIT NONE
           REAL(KIND=8)::A(3), B(3), C(3,3)
           REAL*8 SUM1, SUM2, ANGLE, SUM3
           INTEGER I,J
           SUM1=DSQRT(A(1)**2.0D0+A(2)**2.0D0+A(3)**2.0D0)
           SUM2=DSQRT(B(1)**2.0D0+B(2)**2.0D0+B(3)**2.0D0)
           IF (SUM1.EQ.0.0D0) THEN
              WRITE (6,*) '***ERROR - first vector is zero'
              STOP
           ELSE
              DO I=1,3
                 C(I,1)=A(I)/SUM1
              END DO
           END IF
           IF (SUM2.EQ.0.0D0) THEN
              WRITE (6,*) '***ERROR - second vector is zero'
              STOP
           ELSE
              DO I=1,3
                 C(I,2)=B(I)/SUM2
              END DO
           END IF
           ANGLE=0.0D0
           DO I=1,3
              ANGLE=ANGLE+C(I,1)*C(I,2)
           END DO
           ANGLE=DACOS(ANGLE)
           C(1,3)=C(2,1)*C(3,2)-C(3,1)*C(2,2)
           C(2,3)=C(3,1)*C(1,2)-C(1,1)*C(3,2)
           C(3,3)=C(1,1)*C(2,2)-C(2,1)*C(1,2)
           SUM3=DSQRT(C(1,3)**2.0D0+C(2,3)**2.0D0+C(3,3)**2.0D0)
           IF (SUM3.LT.1.0D-8) THEN
              WRITE (6,*)
     2           '***ERROR - first and second vectors are parallel'
               STOP
            END IF
           RETURN
           END
!----------------------------------------------------------------------!

!----------------------5. 塑性流动子程序STRAINRATE --------------------!
      SUBROUTINE STRAINRATE (TAUSLP, GSLIP, FSLIP, DFDXSP)
      USE PARABANK
      IMPLICIT NONE
      REAL(KIND=8),EXTERNAL::FITP, FITPDX
C-----  This subroutine calculates the slip rates FSLIP and the
C       derivatives of the slip rates with respect to the shear stresses
      INTEGER::I
      REAL(KIND=8)::TAUSLP(ND), GSLIP(ND)
      REAL(KIND=8)::FSLIP(ND), DFDXSP(ND)    
      REAL*8 X
      !! 正常滑移系的滑移率和导数计算
      DO I=1,NSLPTL
        X=TAUSLP(I)/GSLIP(I)
        FSLIP(I)=FITP(X)
        DFDXSP(I)=FITPDX(X)
      END DO
      !! 孪晶滑移系的滑移率和导数计算
      DO I = NSLPTL+1, ND
        IF (TAUSLP(I).GT.0.0D0) THEN
          X=TAUSLP(I)/GSLIP(I)
        ELSE
          X=0.0D0 !!孪晶系必须沿激活方向加载才会滑移
        END IF
        FSLIP(I)=FITP(X)
        DFDXSP(I)=FITPDX(X)
      END DO
      !!需要修改
      RETURN
      END

C-----  Interplotation function
      REAL*8 FUNCTION FITP(TEMPY)
        USE PARABANK
        IMPLICIT NONE
        REAL*8 TEMPY
        FITP=GAMDOT0*(ABS(TEMPY))**qEXP*DSIGN(1.D0,TEMPY)!!检查一下有没有ABS
      RETURN
      END
           
C-----  Interplotation function
      REAL*8 FUNCTION FITPDX(TEMPY)
        USE PARABANK
        IMPLICIT NONE
        REAL*8 TEMPY
        FITPDX=GAMDOT0*qEXP*(ABS(TEMPY))**(qEXP-1.0D0)!!检查一下有没有ABS
      RETURN
      END
!----------------------------------------------------------------------!

!------------------ 6. LU分解与线性方程求解子程序 ---------------------!          
      SUBROUTINE LUDCMP (A, N, NP, INDX, D)
C-----  LU decomposition
c-----  A, output
C-----  Use single precision on cray
      IMPLICIT NONE
      INTEGER NMAX, I, J, K, N, NP, IMAX
      PARAMETER (NMAX=200)
      REAL*8 TINY
      PARAMETER (TINY=1.0D-20)
      REAL*8 AAMAX, D, SUM, DUM
      REAL(KIND=8)::A(NP,NP), VV(NMAX)
      INTEGER::INDX(N)
      D=1.0D0
      DO I=1,N
         AAMAX=0.0D0
         DO J=1,N
            IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))
         END DO
         IF (AAMAX.EQ.0.0D0) THEN 
            WRITE(*,*) 'Singular matrix.'
            AAMAX = TINY
         END IF
         !STOP 'ERROR: Singular matrix.'
         VV(I)=1.0D0/AAMAX
      END DO
      DO J=1,N
         DO I=1,J-1
            SUM=A(I,J)
            DO K=1,I-1
               SUM=SUM-A(I,K)*A(K,J)
            END DO
            A(I,J)=SUM
         END DO
         AAMAX=0.0D0
         DO I=J,N
            SUM=A(I,J)
            DO K=1,J-1
               SUM=SUM-A(I,K)*A(K,J)
            END DO
            A(I,J)=SUM
            DUM=VV(I)*DABS(SUM)
            IF (DUM.GE.AAMAX) THEN
               IMAX=I
               AAMAX=DUM
            END IF
         END DO
         IF (J.NE.IMAX) THEN
            DO K=1,N
               DUM=A(IMAX,K)
               A(IMAX,K)=A(J,K)
               A(J,K)=DUM
            END DO
            D=-D
            VV(IMAX)=VV(J)
         END IF
         INDX(J)=IMAX
         IF (A(J,J).EQ.0.0D0) A(J,J)=TINY
         IF (J.NE.N) THEN
            DUM=1.0D0/A(J,J)
            DO I=J+1,N
               A(I,J)=A(I,J)*DUM
            END DO
         END IF
      END DO
      RETURN
      END
!----------------------------------------------------------------------!

!-------------------- 7. 线性方程求解子程序LUBKSB ---------------------!
      SUBROUTINE LUBKSB (A, N, NP, INDX, B)
C-----  Linear equation solver based on LU decomposition
C-----  A, input; B, output
C-----  Use single precision on cray
      IMPLICIT NONE
      INTEGER II, I, LL, J, N, NP
      REAL*8 SUM
      REAL(KIND=8)::A(NP,NP),B(N)
      INTEGER::INDX(N)
      II=0
      DO I=1,N
         LL=INDX(I)
         SUM=B(LL)
         B(LL)=B(I)
         IF (II.NE.0) THEN
            DO J=II,I-1
               SUM=SUM-A(I,J)*B(J)
            END DO
         ELSE IF (SUM.NE.0.0D0) THEN
            II=I
         END IF
         B(I)=SUM
      END DO
      DO I=N,1,-1
         SUM=B(I)
         IF (I.LT.N) THEN
            DO J=I+1,N
               SUM=SUM-A(I,J)*B(J)
            END DO
         END IF
         B(I)=SUM/A(I,I)
      END DO
      RETURN
      END
!----------------------------------------------------------------------!
!***********************************************************************

!     设置晶粒晶格点整为HCP点阵
      INCLUDE 'HCP_slip_twin.for'
!     设置了NSLIP: 每个滑移族下的滑移系
!           NSLPTL: 滑移系的总数
!           NTWTL: 孪晶系的总数
!           NSET: 滑移族
!           SIGMA0: 每个滑移系族的热相关阻力
!           BURGER: 每个滑移系族的burgers矢量REAL
!           SLPDIR：每个滑移系全局向量
!           SLPNOR：每个滑移系全局法向
      INCLUDE 'DIS_IRRA_HARD.for'
