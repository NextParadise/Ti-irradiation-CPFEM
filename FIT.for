      HRHO=1.0D0
      GMOD=40067.0D0 !Ti的剪切模量，单位MPa
      DGRAIN=14.0D-3 !晶粒尺寸，单位mm

      SVOID=0.0D0
      SLOOP=0.0D0
      SPRE=0.0D0
      HVOID=1.0D0
      HLOOP=1.0D0
      HPRE=1.0D0

C     给每个滑移系和他们的初始位错密度赋值(单位mm^-2)
      DISLOCATION_DENSITY(1)=2.7D6 !Basal
      DISLOCATION_DENSITY(2)=2.7D6 !Prism
      DISLOCATION_DENSITY(3)=2.7D6 !Pyram
      DISLOCATION_DENSITY(4)=2.7D6 !Pyram I
      DISLOCATION_DENSITY(5)=2.7D6 !Pyram II

      PMUL=0.3d0 !位错增殖参数
      PDYN=80.0D0 !位错湮灭参数
      HRHO=1.08D0  !位错
      MPNS=1.0D0  !派纳力参数

C     孪晶参数
      TAU0_TWIN(1) = 389.5D3 !TT
      TAU0_TWIN(2) = 391.5D3 !CT
      
      H_SL_TW=240.0D0                  ! 孪晶对滑移交互硬化系数
      TWIN_D=3.0D0                     ! 孪晶对滑移交互硬化系数
      H_TW_TW=350.0D0                  ! 孪晶对孪晶交互硬化系数
      TWIN_B=2.3D0                     ! 孪晶对孪晶交互硬化系数
      

