C     ==========================================================
C     根据单元编号 (NOEL) 判断所属晶粒 (IGRAIN) 以及是否为层状组织 (LAYERSTRUCTURE )
C     数据来源：基于用户提供的 INP 集合范围
C     ==========================================================

      IF (NOEL .LE. 315) THEN
          IGRAIN = 1
          LAYERSTRUCTURE = 0
      ELSE IF (NOEL .LE. 570) THEN
          IGRAIN = 2
          LAYERSTRUCTURE = 0
      ELSE IF (NOEL .LE. 814) THEN
          IGRAIN = 3
          LAYERSTRUCTURE = 0
      ELSE IF (NOEL .LE. 1205) THEN
          IGRAIN = 4
          LAYERSTRUCTURE = 0
      ELSE IF (NOEL .LE. 1469) THEN
          IGRAIN = 5
          LAYERSTRUCTURE = 0          
      ELSE IF (NOEL .LE. 1841) THEN
          IGRAIN = 6
          LAYERSTRUCTURE = 0
      ELSE IF (NOEL .LE. 2236) THEN
          IGRAIN = 7
          LAYERSTRUCTURE = 0
      ELSE IF (NOEL .LE. 2588) THEN
          IGRAIN = 8
          LAYERSTRUCTURE = 0
      ELSE IF (NOEL .LE. 2856) THEN
          IGRAIN = 9
          LAYERSTRUCTURE = 0
      ELSE IF (NOEL .LE. 3216) THEN
          IGRAIN = 10
          LAYERSTRUCTURE = 0
      ELSE IF (NOEL .LE. 3477) THEN
          IGRAIN = 11
          LAYERSTRUCTURE = 0
      ELSE IF (NOEL .LE. 3843) THEN
          IGRAIN = 12
          LAYERSTRUCTURE = 0
      ELSE !剩下的就是最后一个晶粒
          IGRAIN = 13
          LAYERSTRUCTURE = 0
      END IF