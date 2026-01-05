# -*- coding: utf-8 -*-
from abaqusConstants import *
from odbAccess import *
import numpy as np
import sys

# ==============================================================================
#                                 配N置区域
# ==============================================================================

# 1. ODB 文件名 (不带 .odb)
odb_name = 'Ti13GN' 

# 2. 关注的方向 (1=X, 2=Y, 3=Z)
#    你指定了 Z 方向
target_axis = 3 

# 3. 输出单位转换
#    假设模型单位为 (mm, N) -> 应力直接是 MPa
#    假设模型单位为 (m, N)  -> 应力是 Pa，需要除以 1e6 变 MPa
#    请根据你的实际情况修改 scale_factor
scale_factor = 1.0  # 如果是 Pa 转 MPa，填 1e-6；如果是 MPa 直出，填 1.0

# ==============================================================================
#                                 脚本逻辑
# ==============================================================================

def calculate_volume_average():
    print("Opening ODB: %s.odb..." % odb_name)
    try:
        odb = openOdb(path=odb_name + '.odb')
    except:
        print("Error: ODB file not found.")
        return

    # 获取最后一个分析步
    step_keys = odb.steps.keys()
    last_step = odb.steps[step_keys[-1]]
    print("Processing Step: %s" % last_step.name)

    # 准备存储列表
    macro_strain_list = []
    macro_stress_list = []
    frame_times = []

    # 获取 Instance (假设只有一个)
    instance_name = odb.rootAssembly.instances.keys()[0]
    instance = odb.rootAssembly.instances[instance_name]
    
    # -------------------------------------------------------
    # 预计算单元体积 (IVOL)
    # -------------------------------------------------------
    # 注意：对于小变形，体积变化不大，可以用第一帧的体积
    # 对于大变形，应该每帧更新体积。这里为了精确，尝试每帧读取 IVOL。
    # 如果你的输出里没有 IVOL，我们尝试用 EVOL 或者近似计算。
    
    print("Starting frame processing...")
    
    # 对应 Python 索引 (Z方向: 3 -> index 2)
    idx = target_axis - 1
    
    for i, frame in enumerate(last_step.frames):
        sys.stdout.write("\rProcessing Frame %d / %d" % (i, len(last_step.frames)))
        sys.stdout.flush()
        
        frame_time = frame.frameValue
        frame_times.append(frame_time)

        # 1. 获取场变量
        # S: 应力张量, LE: 对数应变张量, IVOL: 积分点体积
        # 注意：必须在 INP 中输出 IVOL (Element Output)
        try:
            field_S = frame.fieldOutputs['S']
            field_LE = frame.fieldOutputs['LE']
        except KeyError:
            print("\nError: S or LE not found in field outputs.")
            return
            
        # 尝试获取体积数据
        # 优先找 IVOL (Integration Point Volume)，如果没有则找 EVOL (Element Volume)
        has_ivol = 'IVOL' in frame.fieldOutputs.keys()
        has_evol = 'EVOL' in frame.fieldOutputs.keys()
        
        if has_ivol:
            field_V = frame.fieldOutputs['IVOL']
        elif has_evol:
            field_V = frame.fieldOutputs['EVOL']
        else:
            # 如果都没有体积输出，这是一个严重问题，无法加权
            # 此时只能假设体积均匀（通常不正确）或报错
            if i == 0: print("\nWarning: No IVOL/EVOL found. Assuming unit volume (Results will be average of values, not volume-weighted).")
            field_V = None

        # 2. 提取数据 (获取所有单元的数据)
        # values 是一个列表，包含所有积分点的值
        s_vals = field_S.values
        le_vals = field_LE.values
        
        if field_V:
            v_vals = field_V.values
        
        # 3. 执行体积加权平均
        # Sum(S_i * V_i) / Sum(V_i)
        
        sum_stress_vol = 0.0
        sum_strain_vol = 0.0
        total_vol = 0.0
        
        # 这里的 len(s_vals) 是所有积分点的总数
        # 为了效率，使用 numpy 数组操作
        
        # 提取标量分量 (S33 和 LE33)
        # Abaqus S 数据顺序: S11, S22, S33, S12, S13, S23
        # Z方向正应力是 S33 (index 2)
        
        # 转换为 Numpy 数组加速计算
        # data 属性返回 float32 数组
        s_data = np.array([v.data[idx] for v in s_vals])
        le_data = np.array([v.data[idx] for v in le_vals])
        
        if field_V:
            # IVOL 是标量
            vol_data = np.array([v.data for v in v_vals])
        else:
            vol_data = np.ones(len(s_data))

        # 核心计算
        total_vol = np.sum(vol_data)
        sum_stress_vol = np.sum(s_data * vol_data)
        sum_strain_vol = np.sum(le_data * vol_data)
        
        # 平均值
        avg_stress = sum_stress_vol / total_vol
        avg_strain = sum_strain_vol / total_vol
        
        macro_stress_list.append(avg_stress)
        macro_strain_list.append(avg_strain)

    print("\nProcessing complete.")

    # 4. 数据导出与绘图
    stress_arr = np.array(macro_stress_list) * scale_factor
    strain_arr = np.array(macro_strain_list)
    
    # 保存 CSV
    data = np.column_stack((frame_times, strain_arr, stress_arr))
    np.savetxt('vol_avg_stress_strain_Z.csv', data, 
               delimiter=',', 
               header='Time,LE33_Avg,S33_Avg(MPa)', 
               comments='')
    print("Data saved to vol_avg_stress_strain_Z.csv")

    # 简单的 matplotlib 绘图 (需要 matplotlib 库)
    try:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(strain_arr, stress_arr, 'r-', linewidth=2)
        plt.xlabel('True Strain (LE33)')
        plt.ylabel('True Stress (S33) [MPa]')
        plt.title('Volume Averaged Stress-Strain (Z-Direction)')
        plt.grid(True)
        plt.savefig('vol_avg_plot_Z.png')
        print("Plot saved to vol_avg_plot_Z.png")
    except ImportError:
        print("Matplotlib not found, skipping plot.")

    odb.close()

if __name__ == '__main__':
    calculate_volume_average()
