import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker  # [新增] 引入 ticker 模块

# --- 处理实验数据 ---
# 读取实验数据
df_exp = pd.read_csv('Default Dataset.csv')

true_strain = df_exp['True_Strain']
true_stress = df_exp['True_Stress']
# --- 处理模拟数据 ---
# 读取模拟数据
df_sim = pd.read_csv('vol_avg_stress_strain_Z.csv')

# 模拟数据列名通常包含空格，先清理一下
df_sim.columns = df_sim.columns.str.strip()

# 提取真应变和真应力
# 假设 LE33_Avg 是真应变，S33_Avg(MPa) 是真应力
sim_strain = df_sim['LE33_Avg']
sim_stress = df_sim['S33_Avg(MPa)']

# --- 绘图 ---
plt.figure(figsize=(10, 7))

# 绘制实验数据 (红色圆点)
plt.plot(true_strain, true_stress, 'o-', color='#D32F2F', label='Experiment', markersize=8, alpha=0.8, linewidth=3)

# 绘制模拟数据 (蓝色实线)
plt.plot(sim_strain, sim_stress, '-', color='#1976D2', label='Simulation', linewidth=5)

# 设置图表细节
plt.title('Experiment vs Simulation', fontsize=30)
plt.xlabel('True Strain', fontsize=30)
plt.ylabel('True Stress (MPa)', fontsize=30)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=30)

plt.tick_params(axis='both', which='major', labelsize=24)
# 设置坐标轴范围，使其更紧凑美观
plt.xlim(0, 0.1)
plt.ylim(0, max(df_exp['True_Stress'].max(), sim_stress.max()) * 1.05)
# 2. 处理原点 "0" 重叠问题
# 获取当前坐标轴对象
ax = plt.gca()

# 定义 X 轴格式化函数：如果是 0，返回空字符串；否则保留两位小数
# 这样原点处就只剩下 Y 轴的 "0" 了
def x_formatter(x, pos):
    if x == 0:
        return ""
    return f"{x:.2f}"

# 应用格式化器到 X 轴
ax.xaxis.set_major_formatter(ticker.FuncFormatter(x_formatter))
plt.show()
