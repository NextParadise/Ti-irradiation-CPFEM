# Ti-IHCPFEM 参数拟合指南

## 目录
1. [模型概述](#模型概述)
2. [参数分类](#参数分类)
3. [参数对曲线形状的影响](#参数对曲线形状的影响)
4. [拟合策略](#拟合策略)
5. [脚本使用说明](#脚本使用说明)

---

## 模型概述

本模型是钛合金辐照硬化晶体塑性有限元（CPFEM）本构模型，包含以下物理机制：

- **晶体塑性变形**：30个滑移系 + 12个孪晶系
- **位错硬化**：基于位错密度演化的Taylor硬化
- **孪晶变形**：拉伸孪晶(TT)和压缩孪晶(CT)
- **孔隙率演化**：形核和生长两阶段模型
- **辐照硬化**：（预留接口）

---

## 参数分类

### 1. 弹性参数（INP文件中定义）

| 参数 | 符号 | 典型值 | 说明 |
|------|------|--------|------|
| C11 | PROPS(1) | 162.4 GPa | HCP弹性常数 |
| C12 | PROPS(2) | 92.0 GPa | HCP弹性常数 |
| C13 | PROPS(3) | 69.0 GPa | HCP弹性常数 |
| C33 | PROPS(4) | 180.7 GPa | HCP弹性常数 |
| C44 | PROPS(5) | 46.7 GPa | HCP弹性常数 |

### 2. 滑移系参数（FIT.for中定义）

| 参数 | 符号 | 当前值 | 单位 | 说明 |
|------|------|--------|------|------|
| 剪切模量 | GMOD | 40067 | MPa | 用于位错硬化计算 |
| 晶粒尺寸 | DGRAIN | 14e-3 | mm | Hall-Petch效应 |

### 3. 位错演化参数

| 参数 | 符号 | 当前值 | 说明 |
|------|------|--------|------|
| 初始位错密度 | DISLOCATION_DENSITY | 2.7e6 mm⁻² | 各滑移系初始值 |
| 位错增殖系数 | PMUL | 0.3 | 控制加工硬化速率 |
| 位错湮灭系数 | PDYN | 80.0 | 控制动态回复 |
| 位错硬化系数 | HRHO | 1.08 | Taylor硬化强度 |
| 派纳力参数 | MPNS | 1.0 | 热激活阻力 |

### 4. 孪晶参数

| 参数 | 符号 | 当前值 | 单位 | 说明 |
|------|------|--------|------|------|
| TT初始切应力 | TAU0_TWIN(1) | 389.5 | kPa | 拉伸孪晶临界应力 |
| CT初始切应力 | TAU0_TWIN(2) | 391.5 | kPa | 压缩孪晶临界应力 |
| 孪晶-滑移硬化 | H_SL_TW | 240.0 | - | 孪晶对滑移的硬化 |
| 孪晶-滑移指数 | TWIN_D | 3.0 | - | 硬化指数 |
| 孪晶-孪晶硬化 | H_TW_TW | 350.0 | - | 孪晶间交互硬化 |
| 孪晶-孪晶指数 | TWIN_B | 2.3 | - | 硬化指数 |

### 5. 孔隙率参数

| 参数 | 符号 | 当前值 | 说明 |
|------|------|--------|------|
| 初始孔隙率 | POROS0 | 0.00112 | 初始孔隙体积分数 |
| 形核应变均值 | EVALEPN | 0.21 | Gaussian分布均值 |
| 形核应变标准差 | SEVLEPN | 0.036 | Gaussian分布标准差 |
| 夹杂物体积分数 | INCLUSF | 0.17 | 形核位点密度 |
| 孔隙率参数A | PORO_A | 60.0 | 软化函数参数 |
| 孔隙率参数B | PORO_B | 0.173 | 生长速率参数 |

### 6. 流动法则参数（PARABANK模块）

| 参数 | 符号 | 当前值 | 说明 |
|------|------|--------|------|
| 参考滑移率 | GAMDOT0 | 0.001 s⁻¹ | 流动法则参数 |
| 应变率敏感指数 | qEXP | 20.0 | 流动法则指数 |

---

## 参数对曲线形状的影响

### 1. 应力-应变曲线

```
应力 σ
  ↑
  │         ╭────────── 饱和应力（受PDYN控制）
  │       ╱
  │     ╱  ← 加工硬化斜率（受PMUL, HRHO控制）
  │   ╱
  │ ╱ ← 屈服点（受SIGMA0, DISLOCATION_DENSITY控制）
  │╱
  └──────────────────→ 应变 ε
```

| 参数 | 增大时的效果 |
|------|-------------|
| **SIGMA0** | 屈服应力↑，整体曲线上移 |
| **DISLOCATION_DENSITY** | 初始屈服应力↑（Taylor硬化） |
| **PMUL** | 加工硬化速率↑，曲线斜率↑ |
| **PDYN** | 动态回复↑，饱和应力↓，曲线更早趋于平缓 |
| **HRHO** | 位错硬化强度↑，整体硬化效果↑ |
| **qEXP** | 应变率敏感性↓，流动应力对应变率变化更不敏感 |

#### 调参建议：
1. **屈服点偏低**：增大SIGMA0或DISLOCATION_DENSITY
2. **硬化不足**：增大PMUL或HRHO
3. **饱和应力偏高**：增大PDYN
4. **曲线形状不对**：调整PMUL/PDYN比值

### 2. 孪晶体积分数曲线

```
孪晶体积分数 f_tw
  ↑
  │              ╭────── 饱和值（受H_TW_TW控制）
  │            ╱
  │          ╱  ← 增长速率（受TAU0_TWIN, H_SL_TW控制）
  │        ╱
  │      ╱
  │    ╱
  │  ╱ ← 起始点（受TAU0_TWIN控制）
  │╱
  └──────────────────→ 应变 ε
```

| 参数 | 增大时的效果 |
|------|-------------|
| **TAU0_TWIN** | 孪晶启动应变↑，曲线右移 |
| **H_SL_TW** | 滑移对孪晶的抑制↑，孪晶增长变慢 |
| **TWIN_D** | 抑制效果的非线性↑ |
| **H_TW_TW** | 孪晶间交互硬化↑，饱和值↓ |
| **TWIN_B** | 饱和效果的非线性↑ |

#### 调参建议：
1. **孪晶启动太早**：增大TAU0_TWIN
2. **孪晶增长太快**：增大H_SL_TW或H_TW_TW
3. **孪晶饱和值太高**：增大H_TW_TW

### 3. 孔隙率演化曲线

```
孔隙率 f
  ↑
  │                    ╭── 快速增长（生长主导）
  │                  ╱
  │                ╱
  │              ╱
  │            ╱ ← 形核峰值（受EVALEPN, SEVLEPN控制）
  │          ╱
  │        ╱
  │──────╱ ← 初始孔隙率（POROS0）
  └──────────────────→ 应变 ε
```

| 参数 | 增大时的效果 |
|------|-------------|
| **POROS0** | 初始孔隙率↑，曲线整体上移 |
| **EVALEPN** | 形核峰值应变↑，形核延迟 |
| **SEVLEPN** | 形核分布更宽，峰值更平缓 |
| **INCLUSF** | 形核位点↑，形核速率↑ |
| **PORO_A** | 软化效果↑，应力下降更快 |
| **PORO_B** | 孔隙生长速率↑ |

#### 调参建议：
1. **孔隙率增长太早**：增大EVALEPN
2. **孔隙率增长太快**：减小INCLUSF或PORO_B
3. **软化效果不明显**：增大PORO_A

### 4. 位错密度演化曲线

```
位错密度 ρ
  ↑
  │         ╭────────── 饱和密度 ρ_sat = (PMUL/PDYN)²
  │       ╱
  │     ╱  ← 增长速率（受PMUL控制）
  │   ╱
  │ ╱
  │╱ ← 初始密度（DISLOCATION_DENSITY）
  └──────────────────→ 应变 ε
```

位错密度演化方程：
```
dρ/dγ = PMUL × √ρ - PDYN × ρ
```

饱和密度：
```
ρ_sat = (PMUL / PDYN)²
```

| 参数 | 增大时的效果 |
|------|-------------|
| **DISLOCATION_DENSITY** | 初始密度↑ |
| **PMUL** | 增殖速率↑，饱和密度↑ |
| **PDYN** | 湮灭速率↑，饱和密度↓ |

---

## 拟合策略

### 推荐拟合顺序

```
第1步：弹性参数
    ↓
第2步：屈服应力相关参数（SIGMA0, DISLOCATION_DENSITY）
    ↓
第3步：加工硬化参数（PMUL, PDYN, HRHO）
    ↓
第4步：孪晶参数（TAU0_TWIN, H_SL_TW, H_TW_TW）
    ↓
第5步：孔隙率参数（如需要）
```

### 第1步：弹性参数

**目标**：拟合弹性段斜率

**方法**：
- 从文献获取钛合金弹性常数
- 或从实验曲线弹性段计算杨氏模量

### 第2步：屈服应力

**目标**：拟合屈服点位置

**敏感参数**：
- SIGMA0（各滑移系派纳力）
- DISLOCATION_DENSITY（初始位错密度）

**方法**：
1. 固定DISLOCATION_DENSITY为文献值
2. 调整SIGMA0使屈服点匹配

### 第3步：加工硬化

**目标**：拟合塑性段斜率和饱和行为

**敏感参数**：
- PMUL（位错增殖）
- PDYN（位错湮灭）
- HRHO（硬化系数）

**方法**：
1. 先调PMUL/PDYN比值控制饱和应力
2. 再调PMUL绝对值控制硬化速率
3. 微调HRHO优化曲线形状

### 第4步：孪晶参数

**目标**：拟合孪晶体积分数演化

**敏感参数**：
- TAU0_TWIN（孪晶临界应力）
- H_SL_TW, H_TW_TW（交互硬化）

**方法**：
1. 调TAU0_TWIN控制孪晶启动点
2. 调H_SL_TW控制增长速率
3. 调H_TW_TW控制饱和值

### 第5步：孔隙率参数

**目标**：拟合孔隙率演化和软化行为

**敏感参数**：
- EVALEPN, SEVLEPN（形核参数）
- PORO_A, PORO_B（生长参数）

---

## 脚本使用说明

### 1. 提交计算

```bash
# 基本使用
python submit_abaqus.py

# 指定参数
python submit_abaqus.py -j MyJob -c 8 -m 8gb

# 仅检查文件
python submit_abaqus.py --check

# 预览命令（不执行）
python submit_abaqus.py --dry-run

# 创建新工作目录
python submit_abaqus.py --new-dir
```

### 2. 后处理

```bash
# 在Abaqus Python环境中提取数据
abaqus python postprocess.py --odb Ti13GN.odb

# 仅提取数据，不绘图
abaqus python postprocess.py --odb Ti13GN.odb --extract-only

# 使用已提取的数据绘图
python postprocess.py --plot-only --data results/extracted_data.npz

# 与实验数据对比
python postprocess.py --plot-only --data results/extracted_data.npz --exp-data experiment.csv
```

### 3. 参数扫描示例

```python
from submit_abaqus import submit_batch_jobs, DEFAULT_CONFIG

# 定义参数扫描
param_sets = [
    ('PMUL_0.2', {'PMUL': 0.2}),
    ('PMUL_0.3', {'PMUL': 0.3}),
    ('PMUL_0.4', {'PMUL': 0.4}),
]

# 批量提交
results = submit_batch_jobs(
    param_sets,
    DEFAULT_CONFIG,
    src_dir='/path/to/source',
    output_base_dir='/path/to/output'
)
```

### 4. 输出文件说明

| 文件 | 说明 |
|------|------|
| `stress_strain.png` | 应力-应变曲线 |
| `twin_evolution.png` | 孪晶体积分数演化 |
| `porosity_evolution.png` | 孔隙率演化 |
| `dislocation_evolution.png` | 位错密度演化 |
| `slip_activity.png` | 各滑移系活动分析 |
| `comprehensive.png` | 综合分析图 |
| `stress_strain.csv` | 应力-应变数据 |
| `twin_evolution.csv` | 孪晶数据 |
| `porosity_evolution.csv` | 孔隙率数据 |
| `dislocation_evolution.csv` | 位错密度数据 |
| `extracted_data.npz` | 完整提取数据 |

---

## 常见问题

### Q1: 计算不收敛怎么办？

**可能原因**：
1. 时间步长太大
2. 参数设置不合理（如SIGMA0太小）
3. 数值奇异（如除零）

**解决方法**：
1. 减小INP文件中的时间增量
2. 检查参数是否在合理范围
3. 检查代码中的除零保护

### Q2: 应力-应变曲线形状不对？

**检查顺序**：
1. 弹性段斜率 → 弹性常数
2. 屈服点 → SIGMA0, DISLOCATION_DENSITY
3. 硬化斜率 → PMUL, HRHO
4. 饱和行为 → PDYN

### Q3: 孪晶体积分数太高/太低？

**调整参数**：
- 太高：增大H_TW_TW或TAU0_TWIN
- 太低：减小H_TW_TW或TAU0_TWIN
- 增长太快：增大H_SL_TW

---

## 参考文献

1. Huang, Y. (1991). A user-material subroutine incorporating single crystal plasticity in the ABAQUS finite element program.
2. Kalidindi, S.R. (1998). Incorporation of deformation twinning in crystal plasticity models.
3. Gurson, A.L. (1977). Continuum theory of ductile rupture by void nucleation and growth.

---

*文档版本: v1.0*
*最后更新: 2026-01-23*
