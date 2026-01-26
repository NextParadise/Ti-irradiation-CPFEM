# 非局部化GTN孔隙率模型说明文档

## 1. 背景与动机

### 1.1 局部化问题

在传统的局部本构模型中，材料点的行为仅取决于该点自身的状态变量。当材料发生软化（如孔隙率增加导致的强度下降）时，这种局部模型会导致严重的**应变局部化**问题：

1. 某些晶粒因取向"软"，先达到临界应变，孔隙率快速增长
2. 这些晶粒软化后，变形进一步集中
3. 形成正反馈：软化 → 变形集中 → 更快软化 → 数值不收敛
4. 其他晶粒还未开始软化，整体曲线看不到软化效果

### 1.2 非局部化方法

非局部化方法的核心思想是：**材料的损伤/孔隙演化不是孤立发生的，而是受到周围材料约束的**。

这反映了真实的物理现象：
- 位错的长程应力场影响邻近区域
- 晶界的约束作用限制单个晶粒的自由变形
- 应变梯度效应产生额外的几何必需位错

---

## 2. 理论基础

### 2.1 积分型非局部模型

经典的积分型非局部模型将局部变量替换为空间加权平均：

$$
\bar{\varepsilon}^{nl}(\mathbf{x}) = \frac{\int_\Omega w(|\mathbf{x}-\boldsymbol{\xi}|) \bar{\varepsilon}(\boldsymbol{\xi}) d\boldsymbol{\xi}}{\int_\Omega w(|\mathbf{x}-\boldsymbol{\xi}|) d\boldsymbol{\xi}}
$$

其中权重函数通常取高斯型：

$$
w(r) = \exp\left(-\frac{r^2}{2l_c^2}\right)
$$

$l_c$ 是**特征长度**，代表材料的内禀尺度。

### 2.2 隐式梯度型非局部模型

隐式梯度方法将积分型非局部平均转化为求解Helmholtz型偏微分方程：

$$
\bar{\varepsilon}^{nl} - l_c^2 \nabla^2 \bar{\varepsilon}^{nl} = \bar{\varepsilon}^{local}
$$

其中：
- $\bar{\varepsilon}^{local}$ 是局部等效应变
- $\bar{\varepsilon}^{nl}$ 是非局部等效应变
- $l_c$ 是特征长度
- $\nabla^2$ 是拉普拉斯算子

**数学等价性**：在无限域中，该微分方程的解等价于积分型非局部平均。

### 2.3 晶粒尺度离散化

在多晶体材料中，晶界是天然的"屏障"，损伤的非局部影响主要在**晶粒内部**传递。因此，我们采用晶粒平均作为非局部值的近似：

$$
\bar{\varepsilon}^{nl} \approx \bar{\varepsilon}^{grain} = \frac{1}{N_g} \sum_{i=1}^{N_g} \bar{\varepsilon}_i
$$

其中 $N_g$ 是晶粒内的积分点数量。

**物理解释**：当特征长度 $l_c$ 约等于晶粒尺寸，且在晶界处满足自然边界条件 $\nabla\bar{\varepsilon}^{nl} = 0$ 时，隐式梯度方程的解趋于晶粒内均匀化，即晶粒平均值。

---

## 3. 孔隙率演化模型

### 3.1 GTN模型基本形式

孔隙率演化由形核和生长两部分组成：

$$
\dot{f} = \dot{f}_{nucleation} + \dot{f}_{growth}
$$

#### 形核项（Gaussian分布）

$$
\dot{f}_{nucleation} = \frac{f_I}{\sqrt{2\pi}\sigma_N} \exp\left[-\frac{1}{2}\left(\frac{\bar{\varepsilon}^p - \varepsilon_N}{\sigma_N}\right)^2\right] \dot{\bar{\varepsilon}}^p
$$

其中：
- $f_I$ = `INCLUSF`：夹杂物体积分数
- $\varepsilon_N$ = `EVALEPN`：形核应变均值
- $\sigma_N$ = `SEVLEPN`：形核应变标准差

#### 生长项

$$
\dot{f}_{growth} = B \sinh\left(\frac{n-0.5}{n+0.5} \cdot \frac{\sigma_H}{\bar{\sigma}}\right) \left[(1-f)^{-n} - (1-f)\right] \dot{\bar{\varepsilon}}^p
$$

其中：
- $B$ = `PORO_B`：生长速率参数
- $n$ = `PORO_N`：应变率敏感指数
- $\sigma_H/\bar{\sigma}$：应力三轴度

### 3.2 非局部化修改

**核心修改**：形核项中的等效塑性应变使用非局部值（晶粒平均）：

$$
\dot{f}_{nucleation} = \frac{f_I}{\sqrt{2\pi}\sigma_N} \exp\left[-\frac{1}{2}\left(\frac{\textcolor{red}{\bar{\varepsilon}^{nl}} - \varepsilon_N}{\sigma_N}\right)^2\right] \dot{\bar{\varepsilon}}^p
$$

**为什么只对形核项非局部化？**

1. 形核项包含高斯峰值函数，对等效塑性应变非常敏感
2. 当某个积分点的应变先达到 $\varepsilon_N$ 时，形核率会突然增大，导致局部化
3. 生长项相对平缓，主要受应力三轴度控制，局部化倾向较弱

---

## 4. 实现细节

### 4.1 数据结构

在 `PARABANK` 模块中添加：

```fortran
REAL*8 EQVPL_SUM(GRNUM)     ! 各晶粒内等效塑性应变累加和
REAL*8 EQVPL_AVG(GRNUM)     ! 各晶粒的平均等效塑性应变（非局部值）
INTEGER GRAIN_NPTS(GRNUM)   ! 各晶粒内积分点计数
INTEGER KINC_PREV           ! 上一增量步编号
```

### 4.2 计算流程

```
┌─────────────────────────────────────────────────────────────┐
│  每个增量步开始时（KINC ≠ KINC_PREV）：                      │
│  1. 计算上一步各晶粒的平均等效塑性应变                        │
│     EQVPL_AVG(I) = EQVPL_SUM(I) / GRAIN_NPTS(I)            │
│  2. 重置累加器                                              │
│     EQVPL_SUM(I) = 0, GRAIN_NPTS(I) = 0                    │
│  3. 更新 KINC_PREV = KINC                                   │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│  每个积分点计算时：                                          │
│  1. 获取非局部等效塑性应变                                   │
│     EQVPL_NL = EQVPL_AVG(IGRAIN)                           │
│  2. 调用孔隙率演化（使用非局部值）                           │
│     CALL PORO_EVOL(..., EQVPL_NL, ...)                     │
│  3. 累加当前积分点的EQVPL到所属晶粒                          │
│     EQVPL_SUM(IGRAIN) += EQVPL                             │
│     GRAIN_NPTS(IGRAIN) += 1                                │
└─────────────────────────────────────────────────────────────┘
```

### 4.3 时间延迟

由于使用上一增量步的平均值，存在一个增量步的时间延迟。这在实际计算中影响很小，因为：
1. 增量步通常很小
2. 等效塑性应变是连续变化的
3. 延迟一步相当于隐式时间积分的一种形式

---

## 5. 参数说明

### 5.1 非局部化相关参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `GRNUM` | 晶粒数量 | 13 |
| 特征长度 | 隐含为晶粒尺寸 | ~14 μm |

### 5.2 孔隙率本构参数

| 参数 | 符号 | 说明 | 典型值 |
|------|------|------|--------|
| `POROS0` | $f_0$ | 初始孔隙率 | 0.00112 |
| `EVALEPN` | $\varepsilon_N$ | 形核应变均值 | 0.21 |
| `SEVLEPN` | $\sigma_N$ | 形核应变标准差 | 0.036 |
| `INCLUSF` | $f_I$ | 夹杂物体积分数 | 0.17 |
| `PORO_A` | $A$ | 软化函数参数 | 60.0 |
| `PORO_B` | $B$ | 生长速率参数 | 0.173 |

---

## 6. 效果对比

### 6.1 局部模型的问题

```
应变分布（局部模型）：

    晶粒1    晶粒2    晶粒3
    ε=0.30   ε=0.05   ε=0.08   ← 应变高度不均匀
    f=0.15   f=0.001  f=0.002  ← 孔隙率集中在晶粒1

    → 晶粒1过早失效，整体曲线不显示软化
    → 数值不收敛
```

### 6.2 非局部模型的改进

```
应变分布（非局部模型）：

    晶粒1    晶粒2    晶粒3
    ε=0.30   ε=0.05   ε=0.08
    ε_nl=0.14 (晶粒平均)

    → 形核使用 ε_nl=0.14，而非 ε=0.30
    → 孔隙率增长被"稀释"
    → 变形更均匀，数值稳定
    → 整体曲线显示合理的软化行为
```

---

## 7. 参考文献

1. Pijaudier-Cabot, G., & Bažant, Z. P. (1987). Nonlocal damage theory. *Journal of Engineering Mechanics*, 113(10), 1512-1533.

2. Peerlings, R. H. J., et al. (1996). Gradient enhanced damage for quasi-brittle materials. *International Journal for Numerical Methods in Engineering*, 39(19), 3391-3403.

3. Engelen, R. A. B., et al. (2003). Nonlocal implicit gradient-enhanced elasto-plasticity for the modelling of softening behaviour. *International Journal of Plasticity*, 19(4), 403-433.

4. Linse, T., et al. (2012). Simulation of crack propagation using a gradient-enriched ductile damage model based on dilatational strain. *Engineering Fracture Mechanics*, 95, 13-28.

---

## 8. 版本历史

| 版本 | 日期 | 说明 |
|------|------|------|
| v0.8-local-GTN | 2026-01-26 | 局部化GTN模型（修改前） |
| v0.9-nonlocal-GTN | 2026-01-26 | 非局部化GTN模型（当前版本） |

---

*文档版本: v1.0*
*最后更新: 2026-01-26*
*作者: 刘博韬 (SJTU)*
