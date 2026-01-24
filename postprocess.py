#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
CPFEM后处理脚本
从Abaqus ODB文件提取并可视化各种物理量

作者: 刘博韬
日期: 2026-01-23

功能:
    1. 应力-应变曲线
    2. 孪晶体积分数演化
    3. 孔隙率演化
    4. 位错密度演化
    5. 各滑移系活动分析
    6. 参数敏感性对比

使用方法:
    # 在Abaqus Python环境中运行
    abaqus python postprocess.py --odb Ti13GN.odb

    # 或者先提取数据，再用标准Python绘图
    abaqus python postprocess.py --odb Ti13GN.odb --extract-only
    python postprocess.py --plot-only --data results.npz
"""

from __future__ import print_function
import os
import sys
import argparse
import numpy as np

# 尝试导入Abaqus模块（仅在Abaqus Python环境中可用）
try:
    from odbAccess import openOdb
    from abaqusConstants import *
    HAS_ABAQUS = True
except ImportError:
    HAS_ABAQUS = False

# 尝试导入matplotlib（用于绘图）
try:
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.rcParams['font.family'] = 'DejaVu Sans'
    matplotlib.rcParams['axes.unicode_minus'] = False
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


# ============== 状态变量索引定义 ==============
# 根据TI_IHCPFEM.for中的STATEV结构定义
ND = 42  # 滑移系+孪晶系总数
NSLPTL = 30  # 滑移系数量
NTWTL = 12  # 孪晶系数量

# 状态变量索引（Python从0开始，Fortran从1开始，需要-1）
class SDVIndex:
    """状态变量索引类"""
    # 基本变量
    STRENGTH = (0, ND)                    # 0*ND ~ 1*ND: 滑移系强度
    GAMMA = (ND, 2*ND)                    # 1*ND ~ 2*ND: 累积滑移量
    TAU = (2*ND, 3*ND)                    # 2*ND ~ 3*ND: 分切应力

    # 几何信息
    SLPNOR = (3*ND, 6*ND)                 # 滑移面法向
    SLPDIR = (6*ND, 9*ND)                 # 滑移方向

    # 演化变量
    GAMMA_ABS = (9*ND, 10*ND)             # 累积滑移量绝对值
    ANG = (10*ND, 11*ND)                  # 角度
    SMD = (11*ND, 12*ND)                  # 施密特因子

    # 位错密度
    NMD = (14*ND, 15*ND)                  # 可动位错密度

    # 辐照缺陷
    NVOID = (18*ND, 19*ND)                # 空隙密度
    NLOOP = (19*ND, 20*ND)                # 位错环密度
    NPRE = (20*ND, 21*ND)                 # 沉淀相密度
    EMFP = (21*ND, 22*ND)                 # 平均自由程

    # 全局统计量（从末尾计算）
    @staticmethod
    def get_global_indices(nstatv):
        """获取全局统计量的索引"""
        return {
            'POROSITY': nstatv - 31,       # 孔隙率
            'EQVPL': nstatv - 30,          # 等效塑性应变
            'SLIP_SUM': nstatv - 29,       # 滑移系切应变总和
            'TWIN_SUM': nstatv - 28,       # 孪晶体积分数总和
        }


def extract_data_from_odb(odb_path, output_file=None):
    """
    从ODB文件提取数据

    参数:
        odb_path: ODB文件路径
        output_file: 输出文件路径（.npz格式）

    返回:
        包含所有提取数据的字典
    """
    if not HAS_ABAQUS:
        raise RuntimeError("此函数需要在Abaqus Python环境中运行")

    print(f"打开ODB文件: {odb_path}")
    odb = openOdb(path=odb_path, readOnly=True)

    # 获取步骤和帧信息
    step_name = list(odb.steps.keys())[-1]  # 使用最后一个步骤
    step = odb.steps[step_name]
    frames = step.frames

    print(f"步骤: {step_name}")
    print(f"帧数: {len(frames)}")

    # 获取实例
    instance_name = list(odb.rootAssembly.instances.keys())[0]
    instance = odb.rootAssembly.instances[instance_name]

    print(f"实例: {instance_name}")
    print(f"单元数: {len(instance.elements)}")

    # 初始化数据存储
    n_frames = len(frames)
    n_elements = len(instance.elements)

    # 获取单元体积（用于体积平均）
    element_volumes = np.ones(n_elements)  # 简化处理，假设等体积

    # 数据数组
    time_data = np.zeros(n_frames)
    stress_data = np.zeros((n_frames, 6))  # 体积平均应力
    strain_data = np.zeros((n_frames, 6))  # 体积平均应变

    # 状态变量数据
    porosity_data = np.zeros(n_frames)
    twin_sum_data = np.zeros(n_frames)
    slip_sum_data = np.zeros(n_frames)
    eqvpl_data = np.zeros(n_frames)

    # 位错密度数据
    dislocation_data = np.zeros(n_frames)

    # 各滑移系数据（体积平均）
    gamma_data = np.zeros((n_frames, ND))
    strength_data = np.zeros((n_frames, ND))

    print("\n提取数据中...")

    for i_frame, frame in enumerate(frames):
        time_data[i_frame] = frame.frameValue

        if i_frame % 10 == 0:
            print(f"  帧 {i_frame}/{n_frames}, 时间 = {frame.frameValue:.4f}")

        # 提取应力
        if 'S' in frame.fieldOutputs:
            stress_field = frame.fieldOutputs['S']
            stress_values = stress_field.getSubset(region=instance)
            stress_sum = np.zeros(6)
            vol_sum = 0.0
            for val in stress_values.values:
                vol = element_volumes[0]  # 简化
                stress_sum += np.array(val.data) * vol
                vol_sum += vol
            stress_data[i_frame] = stress_sum / vol_sum

        # 提取应变
        if 'E' in frame.fieldOutputs:
            strain_field = frame.fieldOutputs['E']
            strain_values = strain_field.getSubset(region=instance)
            strain_sum = np.zeros(6)
            vol_sum = 0.0
            for val in strain_values.values:
                vol = element_volumes[0]
                strain_sum += np.array(val.data) * vol
                vol_sum += vol
            strain_data[i_frame] = strain_sum / vol_sum

        # 提取状态变量
        sdv_keys = [k for k in frame.fieldOutputs.keys() if k.startswith('SDV')]
        if sdv_keys:
            # 获取NSTATV
            nstatv = len(sdv_keys)
            global_idx = SDVIndex.get_global_indices(nstatv)

            # 孔隙率
            sdv_key = f'SDV{global_idx["POROSITY"] + 1}'
            if sdv_key in frame.fieldOutputs:
                field = frame.fieldOutputs[sdv_key]
                values = field.getSubset(region=instance)
                porosity_data[i_frame] = np.mean([v.data for v in values.values])

            # 孪晶体积分数
            sdv_key = f'SDV{global_idx["TWIN_SUM"] + 1}'
            if sdv_key in frame.fieldOutputs:
                field = frame.fieldOutputs[sdv_key]
                values = field.getSubset(region=instance)
                twin_sum_data[i_frame] = np.mean([v.data for v in values.values])

            # 滑移总量
            sdv_key = f'SDV{global_idx["SLIP_SUM"] + 1}'
            if sdv_key in frame.fieldOutputs:
                field = frame.fieldOutputs[sdv_key]
                values = field.getSubset(region=instance)
                slip_sum_data[i_frame] = np.mean([v.data for v in values.values])

            # 等效塑性应变
            sdv_key = f'SDV{global_idx["EQVPL"] + 1}'
            if sdv_key in frame.fieldOutputs:
                field = frame.fieldOutputs[sdv_key]
                values = field.getSubset(region=instance)
                eqvpl_data[i_frame] = np.mean([v.data for v in values.values])

            # 位错密度（取第一个滑移系作为代表）
            sdv_key = f'SDV{SDVIndex.NMD[0] + 1}'
            if sdv_key in frame.fieldOutputs:
                field = frame.fieldOutputs[sdv_key]
                values = field.getSubset(region=instance)
                dislocation_data[i_frame] = np.mean([v.data for v in values.values])

            # 各滑移系累积滑移量
            for j in range(ND):
                sdv_key = f'SDV{SDVIndex.GAMMA_ABS[0] + j + 1}'
                if sdv_key in frame.fieldOutputs:
                    field = frame.fieldOutputs[sdv_key]
                    values = field.getSubset(region=instance)
                    gamma_data[i_frame, j] = np.mean([v.data for v in values.values])

    odb.close()
    print("ODB文件已关闭")

    # 组织数据
    data = {
        'time': time_data,
        'stress': stress_data,
        'strain': strain_data,
        'porosity': porosity_data,
        'twin_sum': twin_sum_data,
        'slip_sum': slip_sum_data,
        'eqvpl': eqvpl_data,
        'dislocation': dislocation_data,
        'gamma': gamma_data,
        'strength': strength_data,
    }

    # 保存数据
    if output_file:
        np.savez(output_file, **data)
        print(f"数据已保存到: {output_file}")

    return data


def load_data(data_file):
    """加载保存的数据"""
    data = np.load(data_file)
    return {key: data[key] for key in data.files}


def plot_stress_strain(data, direction=2, save_path=None, exp_data=None):
    """
    绘制应力-应变曲线

    参数:
        data: 数据字典
        direction: 方向索引 (0=X, 1=Y, 2=Z)
        save_path: 保存路径
        exp_data: 实验数据 (strain, stress) 元组
    """
    if not HAS_MATPLOTLIB:
        print("警告: matplotlib不可用，跳过绘图")
        return

    fig, ax = plt.subplots(figsize=(8, 6))

    strain = data['strain'][:, direction]
    stress = data['stress'][:, direction]

    # 转换单位（假设应力为Pa，转换为MPa）
    stress_mpa = stress / 1e6

    ax.plot(strain, stress_mpa, 'b-', linewidth=2, label='Simulation')

    if exp_data is not None:
        exp_strain, exp_stress = exp_data
        ax.plot(exp_strain, exp_stress, 'ro', markersize=4, label='Experiment')

    ax.set_xlabel('Strain', fontsize=12)
    ax.set_ylabel('Stress (MPa)', fontsize=12)
    ax.set_title('Stress-Strain Curve', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"图片已保存: {save_path}")

    plt.show()


def plot_twin_evolution(data, save_path=None):
    """绘制孪晶体积分数演化曲线"""
    if not HAS_MATPLOTLIB:
        return

    fig, ax = plt.subplots(figsize=(8, 6))

    strain = data['strain'][:, 2]  # Z方向应变
    twin_sum = data['twin_sum']

    ax.plot(strain, twin_sum, 'g-', linewidth=2)

    ax.set_xlabel('Strain', fontsize=12)
    ax.set_ylabel('Twin Volume Fraction', fontsize=12)
    ax.set_title('Twin Evolution', fontsize=14)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"图片已保存: {save_path}")

    plt.show()


def plot_porosity_evolution(data, save_path=None):
    """绘制孔隙率演化曲线"""
    if not HAS_MATPLOTLIB:
        return

    fig, ax = plt.subplots(figsize=(8, 6))

    strain = data['strain'][:, 2]
    porosity = data['porosity']

    ax.plot(strain, porosity, 'r-', linewidth=2)

    ax.set_xlabel('Strain', fontsize=12)
    ax.set_ylabel('Porosity', fontsize=12)
    ax.set_title('Porosity Evolution', fontsize=14)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"图片已保存: {save_path}")

    plt.show()


def plot_dislocation_evolution(data, save_path=None):
    """绘制位错密度演化曲线"""
    if not HAS_MATPLOTLIB:
        return

    fig, ax = plt.subplots(figsize=(8, 6))

    strain = data['strain'][:, 2]
    dislocation = data['dislocation']

    ax.plot(strain, dislocation, 'm-', linewidth=2)

    ax.set_xlabel('Strain', fontsize=12)
    ax.set_ylabel('Dislocation Density (mm$^{-2}$)', fontsize=12)
    ax.set_title('Dislocation Density Evolution', fontsize=14)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"图片已保存: {save_path}")

    plt.show()


def plot_slip_activity(data, save_path=None):
    """绘制各滑移系活动分析"""
    if not HAS_MATPLOTLIB:
        return

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    strain = data['strain'][:, 2]
    gamma = data['gamma']

    # 滑移系分组
    slip_groups = {
        'Basal (1-3)': (0, 3),
        'Prism (4-6)': (3, 6),
        'Pyram (7-12)': (6, 12),
        'Pyram I (13-24)': (12, 24),
        'Pyram II (25-30)': (24, 30),
        'Twin (31-42)': (30, 42),
    }

    for ax, (name, (start, end)) in zip(axes.flat, slip_groups.items()):
        group_gamma = gamma[:, start:end]
        total_gamma = np.sum(group_gamma, axis=1)

        ax.plot(strain, total_gamma, linewidth=2)
        ax.set_xlabel('Strain', fontsize=10)
        ax.set_ylabel('Cumulative Shear', fontsize=10)
        ax.set_title(name, fontsize=12)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"图片已保存: {save_path}")

    plt.show()


def plot_comprehensive(data, save_path=None, exp_data=None):
    """绘制综合分析图"""
    if not HAS_MATPLOTLIB:
        return

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    strain = data['strain'][:, 2]

    # 1. 应力-应变曲线
    ax = axes[0, 0]
    stress_mpa = data['stress'][:, 2] / 1e6
    ax.plot(strain, stress_mpa, 'b-', linewidth=2, label='Simulation')
    if exp_data is not None:
        ax.plot(exp_data[0], exp_data[1], 'ro', markersize=4, label='Experiment')
        ax.legend()
    ax.set_xlabel('Strain')
    ax.set_ylabel('Stress (MPa)')
    ax.set_title('(a) Stress-Strain Curve')
    ax.grid(True, alpha=0.3)

    # 2. 孪晶体积分数
    ax = axes[0, 1]
    ax.plot(strain, data['twin_sum'], 'g-', linewidth=2)
    ax.set_xlabel('Strain')
    ax.set_ylabel('Twin Volume Fraction')
    ax.set_title('(b) Twin Evolution')
    ax.grid(True, alpha=0.3)

    # 3. 孔隙率
    ax = axes[1, 0]
    ax.plot(strain, data['porosity'], 'r-', linewidth=2)
    ax.set_xlabel('Strain')
    ax.set_ylabel('Porosity')
    ax.set_title('(c) Porosity Evolution')
    ax.grid(True, alpha=0.3)

    # 4. 位错密度
    ax = axes[1, 1]
    ax.plot(strain, data['dislocation'], 'm-', linewidth=2)
    ax.set_xlabel('Strain')
    ax.set_ylabel('Dislocation Density (mm$^{-2}$)')
    ax.set_title('(d) Dislocation Density Evolution')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"图片已保存: {save_path}")

    plt.show()


def compare_parameter_sensitivity(data_list, param_name, param_values, save_path=None):
    """
    参数敏感性对比分析

    参数:
        data_list: 数据字典列表
        param_name: 参数名称
        param_values: 参数值列表
        save_path: 保存路径
    """
    if not HAS_MATPLOTLIB:
        return

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    colors = plt.cm.viridis(np.linspace(0, 1, len(data_list)))

    for i, (data, param_val) in enumerate(zip(data_list, param_values)):
        strain = data['strain'][:, 2]
        label = f'{param_name}={param_val}'

        # 应力-应变
        axes[0, 0].plot(strain, data['stress'][:, 2]/1e6, color=colors[i],
                        linewidth=2, label=label)

        # 孪晶
        axes[0, 1].plot(strain, data['twin_sum'], color=colors[i],
                        linewidth=2, label=label)

        # 孔隙率
        axes[1, 0].plot(strain, data['porosity'], color=colors[i],
                        linewidth=2, label=label)

        # 位错密度
        axes[1, 1].plot(strain, data['dislocation'], color=colors[i],
                        linewidth=2, label=label)

    titles = ['Stress-Strain', 'Twin Evolution', 'Porosity', 'Dislocation Density']
    ylabels = ['Stress (MPa)', 'Twin Volume Fraction', 'Porosity', 'Dislocation Density']

    for ax, title, ylabel in zip(axes.flat, titles, ylabels):
        ax.set_xlabel('Strain')
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.suptitle(f'Parameter Sensitivity: {param_name}', fontsize=14)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"图片已保存: {save_path}")

    plt.show()


def export_to_csv(data, output_dir):
    """导出数据到CSV文件"""
    os.makedirs(output_dir, exist_ok=True)

    # 应力-应变数据
    ss_data = np.column_stack([
        data['strain'][:, 2],
        data['stress'][:, 2] / 1e6,
    ])
    np.savetxt(os.path.join(output_dir, 'stress_strain.csv'), ss_data,
               delimiter=',', header='Strain,Stress_MPa', comments='')

    # 孪晶数据
    twin_data = np.column_stack([
        data['strain'][:, 2],
        data['twin_sum'],
    ])
    np.savetxt(os.path.join(output_dir, 'twin_evolution.csv'), twin_data,
               delimiter=',', header='Strain,Twin_Volume_Fraction', comments='')

    # 孔隙率数据
    poro_data = np.column_stack([
        data['strain'][:, 2],
        data['porosity'],
    ])
    np.savetxt(os.path.join(output_dir, 'porosity_evolution.csv'), poro_data,
               delimiter=',', header='Strain,Porosity', comments='')

    # 位错密度数据
    disl_data = np.column_stack([
        data['strain'][:, 2],
        data['dislocation'],
    ])
    np.savetxt(os.path.join(output_dir, 'dislocation_evolution.csv'), disl_data,
               delimiter=',', header='Strain,Dislocation_Density', comments='')

    print(f"CSV文件已导出到: {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description='CPFEM后处理脚本',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument('--odb', type=str, help='ODB文件路径')
    parser.add_argument('--data', type=str, help='已保存的数据文件(.npz)')
    parser.add_argument('--output', type=str, default='results',
                        help='输出目录 (默认: results)')
    parser.add_argument('--extract-only', action='store_true',
                        help='仅提取数据，不绘图')
    parser.add_argument('--plot-only', action='store_true',
                        help='仅绘图（需要--data参数）')
    parser.add_argument('--exp-data', type=str,
                        help='实验数据CSV文件（两列：应变，应力MPa）')

    args = parser.parse_args()

    # 创建输出目录
    os.makedirs(args.output, exist_ok=True)

    # 加载实验数据
    exp_data = None
    if args.exp_data and os.path.exists(args.exp_data):
        exp = np.loadtxt(args.exp_data, delimiter=',', skiprows=1)
        exp_data = (exp[:, 0], exp[:, 1])
        print(f"已加载实验数据: {args.exp_data}")

    # 提取或加载数据
    if args.plot_only:
        if not args.data:
            print("错误: --plot-only需要--data参数")
            sys.exit(1)
        data = load_data(args.data)
    elif args.odb:
        data_file = os.path.join(args.output, 'extracted_data.npz')
        data = extract_data_from_odb(args.odb, data_file)
    else:
        print("错误: 需要--odb或--data参数")
        sys.exit(1)

    if args.extract_only:
        print("数据提取完成")
        return

    # 绘图
    print("\n生成图表...")

    plot_stress_strain(data,
                       save_path=os.path.join(args.output, 'stress_strain.png'),
                       exp_data=exp_data)

    plot_twin_evolution(data,
                        save_path=os.path.join(args.output, 'twin_evolution.png'))

    plot_porosity_evolution(data,
                            save_path=os.path.join(args.output, 'porosity_evolution.png'))

    plot_dislocation_evolution(data,
                               save_path=os.path.join(args.output, 'dislocation_evolution.png'))

    plot_slip_activity(data,
                       save_path=os.path.join(args.output, 'slip_activity.png'))

    plot_comprehensive(data,
                       save_path=os.path.join(args.output, 'comprehensive.png'),
                       exp_data=exp_data)

    # 导出CSV
    export_to_csv(data, args.output)

    print("\n后处理完成!")


if __name__ == '__main__':
    main()
