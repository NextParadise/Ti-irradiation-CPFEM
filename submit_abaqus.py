#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Abaqus作业提交脚本
用于自动化提交CPFEM计算任务

作者: 刘博韬
日期: 2026-01-23
"""

import os
import sys
import subprocess
import argparse
import shutil
from datetime import datetime

# ============== 配置参数 ==============
DEFAULT_CONFIG = {
    'job_name': 'Ti13GN',           # 作业名称（不含扩展名）
    'inp_file': 'Ti13GN.inp',       # 输入文件
    'user_subroutine': 'TI_IHCPFEM.for',  # 用户子程序
    'cpus': 4,                      # CPU核心数
    'memory': '4gb',                # 内存限制
    'scratch': '',                  # 临时文件目录（空则使用默认）
    'abaqus_cmd': 'abaqus',         # Abaqus命令（可能是abq2020, abaqus等）
}

# 需要复制到工作目录的文件列表
REQUIRED_FILES = [
    'TI_IHCPFEM.for',
    'HCP_slip_twin.for',
    'DIS_IRRA_HARD.for',
    'POROSITY.for',
    'FIT.for',
    'GEULER.for',
    'ELEGRN.for',
]


def check_files_exist(work_dir, files):
    """检查必要文件是否存在"""
    missing = []
    for f in files:
        if not os.path.exists(os.path.join(work_dir, f)):
            missing.append(f)
    return missing


def create_work_directory(base_dir, job_name):
    """创建带时间戳的工作目录"""
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    work_dir = os.path.join(base_dir, f'{job_name}_{timestamp}')
    os.makedirs(work_dir, exist_ok=True)
    return work_dir


def copy_files_to_workdir(src_dir, work_dir, files):
    """复制文件到工作目录"""
    for f in files:
        src = os.path.join(src_dir, f)
        dst = os.path.join(work_dir, f)
        if os.path.exists(src):
            shutil.copy2(src, dst)
            print(f'  复制: {f}')
        else:
            print(f'  警告: 文件不存在 {f}')


def submit_job(config, work_dir=None, dry_run=False):
    """
    提交Abaqus作业

    参数:
        config: 配置字典
        work_dir: 工作目录（None则在当前目录运行）
        dry_run: 仅打印命令，不实际执行
    """
    # 构建命令
    cmd_parts = [
        config['abaqus_cmd'],
        'job=' + config['job_name'],
        'input=' + config['inp_file'],
        'user=' + config['user_subroutine'],
        'cpus=' + str(config['cpus']),
        'memory=' + config['memory'],
        'interactive',  # 交互模式，等待完成
    ]

    if config['scratch']:
        cmd_parts.append('scratch=' + config['scratch'])

    cmd = ' '.join(cmd_parts)

    print('=' * 60)
    print('Abaqus作业提交')
    print('=' * 60)
    print(f'作业名称: {config["job_name"]}')
    print(f'输入文件: {config["inp_file"]}')
    print(f'用户子程序: {config["user_subroutine"]}')
    print(f'CPU核心数: {config["cpus"]}')
    print(f'内存: {config["memory"]}')
    if work_dir:
        print(f'工作目录: {work_dir}')
    print('-' * 60)
    print(f'执行命令: {cmd}')
    print('=' * 60)

    if dry_run:
        print('[DRY RUN] 命令未执行')
        return 0

    # 切换到工作目录并执行
    original_dir = os.getcwd()
    try:
        if work_dir:
            os.chdir(work_dir)

        print(f'\n开始时间: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
        print('正在运行Abaqus...\n')

        result = subprocess.run(cmd, shell=True)

        print(f'\n结束时间: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')

        if result.returncode == 0:
            print('\n作业完成!')
        else:
            print(f'\n作业失败，返回码: {result.returncode}')

        return result.returncode

    finally:
        os.chdir(original_dir)


def submit_batch_jobs(param_sets, base_config, src_dir, output_base_dir):
    """
    批量提交参数扫描作业

    参数:
        param_sets: 参数集列表，每个元素是(名称, 参数字典)
        base_config: 基础配置
        src_dir: 源文件目录
        output_base_dir: 输出基础目录
    """
    results = []

    for name, params in param_sets:
        print(f'\n{"="*60}')
        print(f'参数集: {name}')
        print(f'{"="*60}')

        # 创建工作目录
        work_dir = create_work_directory(output_base_dir, name)

        # 复制文件
        files_to_copy = REQUIRED_FILES + [base_config['inp_file']]
        copy_files_to_workdir(src_dir, work_dir, files_to_copy)

        # 修改FIT.for中的参数（如果需要）
        if params:
            modify_fit_parameters(work_dir, params)

        # 提交作业
        config = base_config.copy()
        config['job_name'] = name
        returncode = submit_job(config, work_dir)

        results.append((name, work_dir, returncode))

    # 打印汇总
    print('\n' + '=' * 60)
    print('批量作业汇总')
    print('=' * 60)
    for name, work_dir, returncode in results:
        status = '成功' if returncode == 0 else f'失败({returncode})'
        print(f'{name}: {status}')
        print(f'  目录: {work_dir}')

    return results


def modify_fit_parameters(work_dir, params):
    """
    修改FIT.for中的参数

    参数:
        work_dir: 工作目录
        params: 参数字典，如 {'PMUL': 0.5, 'PDYN': 100.0}
    """
    fit_file = os.path.join(work_dir, 'FIT.for')

    if not os.path.exists(fit_file):
        print(f'  警告: FIT.for不存在，跳过参数修改')
        return

    with open(fit_file, 'r') as f:
        content = f.read()

    for param_name, param_value in params.items():
        # 匹配形如 "PARAM=value" 或 "PARAM = value" 的模式
        import re
        pattern = rf'({param_name}\s*=\s*)[\d.eEdD+-]+'
        replacement = rf'\g<1>{param_value}'
        content, count = re.subn(pattern, replacement, content, flags=re.IGNORECASE)
        if count > 0:
            print(f'  修改参数: {param_name} = {param_value}')
        else:
            print(f'  警告: 未找到参数 {param_name}')

    with open(fit_file, 'w') as f:
        f.write(content)


def main():
    parser = argparse.ArgumentParser(
        description='Abaqus CPFEM作业提交脚本',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例:
  # 使用默认配置提交
  python submit_abaqus.py

  # 指定作业名和CPU数
  python submit_abaqus.py -j MyJob -c 8

  # 仅打印命令，不执行
  python submit_abaqus.py --dry-run

  # 创建新工作目录并提交
  python submit_abaqus.py --new-dir
        '''
    )

    parser.add_argument('-j', '--job', default=DEFAULT_CONFIG['job_name'],
                        help=f'作业名称 (默认: {DEFAULT_CONFIG["job_name"]})')
    parser.add_argument('-i', '--input', default=DEFAULT_CONFIG['inp_file'],
                        help=f'输入文件 (默认: {DEFAULT_CONFIG["inp_file"]})')
    parser.add_argument('-u', '--user', default=DEFAULT_CONFIG['user_subroutine'],
                        help=f'用户子程序 (默认: {DEFAULT_CONFIG["user_subroutine"]})')
    parser.add_argument('-c', '--cpus', type=int, default=DEFAULT_CONFIG['cpus'],
                        help=f'CPU核心数 (默认: {DEFAULT_CONFIG["cpus"]})')
    parser.add_argument('-m', '--memory', default=DEFAULT_CONFIG['memory'],
                        help=f'内存限制 (默认: {DEFAULT_CONFIG["memory"]})')
    parser.add_argument('--abaqus', default=DEFAULT_CONFIG['abaqus_cmd'],
                        help=f'Abaqus命令 (默认: {DEFAULT_CONFIG["abaqus_cmd"]})')
    parser.add_argument('--dry-run', action='store_true',
                        help='仅打印命令，不实际执行')
    parser.add_argument('--new-dir', action='store_true',
                        help='创建新的工作目录')
    parser.add_argument('--check', action='store_true',
                        help='仅检查文件是否存在')

    args = parser.parse_args()

    # 构建配置
    config = {
        'job_name': args.job,
        'inp_file': args.input,
        'user_subroutine': args.user,
        'cpus': args.cpus,
        'memory': args.memory,
        'scratch': '',
        'abaqus_cmd': args.abaqus,
    }

    # 获取脚本所在目录
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # 检查文件
    files_to_check = REQUIRED_FILES + [config['inp_file']]
    missing = check_files_exist(script_dir, files_to_check)

    if missing:
        print('错误: 以下文件缺失:')
        for f in missing:
            print(f'  - {f}')
        sys.exit(1)

    if args.check:
        print('所有必要文件存在')
        sys.exit(0)

    # 确定工作目录
    work_dir = None
    if args.new_dir:
        work_dir = create_work_directory(script_dir, config['job_name'])
        print(f'创建工作目录: {work_dir}')
        files_to_copy = REQUIRED_FILES + [config['inp_file']]
        copy_files_to_workdir(script_dir, work_dir, files_to_copy)

    # 提交作业
    returncode = submit_job(config, work_dir, args.dry_run)
    sys.exit(returncode)


if __name__ == '__main__':
    main()
