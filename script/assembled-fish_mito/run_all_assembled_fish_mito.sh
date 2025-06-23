#!/bin/bash

# 获取当前脚本所在目录
assembled_fish_mito_path="$(dirname "$0")"

# 遍历目录中的所有 Python 文件并执行
for python_file in "$assembled_fish_mito_path"/*.py; do
    echo "Running $python_file..."
    python3 "$python_file"
done