#!/bin/bash

# Convert script file to LF if necessary
[ "$(head -n 1 "$0" | grep -o $'\r')" ] && exec 2> /dev/null < <(tr -d '\r' < "$0") || true

# このスクリプトがあるディレクトリを取得
script_dir=$(dirname "$(realpath "$0")")

# ターゲットディレクトリをスクリプトのあるディレクトリに設定
target_dir="$script_dir"

# 引数チェック
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <hostname> <operation_number>"
    exit 1
fi

# 引数を変数に格納
HOSTNAME=$1
OPERATION_NUMBER=$2

# 操作番号に基づいてSCPでファイルを転送
case $OPERATION_NUMBER in
    1)
        echo "Operation 1: Copying all files to $HOSTNAME:~/"
        rsync -avzP "$target_dir" "$HOSTNAME:~/"
        ;;
    2)
        echo "Operation 2: Copying src to $HOSTNAME:~/$target_dir"
        rsync -avzP "$target_dir/src" "$HOSTNAME:~/$target_dir"
        rsync -avzP "$target_dir/analysis/input/"*.txt "$HOSTNAME:~/$target_dir/analysis/input"
        rsync -avzP "$target_dir/analysis/material/"*.txt "$HOSTNAME:~/$target_dir/analysis/material"
        rsync -avzP "$target_dir/analysis/material/"*.ini "$HOSTNAME:~/$target_dir/analysis/material"
        rsync -avzP "$target_dir/analysis/material/"*.sh "$HOSTNAME:~/$target_dir/analysis/material"
        rsync -avzP "$target_dir/Makefile" "$HOSTNAME:~/$target_dir"
        ;;
    3)
        echo "Operation 3: Copying output"
        rm -rf "$target_dir/analysis/output"
        mkdir -p "$target_dir/analysis/output/bin"
        rsync -avzP "$HOSTNAME:~/$target_dir/analysis/output" "$target_dir/analysis"
        ;;
    *)
        echo "Invalid operation number. Please use 1, 2, or 3."
        exit 1
        ;;
esac

echo "Operation completed."