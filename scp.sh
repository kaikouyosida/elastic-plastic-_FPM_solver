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
        echo "Operation 1: Copying all files to $HOSTNAME:~/elastic-plastic-_FPM_solver_5/"
        rsync -avzP "$target_dir" "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"
        ;;
    2)
        echo "Operation 2: Copying src to $HOSTNAME:~$target_dir"
        rsync -avzP "$target_dir/Data_Files_Input/"*.dat "$HOSTNAME:~/elastic-plastic-_FPM_solver_2/Data_Files_Input/"
        ;;
    3)
        # echo "Operation 3: Copying output"
        # rm -rf "$target_dir/Data_Files_Output"
        # mkdir -p "$target_dir/Data_Files_Output"
        # rsync -avzP "$HOSTNAME:~/elastic-plastic-_FPM_solver_2/Data_Files_Output/" "$target_dir/Data_Files_Output"
        # rm -rf "$target_dir/paraview"
        # mkdir -p "$target_dir/paraview"
        # rsync -avzP "$HOSTNAME:~/elastic-plastic-_FPM_solver_2/paraview/" "$target_dir/paraview"
        rm -rf "$target_dir/debug_for_residual"
        mkdir -p "$target_dir/debug_for_residual"
        rsync -avzP "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/debug_for_residual/" "$target_dir/debug_for_residual"
        ;;
    4)
        echo "Operation 4: Copying src to $HOSTNAME:~$target_dir sending any file"
        
        rsync -avzP "$target_dir/"Analysis_condition_setting.dat "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"
        rsync -avzP "$target_dir/"material_constant_setting.dat "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"

        rsync -avzP "$target_dir/Data_Files_Input/"Set_boundary_condition.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_2/Data_Files_Input/"
        rsync -avzP "$target_dir/Data_Files_Input/"Set_point_coordinate.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_2/Data_Files_Input/"

        rsync -avzP "$target_dir/"type.h "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"

        rsync -avzP "$target_dir/"b_matrix.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"
        rsync -avzP "$target_dir/"coefficient_matrix.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"
        rsync -avzP "$target_dir/"d_matrix.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"
        rsync -avzP "$target_dir/"external_force.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"
        rsync -avzP "$target_dir/"field.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"
        rsync -avzP "$target_dir/"fpm.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"
        rsync -avzP "$target_dir/"GetGaussPoints.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"
        rsync -avzP "$target_dir/"ImposeDirichretCondition.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"
        rsync -avzP "$target_dir/"internal_force.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"
        rsync -avzP "$target_dir/"external_force.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"
        rsync -avzP "$target_dir/"main.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"
        rsync -avzP "$target_dir/"matrix.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"
        rsync -avzP "$target_dir/"MKL_solver.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"
        rsync -avzP "$target_dir/"model.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"        
        rsync -avzP "$target_dir/"Output_data.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"
        rsync -avzP "$target_dir/"s_matrix.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"
        rsync -avzP "$target_dir/"scalar.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"        
        rsync -avzP "$target_dir/"ss_curve.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"
        rsync -avzP "$target_dir/"ss_curve.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"
        rsync -avzP "$target_dir/"ss_curve.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"
        rsync -avzP "$target_dir/"stress.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"
        rsync -avzP "$target_dir/"tensor.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"        
        rsync -avzP "$target_dir/"vector.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_5/"
        
        rsync -avzP "$target_dir/Data_Files_Output/"Post_process.c "$HOSTNAME:~/elastic-plastic-_FPM_solver_2/Data_Files_Output/"
        
        ;;
    *)
        echo "Invalid operation number. Please use 1, 2, or 3."
        exit 1
        ;;
esac

echo "Operation completed."