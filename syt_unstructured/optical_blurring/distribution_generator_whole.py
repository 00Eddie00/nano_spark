import numpy as np
import os
from concentration_matrix_generator import *
from tool.shrink_file import *


# 荧光染料卷积后空间分布
def spatial_distribution(dirname, nums, kernel, xy_c_val, z_c_val):
    path = "../result/NANO/"
    fluo_file_name = f"{path}{dirname}/000{nums}.csv"
    fluo_concentration = np.loadtxt(fluo_file_name)
    processed_con_matrix = nano_process_concentration(fluo_concentration, xy_c_val)  # 处理得到荧光钙离子的浓度矩阵
    matrix = convolve3d(processed_con_matrix, kernel, xy_c_val, z_c_val)
    return matrix


# 荧光染料卷积后时间分布
def temporal_distribution(dirname, kernel, xy_c_val, z_c_val, position_list):
    print("temporal_distribution开始")
    nano_path = "../result/NANO/"
    open_path = "../result/OPEN/"
    # 初始化各点Ca浓度
    nano_filenames = os.listdir(f"{nano_path}{dirname}/")  # 纳米空间目前所有步数的浓度文件
    open_filenames = os.listdir(f"{open_path}{dirname}/")  # 纳米空间目前所有步数的浓度文件
    fluo_dir_list_len = len(nano_filenames)  # 当前生成文件数
    position_list_len = len(position_list)
    # 哪个文件，该文件取四个点 [[0, 0], [100, 0], [300, 0], [400, 0]]
    position_con = np.empty((position_list_len, fluo_dir_list_len))
    # 文件数，每隔10个（或100个）文件计算一次
    for fluo_file_index in range(fluo_dir_list_len):
        nano_file = nano_filenames[fluo_file_index]
        open_file = open_filenames[fluo_file_index]
        print(f"{nano_file}{open_file}开始")
        nano_c = np.loadtxt(f"{nano_path}{dirname}/{nano_file}")
        open_c = np.loadtxt(f"{open_path}{dirname}/{open_file}")
        processed_con_matrix = process_concentration(nano_c, open_c, xy_c_val)
        matrix = convolve3d(processed_con_matrix, kernel, xy_c_val, z_c_val)
        # 4
        for position_index in range(position_list_len):
            a = position_list[position_index][0]
            b = position_list[position_index][1]
            i = a + 500  # x:0+300/300+300
            j = b + 500  # y:0+300/0+300
            # 在纳米空间
            if np.sqrt(a ** 2 + b ** 2) <= 300:

                position_con[position_index, fluo_file_index] = matrix[i, j, 15]
            # 在开放空间
            else:
                position_con[position_index, fluo_file_index] = matrix[i, j, 500]

        print(f"{nano_file}{open_file}结束")
    return position_con


def optical_blurring(dirname, position_list):
    kernel = np.load("../optical_blurring/kernel_v1.npy", allow_pickle=True)
    result_set = "../result/NANO"
    # C_VAL = 0.0
    # 开放空间中没有
    if dirname == "CaG":
        C_VAL = 0.0  # 用于纳米钙火花
    elif dirname == "CaF":
        C_VAL = 4.081632653061224858e-03  # 用于钙火花
    xy_c_val = C_VAL
    z_c_val = C_VAL
    multiple = 100
    position_con = temporal_distribution(dirname, kernel, xy_c_val, 0, position_list)
    new_position_con = shrink_files(position_con, multiple)
    # 2
    for position_index in range(len(new_position_con)):
        np.savetxt(
            f"{result_set}/{dirname}_psf_({position_list[position_index][0]},{position_list[position_index][1]}).csv"
            , new_position_con[position_index])


def main():
    dirname = "CaF"
    my_position_list = [[0, 0], [100, 0], [300, 0], [400, 0]]
    optical_blurring(dirname, my_position_list)


if __name__ == "__main__":
    main()
