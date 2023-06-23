import numpy as np
import os
from concentration_matrix_generator import *


# 荧光染料卷积后空间分布
def spatial_distribution(dirname, nums, kernel, xy_c_val, z_c_val):
    path = "../result/NANO/"
    fluo_file_name = f"{path}{dirname}/000{nums}.csv"
    fluo_concentration = np.loadtxt(fluo_file_name)
    processed_con_matrix = process_concentration(fluo_concentration, xy_c_val)  # 处理得到荧光钙离子的浓度矩阵
    matrix = convolve3d(processed_con_matrix, kernel, xy_c_val, z_c_val)
    return matrix


# 荧光染料卷积后时间分布
def temporal_distribution(dirname, kernel, xy_c_val, z_c_val, position_list):
    print("temporal_distribution开始")
    path = "../result/NANO/"
    # 初始化各点Ca浓度
    filenames = os.listdir(f"{path}{dirname}/")  # 目前所有步数的浓度文件

    fluo_dir_list_len = len(filenames)  # 当前生成文件数
    # my_position_list = [[0, 0], [300, 0]]
    position_list_len = len(position_list)
    # 哪个文件，该文件取两个点 (0, 0) (300, 0)
    position_con = np.empty((fluo_dir_list_len, position_list_len))
    # 文件数，每隔十个文件计算一次
    for fluo_file_index in range(0, fluo_dir_list_len):
        filename = filenames[fluo_file_index]
        print(f"{filename}开始")
        nano_c = np.loadtxt(f"{path}{dirname}/{filename}")
        processed_con_matrix = process_concentration(nano_c, xy_c_val)
        matrix = convolve3d(processed_con_matrix, kernel, xy_c_val, z_c_val)
        # 2
        for position_index in range(position_list_len):
            i = position_list[position_index][0] + 300  # x:0+300/300+300
            j = position_list[position_index][1] + 300  # y:0+300/0+300
            position_con[fluo_file_index, position_index] = matrix[i, j, 15]
        print(f"{filename}结束")
    return position_con


def optical_blurring(dirname, position_list):
    kernel = np.load("../optical_blurring/kernel_v1.npy", allow_pickle=True)
    result_set = "../result/NANO"
    C_VAL = 0.0
    if dirname == "CaG":
        C_VAL = 0.0  # 用于纳米钙火花
    elif dirname == "CaF":
        C_VAL = 4.081632653061224858e-03  # 用于钙火花
    xy_c_val = C_VAL
    z_c_val = C_VAL
    position_con = temporal_distribution(dirname, kernel, xy_c_val, 0, position_list)
    # position_con.shape[1]=2
    for position_index in range(position_con.shape[1]):
        np.savetxt(
            f"{result_set}/{dirname}_psf_({position_list[position_index][0]},{position_list[position_index][1]}).csv"
            , position_con[:, position_index])


def main():
    dirname = "CaF"
    my_position_list = [[0, 0], [300, 0]]
    optical_blurring(dirname, my_position_list)


if __name__ == "__main__":
    main()
