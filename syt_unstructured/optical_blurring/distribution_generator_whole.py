import numpy as np
import os
from concentration_matrix_generator import *
from tool.shrink_file import *
import time


# 荧光染料卷积后空间分布
def spatial_distribution(dirname, nums, kernel, xy_c_val, z_c_val):
    path = "../result/NANO/"
    fluo_file_name = f"{path}{dirname}/000{nums}.csv"
    fluo_concentration = np.loadtxt(fluo_file_name)
    processed_con_matrix = nano_process_concentration(fluo_concentration, xy_c_val)  # 处理得到荧光钙离子的浓度矩阵
    matrix = convolve3d(processed_con_matrix, kernel, xy_c_val, z_c_val)
    return matrix


# 荧光染料卷积后时间分布
# 点太多数据量太大，不好用
def temporal_distribution(dirname, kernel, xy_c_val, z_c_val, position_list, multiple):
    print("temporal_distribution开始")
    nano_path = "../result/NANO/"
    open_path = "../result/OPEN/"
    # 初始化各点Ca浓度
    nano_filenames = os.listdir(f"{nano_path}{dirname}/")  # 纳米空间目前所有步数的浓度文件
    open_filenames = os.listdir(f"{open_path}{dirname}/")  # 纳米空间目前所有步数的浓度文件
    fluo_dir_list_len = len(nano_filenames)  # 当前生成文件数
    position_list_len = len(position_list)
    # 哪个文件，该文件取四个点 [[0, 0], [100, 0], [300, 0], [400, 0]]
    position_con = np.empty((position_list_len, fluo_dir_list_len // multiple + 1))
    # 文件数，每隔10个（或100个）文件计算一次
    for fluo_file_index in range(0, fluo_dir_list_len, multiple):
        nano_file = nano_filenames[fluo_file_index]
        open_file = open_filenames[fluo_file_index]
        print(f"{nano_file}开始")
        nano_c = np.loadtxt(f"{nano_path}{dirname}/{nano_file}")
        open_c = np.loadtxt(f"{open_path}{dirname}/{open_file}")
        nano_processed = nano_process_concentration(nano_c, xy_c_val)
        # 4
        for position_index in range(position_list_len):
            processed_con_matrix = process_concentration(nano_processed, open_c, xy_c_val,
                                                         position_list[position_index])
            matrix = convolve3d(processed_con_matrix, kernel, xy_c_val, z_c_val)
            position_con[position_index, fluo_file_index] = matrix[500, 500, 500]
        print(f"{nano_file}{open_file}结束")
    return position_con


# 荧光染料卷积后时间分布
# 点太多数据量太大，不好用
def temporal_distribution_v2(dirname, kernel, xy_c_val, z_c_val, position_list, multiple):
    print("temporal_distribution开始")
    nano_path = "../result/NANO/"
    open_path = "../result/OPEN/"
    # 初始化各点Ca浓度
    nano_filenames = os.listdir(f"{nano_path}{dirname}/")  # 纳米空间目前所有步数的浓度文件
    open_filenames = os.listdir(f"{open_path}{dirname}/")  # 纳米空间目前所有步数的浓度文件
    fluo_dir_list_len = len(nano_filenames)  # 当前生成文件数
    position_list_len = len(position_list)  # 需要检测的点个数
    # 哪个文件，该文件取四个点 [[0, 0], [100, 0], [300, 0], [400, 0]]
    position_con = np.empty((position_list_len, fluo_dir_list_len // multiple + 1))
    max_relations = position_list[position_list_len - 1]
    grids_rz = np.load(f"grids_rz_v2_({max_relations[0]},{max_relations[1]}).npy")  # 400,0
    points_rz = np.load(f"points_rz_v2_({max_relations[0]},{max_relations[1]}).npy")

    # 文件数，每隔10个（或100个）文件计算一次
    for fluo_file_index in range(0, fluo_dir_list_len, multiple):
        nano_file = nano_filenames[fluo_file_index]
        open_file = open_filenames[fluo_file_index]
        print(f"{nano_file}开始")
        nano_c = np.loadtxt(f"{nano_path}{dirname}/{nano_file}")
        open_c = np.loadtxt(f"{open_path}{dirname}/{open_file}")
        nano_processed = nano_process_concentration(nano_c, xy_c_val)
        open_processed_whole = open_process_concentration_v2(nano_processed, open_c, xy_c_val,
                                                             max_relations, grids_rz, points_rz)  # 901 501 501
        # 4
        for position_index in range(position_list_len):
            position_list_i = position_list[position_index]
            processed_con_matrix = process_concentration_v2(nano_processed, open_processed_whole, xy_c_val,
                                                            position_list_i, points_rz)
            matrix = convolve3d(processed_con_matrix, kernel, xy_c_val, z_c_val)
            position_con[position_index, fluo_file_index] = matrix[500, 500, 500]
        print(f"{nano_file}结束")
    return position_con


# 比V4慢可能也没有V4准确
def temporal_distribution_v3(dirname, kernel, xy_c_val, z_c_val, position_list, multiple):
    print("temporal_distribution开始")
    nano_path = "../result/NANO/"
    open_path = "../result/OPEN/"
    # 初始化各点Ca浓度
    nano_filenames = os.listdir(f"{nano_path}{dirname}/")  # 纳米空间目前所有步数的浓度文件
    open_filenames = os.listdir(f"{open_path}{dirname}/")  # 纳米空间目前所有步数的浓度文件
    fluo_dir_list_len = len(nano_filenames)  # 当前生成文件数
    position_list_len = len(position_list)  # 需要检测的点个数
    # 哪个文件，该文件取四个点 [[0, 0], [100, 0], [300, 0], [400, 0]]
    position_con = np.empty((position_list_len, fluo_dir_list_len // multiple + 1))
    max_relations = position_list[position_list_len - 1]
    grids_rz = np.load(f"grids_rz_v2_({max_relations[0]},{max_relations[1]}).npy")  # 400,0
    # points_rz = np.load(f"points_rz_v2_({max_relations[0]},{max_relations[1]}).npy")

    # 文件数，每隔10个（或100个）文件计算一次
    for fluo_file_index in range(0, fluo_dir_list_len, multiple):
        nano_file = nano_filenames[fluo_file_index]
        open_file = open_filenames[fluo_file_index]
        print(f"{nano_file}开始")
        nano_c = np.loadtxt(f"{nano_path}{dirname}/{nano_file}")
        open_c = np.loadtxt(f"{open_path}{dirname}/{open_file}")
        nano_processed = nano_process_concentration(nano_c, xy_c_val)

        # 4
        for position_index in range(position_list_len):
            position_list_i = position_list[position_index]
            processed_con_matrix, mid = process_concentration_v3(nano_processed, open_c, xy_c_val,
                                                                 position_list_i, grids_rz)
            matrix = convolve3d(processed_con_matrix, kernel, xy_c_val, z_c_val)
            position_con[position_index, fluo_file_index // multiple] = matrix[mid - 1, mid - 1, mid - 1]
        print(f"{nano_file}结束")
    return position_con


# 一般用这个，最准确但是很慢60000步要跑6、7天
def temporal_distribution_v4(dirname, kernel, xy_c_val, z_c_val, position_list, multiple):
    print("temporal_distribution开始")
    nano_path = "../result/NANO/"
    open_path = "../result/OPEN/"
    nano_path = "../result/NANO/"
    open_path = "../result/OPEN/"
    # 初始化各点Ca浓度
    nano_filenames = os.listdir(f"{nano_path}{dirname}/")  # 纳米空间目前所有步数的浓度文件
    open_filenames = os.listdir(f"{open_path}{dirname}/")  # 纳米空间目前所有步数的浓度文件
    fluo_dir_list_len = len(nano_filenames)  # 当前生成文件数
    position_list_len = len(position_list)  # 需要检测的点个数
    # 哪个文件，该文件取四个点 [[0, 0], [100, 0], [300, 0], [400, 0]]
    position_con = np.empty((position_list_len, fluo_dir_list_len // multiple + 1))
    max_relations = position_list[position_list_len - 1]
    grids_rz = np.load(f"grids_rz_v2_({max_relations[0]},{max_relations[1]}).npy")  # 400,0
    xy_index = np.load(f"xy_list/xy_index.npy")
    new_open_r = np.flipud(point_scatter(295, 0, 5, k=2, positive=False))
    radius_list = np.concatenate((new_open_r, (300,)))
    grid_coordinates = np.loadtxt("../config/open/open_grid_coordinates.csv", delimiter=",")
    # points_rz = np.load(f"points_rz_v2_({max_relations[0]},{max_relations[1]}).npy")
    # 文件数，每隔10个（或100个）文件计算一次
    for fluo_file_index in range(0, fluo_dir_list_len, multiple):
        nano_file = nano_filenames[fluo_file_index]
        open_file = open_filenames[fluo_file_index]
        print(f"{nano_file}开始")
        nano_c = np.loadtxt(f"{nano_path}{dirname}/{nano_file}")
        open_c = np.loadtxt(f"{open_path}{dirname}/{open_file}")
        nano_processed = nano_process_concentration(nano_c, xy_c_val)
        # 4
        for position_index in range(position_list_len):
            position_list_i = position_list[position_index]
            processed_con_matrix, mid = process_concentration_v4(nano_processed, open_c, xy_c_val,
                                                                 position_list_i, grids_rz, xy_index, radius_list,
                                                                 grid_coordinates)
            matrix = convolve3d(processed_con_matrix, kernel, xy_c_val, z_c_val)
            position_con[position_index, fluo_file_index // multiple] = matrix[mid - 1, mid - 1, mid - 1]
        print(f"{nano_file}结束")
    return position_con

# 用这个，最快了
def temporal_distribution_v5(dirname, kernel, xy_c_val, z_c_val, position_list, multiple):
    print("temporal_distribution开始")
    nano_path = "E:\\github_test\\nano_spark\\DATA\\result\\NANO"
    open_path = "E:\\github_test\\nano_spark\\DATA\\result\\OPEN"
    # 初始化各点Ca浓度
    nano_filenames = os.listdir(f"{nano_path}\\{dirname}")  # 纳米空间目前所有步数的浓度文件
    open_filenames = os.listdir(f"{open_path}\\{dirname}")  # 开放空间目前所有步数的浓度文件
    fluo_dir_list_len = len(nano_filenames)  # 当前生成文件数
    position_list_len = len(position_list)  # 需要检测的点个数
    # 哪个文件，该文件取四个点 [[0, 0], [100, 0], [300, 0], [400, 0]]
    position_con = np.empty((position_list_len, fluo_dir_list_len // multiple + 1))
    max_relations = position_list[position_list_len - 1]
    grids_rz = np.load(f"grids_rz/grids_rz_v2_({max_relations[0]},{max_relations[1]}).npy")  # 400,0
    xy_index = np.load(f"xy_list/xy_index.npy")
    new_open_r = np.flipud(point_scatter(295, 0, 5, k=2, positive=False))
    radius_list = np.concatenate((new_open_r, (300,)))
    grid_coordinates = np.loadtxt("../config/open/open_grid_coordinates.csv", delimiter=",")
    # 创建一个哈希表（字典）用于存储坐标值与索引的映射
    coordinates_dict = {(z, r): i for i, (z, r) in enumerate(grid_coordinates)}
    # 文件数，每隔10个（或100个）文件计算一次
    for fluo_file_index in range(0, fluo_dir_list_len, multiple):
        nano_file = nano_filenames[fluo_file_index]
        open_file = open_filenames[fluo_file_index]
        print(f"{nano_file}开始")
        nano_c = np.loadtxt(f"{nano_path}\\{dirname}\\{nano_file}")
        open_c = np.loadtxt(f"{open_path}\\{dirname}\\{open_file}")
        nano_processed = nano_process_concentration(nano_c, xy_c_val)
        # 4
        for position_index in range(position_list_len):
            position_list_i = position_list[position_index]
            # 记录程序开始时间
            start_time = time.time()
            processed_con_matrix, mid = process_concentration_v5(nano_processed, open_c, xy_c_val,
                                                                 position_list_i, grids_rz, xy_index, radius_list,
                                                                 coordinates_dict)
            # 记录程序结束时间
            end_time = time.time()
            # 计算程序运行时间
            elapsed_time = end_time - start_time
            # print(f"程序运行时间: {elapsed_time:.6f} 秒")
            elapsed_time_ms = elapsed_time * 1000  # 将秒转换为毫秒
            print(f"程序process_concentration_v4运行时间: {elapsed_time_ms:.3f} 毫秒")

            # 记录程序开始时间
            start_time = time.time()
            matrix = convolve3d(processed_con_matrix, kernel, xy_c_val, z_c_val)
            # 记录程序结束时间
            end_time = time.time()
            # 计算程序运行时间
            elapsed_time = end_time - start_time
            # print(f"程序运行时间: {elapsed_time:.6f} 秒")
            elapsed_time_ms = elapsed_time * 1000  # 将秒转换为毫秒
            print(f"程序convolve3d运行时间: {elapsed_time_ms:.3f} 毫秒")

            position_con[position_index, fluo_file_index // multiple] = matrix[mid - 1, mid - 1, mid - 1]
        print(f"{nano_file}结束")
    return position_con


def optical_blurring(dirname, position_list):
    kernel = np.load("kernel/kernel_v1.npy", allow_pickle=True)
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
    position_con = temporal_distribution_v5(dirname, kernel, xy_c_val, 0, position_list, multiple)
    # new_position_con = shrink_files(position_con, multiple)
    # 2
    # for position_index in range(len(position_con)):
    #     np.savetxt(
    #         f"{result_set}/{dirname}_psf_({position_list[position_index][0]},{position_list[position_index][1]}).csv"
    #         , position_con[position_index])
    # del position_con


def main():
    dirname = "CaF"
    my_position_list = [[0, 0], [100, 0], [300, 0], [400, 0]]
    optical_blurring(dirname, my_position_list)


if __name__ == "__main__":
    main()
