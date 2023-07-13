import numpy as np
from nano_spark.nano_parameters import *
from tool.cal_bcnl import cal_elements
from scipy.ndimage import convolve1d
from nano_spark.open.open_parameters import *
from tool.generate_open_grid import point_scatter


# 三维卷积
def convolve3d(original_matrix, kernel, xy_c_val, z_c_val):
    matrix = original_matrix.copy()
    c_val = xy_c_val
    i = 0
    print("convolve3d开始")
    for kernel_i in kernel:
        if i == 2:
            c_val = z_c_val
        # 在第一个维度上进行一维卷积操作
        result = convolve1d(matrix, kernel_i, axis=i, mode='constant', cval=c_val)
        matrix = result
        i = i + 1
    return matrix


def nano_process_concentration(original_concentration, c_val):
    concentration = np.full((601, 601, 16), c_val, dtype=float)
    nods = np.loadtxt(nod_file_name, dtype=int) - 1
    relations = np.load("../optical_blurring/refined_relations.npy")
    single_area, control_area, near_triangle, index_in_triangle, nix_multiply_l, niy_multiply_l, a_arr, b_arr, c_arr, nmax, total_area = cal_elements(
        nano_grid_file_name, nod_file_name)
    print("nano_process_concentration开始")
    for i in range(601):
        for j in range(601):
            x = i - 300
            y = j - 300
            relation = relations[i, j, 0]
            c_id = relations[i, j, 1]
            if relation == 1:
                concentration[i, j, :] = original_concentration[c_id]
            elif relation == 2:
                approximate_concentration = 0.0
                for k in range(0, 3):
                    nod_k = nods[c_id, k]
                    f_k = original_concentration[nod_k]
                    a_k = a_arr[c_id, k]
                    b_k = b_arr[c_id, k]
                    c_k = c_arr[c_id, k]
                    approximate_concentration = approximate_concentration + f_k * (a_k + b_k * x + c_k * y)
                concentration[i, j, :] = approximate_concentration
    return concentration


def find_point_index(array, point):
    # 在数组中寻找匹配的索引
    index = np.where((array == point).all(axis=1))

    if index[0].size > 0:
        return index[0][0]
    # else:
    #     return None


def find_multiple_dimensions_index(array, point):
    r, z = point[0], point[1]

    # 在数组中寻找匹配的索引
    indices = np.where((array[:, :, :, 0] == r) & (array[:, :, :, 1] == z))

    i, j, k = indices[0][0], indices[1][0], indices[2][0]
    return i, j, k


# 插值计算
def interpolation_calculation(lower_r, upper_r, lower_z, upper_z, radius, height, original_concentration,
                              grid_coordinates):
    r1, z1 = lower_r, lower_z  # 左下
    r2, z2 = upper_r, lower_z  # 右下
    r3, z3 = upper_r, upper_z  # 右上
    r4, z4 = lower_r, upper_z  # 左上
    alpha = (radius - lower_r) / (upper_r - lower_r)
    beta = (height - lower_z) / (upper_z - lower_z)
    index1 = find_point_index(grid_coordinates, [z1, r1])
    index2 = find_point_index(grid_coordinates, [z2, r2])
    index3 = find_point_index(grid_coordinates, [z3, r3])
    index4 = find_point_index(grid_coordinates, [z4, r4])
    v1, v2, v3, v4 = original_concentration[index1], original_concentration[index2], original_concentration[index3], \
                     original_concentration[index4],
    interpolated_value = (1 - alpha) * (1 - beta) * v1 + alpha * (1 - beta) * v2 + alpha * beta * v3 + (
            1 - alpha) * beta * v4
    return interpolated_value


def open_process_concentration(original_concentration, c_val):
    grid_coordinates = np.loadtxt("../config/open/open_grid_coordinates.csv", delimiter=",")
    grids_rz = np.load("grids_rz.npy")
    points_rz = np.load("points_rz.npy")

    xy_start = -500
    z_start = -500
    xy_end = 500
    z_end = 500
    # 每隔1nm取一个点
    interval = 1
    d1, d2, d3 = (xy_end - xy_start) // interval + 1, (xy_end - xy_start) // interval + 1, (
            z_end - z_start) // interval + 1  # 1001 1001 1001
    concentration = np.full((d1, d2, d3), dtype=float, fill_value=c_val)  # 1001 1001 1001
    for i in range(d1):
        for j in range(d2):
            for k in range(d3):
                lower_r, upper_r, lower_z, upper_z = grids_rz[i, j, k]
                if lower_r != -1.0 and upper_r != -1.0 and lower_z != -1.0 and upper_z != -1.0:
                    radius, height = points_rz[i, j, k]
                    # 对浓进行插值计算
                    concentration[i, j, k] = interpolation_calculation(lower_r, upper_r, lower_z, upper_z, radius,
                                                                       height,
                                                                       original_concentration,
                                                                       grid_coordinates)
    return concentration


def process_concentration(nano_original_concentration, open_original_concentration, c_val):
    points_rz = np.load("points_rz.npy")
    nano_concentration = nano_process_concentration(nano_original_concentration, c_val)
    open_concentration = open_process_concentration(open_original_concentration, c_val)

    xy_start = -500
    z_start = -500
    xy_end = 500
    z_end = 500
    # 每隔1nm取一个点
    interval = 1
    d1, d2, d3 = (xy_end - xy_start) // interval + 1, (xy_end - xy_start) // interval + 1, (
            z_end - z_start) // interval + 1  # 1001 1001 1001
    concentration = np.full((d1, d2, d3), dtype=float, fill_value=c_val)  # 1001 1001 1001

    xy_arange = np.arange(xy_start, xy_end + 1)  # -500,501
    xy_squared = np.square(xy_arange)
    z_arange = np.abs(np.arange(z_start, z_end + 1))  # -500,501

    for i in range(d1):
        x = xy_arange[i]
        x2 = xy_squared[i]
        for j in range(d2):
            y = xy_arange[j]
            y2 = xy_squared[j]
            radius = np.sqrt(x2 + y2)
            for k in range(d3):
                height = z_arange[k]
                # 点在开放空间中
                if height >= 7.5 or radius >= 300:
                    a, b, c = find_multiple_dimensions_index(points_rz, [radius, height])
                    concentration[i, j, k] = open_concentration[a, b, c]
                # 点在纳米空间中
                else:
                    concentration[i, j, k] = nano_concentration[x + 300, y + 300, height]
    return concentration
