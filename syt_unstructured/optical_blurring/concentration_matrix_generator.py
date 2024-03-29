import numpy as np
from nano_spark.nano_parameters import *
from tool.cal_bcnl import cal_elements
from scipy.ndimage import convolve1d
from tool.generate_nano_con import *
import math
from ob_parameters import half_length


# 三维卷积
def convolve3d(original_matrix, kernel, xy_c_val, z_c_val):
    # matrix = original_matrix.copy()
    c_val = xy_c_val
    i = 0
    print("convolve3d开始")
    for kernel_i in kernel:
        if i == 2:
            c_val = z_c_val
        # 在第一个维度上进行一维卷积操作
        result = convolve1d(original_matrix, kernel_i, axis=i, mode='constant', cval=c_val)
        original_matrix = result
        i = i + 1
    return original_matrix


def nano_process_concentration(original_concentration, c_val):
    concentration = np.full((601, 601, 16), c_val, dtype=float)
    nods = np.loadtxt(nod_file_name, dtype=int) - 1
    relations = np.load("relation/refined_relations.npy")
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


def find_multiple_dimensions_index(array, point):
    r, z = point[0], point[1]

    # 在数组中寻找匹配的索引
    indices = np.where((array[:, :, :, 0] == r) & (array[:, :, :, 1] == z))

    i, j, k = indices[0][0], indices[1][0], indices[2][0]
    return i, j, k


# 插值计算
def interpolation_calculation(lower_r, upper_r, lower_z, upper_z, radius, height, open_original_concentration,
                              grid_coordinates, nano_processed):
    r1, z1 = lower_r, lower_z  # 左下
    r2, z2 = upper_r, lower_z  # 右下
    r3, z3 = upper_r, upper_z  # 右上
    r4, z4 = lower_r, upper_z  # 左上
    alpha = (radius - lower_r) / (upper_r - lower_r)
    beta = (height - lower_z) / (upper_z - lower_z)
    if lower_z == 7.5 and radius < 300:
        index1 = np.int32(math.floor(radius))
        index2 = np.int32(math.floor(radius)) + 1
        index3 = find_point_index(grid_coordinates, [z3, r3])
        index4 = find_point_index(grid_coordinates, [z4, r4])
        # 有一丢丢问题
        v1, v2, v3, v4 = nano_processed[index1 + 300, 0, 8], nano_processed[index2 + 300, 0, 8], \
                         open_original_concentration[
                             index3], open_original_concentration[index4]
    elif radius == 300 and height <= 13.0:
        index1 = np.int32(r1)
        index2 = find_point_index(grid_coordinates, [z2, r2])
        index3 = find_point_index(grid_coordinates, [z3, r3])
        v1, v2, v3 = nano_processed[index1 + 300, 0, 8], open_original_concentration[index2], \
                     open_original_concentration[index3]
        if height <= 7.0:
            index4 = np.int32(r4)
            v4 = nano_processed[index4 + 300, 0, 8]
        else:
            index4 = find_point_index(grid_coordinates, [z4, r4])
            v4 = open_original_concentration[index4]
    else:
        index1 = find_point_index(grid_coordinates, [z1, r1])
        index2 = find_point_index(grid_coordinates, [z2, r2])
        index3 = find_point_index(grid_coordinates, [z3, r3])
        index4 = find_point_index(grid_coordinates, [z4, r4])
        v1, v2, v3, v4 = open_original_concentration[index1], open_original_concentration[index2], \
                         open_original_concentration[index3], \
                         open_original_concentration[index4]

    interpolated_value = (1 - alpha) * (1 - beta) * v1 + alpha * (1 - beta) * v2 + alpha * beta * v3 + (
            1 - alpha) * beta * v4
    return interpolated_value


# 点太多数据量太大
def open_process_concentration(original_concentration, c_val, position_list_i):
    grid_coordinates = np.loadtxt("../config/open/open_grid_coordinates.csv", delimiter=",")
    grids_rz = np.load(f"grids_rz_v2_({position_list_i[0]},{position_list_i[1]}).npy")
    points_rz = np.load(f"points_rz_v2_({position_list_i[0]},{position_list_i[1]}).npy")
    print("open_process_concentration 开始")
    dis = position_list_i[0]  # 0,100,300,400
    x_start = dis - 500
    y_start = -500
    z_start = -500
    x_end = dis + 500
    y_end = 500
    z_end = 500
    # 每隔1nm取一个点
    interval = 1
    d1, d2, d3 = (x_end - x_start) // interval + 1, (y_end - y_start) // interval + 1, (
            z_end - z_start) // interval + 1  # 1001 1001 1001
    concentration = np.full((d1, d2, d3), fill_value=c_val)  # 1001 1001 1001
    for i in range(d1):
        a = np.abs(i - (500 - dis))
        for j in range(d2):
            b = np.abs(j - 500)
            for k in range(d3):
                c = np.abs(k - 500)
                lower_r, upper_r, lower_z, upper_z = grids_rz[a, b, c]
                if lower_r != -1.0 and upper_r != -1.0 and lower_z != -1.0 and upper_z != -1.0:
                    radius, height = points_rz[a, b, c]
                    # 对浓度进行插值计算
                    concentration[i, j, k] = interpolation_calculation(lower_r, upper_r, lower_z, upper_z, radius,
                                                                       height,
                                                                       original_concentration,
                                                                       grid_coordinates)
    return concentration


# 点太多数据量太大
def open_process_concentration_v2(nano_processed, original_concentration, c_val, position_list_i, grids_rz, points_rz):
    grid_coordinates = np.loadtxt("../config/open/open_grid_coordinates.csv", delimiter=",")
    print("open_process_concentration 开始")
    dis = position_list_i[0]  # 400
    x_start = 0
    y_start = 0
    z_start = 0
    x_end = dis + 500
    y_end = 500
    z_end = 500
    # 每隔1nm取一个点
    interval = 1
    d1, d2, d3 = (x_end - x_start) // interval + 1, (y_end - y_start) // interval + 1, (
            z_end - z_start) // interval + 1  # 901 501 501
    concentration = np.full((d1, d2, d3), fill_value=c_val)  # 901 501 501
    for i in range(d1):
        for j in range(d2):
            for k in range(d3):
                lower_r, upper_r, lower_z, upper_z = grids_rz[i, j, k]
                if lower_r != -1.0 and upper_r != -1.0 and lower_z != -1.0 and upper_z != -1.0:
                    radius, height = points_rz[i, j, k]
                    # 对浓度进行插值计算
                    concentration[i, j, k] = interpolation_calculation(lower_r, upper_r, lower_z, upper_z, radius,
                                                                       height,
                                                                       original_concentration,
                                                                       grid_coordinates, nano_processed)
    return concentration


def interpolation_calculation_v3(lower_r, upper_r, lower_z, upper_z, radius, height, open_original_concentration,
                                 grid_coordinates, nano_processed):
    r1, z1 = lower_r, lower_z  # 左下
    r2, z2 = upper_r, lower_z  # 右下
    r3, z3 = upper_r, upper_z  # 右上
    r4, z4 = lower_r, upper_z  # 左上
    alpha = (radius - lower_r) / (upper_r - lower_r)
    beta = (height - lower_z) / (upper_z - lower_z)
    # 该网格底部为上侧交界处（开放网格中没有）
    if lower_z == 7.5 and radius < 300:
        index1 = np.int32(lower_r)
        index2 = np.int32(upper_r)
        index3 = find_point_index(grid_coordinates, [z3, r3])
        index4 = find_point_index(grid_coordinates, [z4, r4])
        v1, v2, v3, v4 = nano_processed[index1 + 300, 300, 8], nano_processed[index2 + 300, 300, 8], \
                         open_original_concentration[index3], open_original_concentration[index4]
        # 该网格左侧在纳米空间内，右侧在交界处
    elif radius == 300 and height <= 7.5:
        index1 = np.int32(lower_r)
        index2 = find_point_index(grid_coordinates, [z2, r2])
        index3 = find_point_index(grid_coordinates, [z3, r3])
        index4 = np.int32(lower_r)
        v1, v2, v3, v4 = nano_processed[index1 + 300, 300, 8], open_original_concentration[index2], \
                         open_original_concentration[index3], nano_processed[index4 + 300, 300, 8]
    # 该网格左下在上侧交界处
    elif radius == 300 and height <= 13.0:
        index1 = np.int32(lower_r)
        index2 = find_point_index(grid_coordinates, [z2, r2])
        index3 = find_point_index(grid_coordinates, [z3, r3])
        index4 = find_point_index(grid_coordinates, [z4, r4])
        v1, v2, v3, v4 = nano_processed[index1 + 300, 300, 8], open_original_concentration[index2], \
                         open_original_concentration[index3], open_original_concentration[index4]
    # 四个点的浓度都在开放空间文件中
    else:
        index1 = find_point_index(grid_coordinates, [z1, r1])
        index2 = find_point_index(grid_coordinates, [z2, r2])
        index3 = find_point_index(grid_coordinates, [z3, r3])
        index4 = find_point_index(grid_coordinates, [z4, r4])
        v1, v2, v3, v4 = open_original_concentration[index1], open_original_concentration[index2], \
                         open_original_concentration[index3], \
                         open_original_concentration[index4]

    interpolated_value = (1 - alpha) * (1 - beta) * v1 + alpha * (1 - beta) * v2 + alpha * beta * v3 + (
            1 - alpha) * beta * v4
    return interpolated_value


# 比V4慢可能也没有V4准确
def process_concentration_v3(nano_processed, original_concentration, c_val, position_list_i, grids_rz):
    grid_coordinates = np.loadtxt("../config/open/open_grid_coordinates.csv", delimiter=",")
    print("process_concentration 开始")
    dis = position_list_i[0]  # 0,100,300,400
    x_arr, y_arr, z_arr, process_points = generate_interval(dis, half_length)
    # xyz轴方向上点的个数d1<=900, d2<=500, d3<=500
    d1, d2, d3 = len(x_arr), len(y_arr), len(z_arr)
    # xy的平方，方便后续计算半径
    x_square, y_square = np.square(x_arr), np.square(y_arr)
    # xyz轴方向上点的坐标的绝对值，因为对称性
    x_abs, y_abs, z_abs = np.abs(x_arr), np.abs(y_arr), np.abs(z_arr)
    # 最终生成的三维矩阵
    concentration = np.full((d1, d2, d3), fill_value=c_val)
    for i in range(d1):
        x2 = x_square[i]
        a = x_abs[i]
        x_i = x_arr[i]
        for j in range(d2):
            y2 = y_square[j]
            b = y_abs[j]
            y_j = y_arr[j]
            for k in range(d3):
                height = z_abs[k]
                radius = np.sqrt(x2 + y2)
                c = z_abs[k]
                if height >= 7.5 or radius >= 300:
                    lower_r, upper_r, lower_z, upper_z = grids_rz[np.int32(a), np.int32(b), np.int32(c)]
                    concentration[i, j, k] = interpolation_calculation_v3(lower_r, upper_r, lower_z, upper_z,
                                                                          radius,
                                                                          height,
                                                                          original_concentration,
                                                                          grid_coordinates, nano_processed)
                else:
                    concentration[i, j, k] = nano_processed[np.int32(x_i) + 300, np.int32(y_j) + 300, 8]
    return concentration, len(process_points)


# 不对，点太多数据量太大
def process_concentration(nano_processed, open_original_concentration, c_val, position_list_i):
    # points_rz = np.load(f"points_rz_v2_({position_list_i[0]},{position_list_i[1]}).npy")
    # nano_concentration = nano_process_concentration(nano_original_concentration, c_val)
    open_concentration = open_process_concentration(open_original_concentration, c_val, position_list_i)
    dis = position_list_i[0]  # 0,100,300,400
    x_start = dis - 500
    y_start = -500
    z_start = -500
    x_end = dis + 500
    y_end = 500
    z_end = 500
    # 每隔1nm取一个点
    interval = 1
    d1, d2, d3 = (x_end - x_start) // interval + 1, (y_end - y_start) // interval + 1, (
            z_end - z_start) // interval + 1  # 1001 1001 1001
    concentration = np.full((d1, d2, d3), dtype=np.float32, fill_value=c_val)  # 1001 1001 1001

    interval = 1
    x_arange = np.arange(x_start, x_end + 1, interval)
    y_arange = np.arange(y_start, y_end + 1, interval)  # -500,501
    z_arange = np.arange(z_start, z_end + 1, interval)
    x_squared = np.square(x_arange)
    y_squared = np.square(y_arange)
    z_arange_abs = np.abs(z_arange)  # -500,501

    for i in range(d1):
        x = x_arange[i]
        x2 = x_squared[i]
        for j in range(d2):
            y = y_arange[j]
            y2 = y_squared[j]
            radius = np.sqrt(x2 + y2)
            for k in range(d3):
                height = z_arange_abs[k]
                # 点在开放空间中
                if height >= 7.5 or radius >= 300:
                    # a, b, c = find_multiple_dimensions_index(points_rz, [radius, height])
                    concentration[i, j, k] = open_concentration[i, j, k]
                # 点在纳米空间中
                else:
                    concentration[i, j, k] = nano_processed[x + 300, y + 300, height]
    del open_concentration
    return concentration


# 不对，点太多数据量太大
def process_concentration_v2(nano_processed, open_processed_whole, c_val, position_list_i, points_rz):
    dis = position_list_i[0]  # 0,100,300,400
    x_start = dis - 500
    y_start = -500
    z_start = -500
    x_end = dis + 500
    y_end = 500
    z_end = 500
    # 每隔1nm取一个点
    interval = 1
    d1, d2, d3 = (x_end - x_start) // interval + 1, (y_end - y_start) // interval + 1, (
            z_end - z_start) // interval + 1  # 1001 1001 1001
    concentration = np.full((d1, d2, d3), fill_value=c_val)  # 1001 1001 1001

    interval = 1
    x_arange = np.arange(x_start, x_end + 1, interval)
    y_arange = np.arange(y_start, y_end + 1, interval)  # -500,501
    z_arange = np.arange(z_start, z_end + 1, interval)
    x_squared = np.square(x_arange)
    y_squared = np.square(y_arange)
    z_arange_abs = np.abs(z_arange)  # -500,501

    for i in range(d1):
        x = x_arange[i]
        x2 = x_squared[i]
        for j in range(d2):
            y = y_arange[j]
            y2 = y_squared[j]
            radius = np.sqrt(x2 + y2)
            for k in range(d3):
                height = z_arange_abs[k]
                # 点在开放空间中
                if height >= 7.5 or radius >= 300:
                    a, b, c = find_multiple_dimensions_index(points_rz, [radius, height])
                    concentration[i, j, k] = open_processed_whole[a, b, c]
                # 点在纳米空间中
                else:
                    concentration[i, j, k] = nano_processed[x + 300, y + 300, height]
    return concentration


def same_radius_avg(nano_processed, r, xy_index, radius_list):
    indexes = np.where(radius_list == r)
    a = indexes[0][0]
    r_xy_index = xy_index[a]
    con_sum = 0.0
    count = 0
    for q in range(len(r_xy_index)):
        if r_xy_index[q][0] == -1:
            break
        i = r_xy_index[q][0]
        j = r_xy_index[q][1]
        con_sum = con_sum + nano_processed[i, j, 8]
        count = count + 1
    con_avg = con_sum / count
    return con_avg


def interpolation_calculation_v4(lower_r, upper_r, lower_z, upper_z, radius, height, open_original_concentration,
                                 grid_coordinates, nano_processed, xy_index, radius_list):
    r1, z1 = lower_r, lower_z  # 左下
    r2, z2 = upper_r, lower_z  # 右下
    r3, z3 = upper_r, upper_z  # 右上
    r4, z4 = lower_r, upper_z  # 左上
    alpha = (radius - lower_r) / (upper_r - lower_r)
    beta = (height - lower_z) / (upper_z - lower_z)
    # 该网格底部为上侧交界处（开放网格中没有）
    if lower_z == 7.5 and radius < 300:
        # index1 = np.int32(lower_r)
        # index2 = np.int32(upper_r)
        index3 = find_point_index(grid_coordinates, [z3, r3])
        index4 = find_point_index(grid_coordinates, [z4, r4])
        v1, v2, v3, v4 = same_radius_avg(nano_processed, lower_r, xy_index, radius_list), same_radius_avg(
            nano_processed, upper_r, xy_index, radius_list), \
                         open_original_concentration[index3], open_original_concentration[index4]
    # 该网格左侧在纳米空间内，右侧在交界处
    elif radius == 300 and height <= 7.5:
        # index1 = np.int32(lower_r)
        index2 = find_point_index(grid_coordinates, [z2, r2])
        index3 = find_point_index(grid_coordinates, [z3, r3])
        # index4 = np.int32(lower_r)
        vv = same_radius_avg(nano_processed, lower_r, xy_index, radius_list)
        v1, v2, v3, v4 = vv, open_original_concentration[index2], open_original_concentration[index3], vv
    # 该网格左下在上侧交界处
    elif radius == 300 and height <= 13.0:
        # index1 = np.int32(lower_r)
        index2 = find_point_index(grid_coordinates, [z2, r2])
        index3 = find_point_index(grid_coordinates, [z3, r3])
        index4 = find_point_index(grid_coordinates, [z4, r4])
        v1, v2, v3, v4 = same_radius_avg(nano_processed, lower_r, xy_index, radius_list), open_original_concentration[
            index2], \
                         open_original_concentration[index3], open_original_concentration[index4]
    # 四个点的浓度都在开放空间文件中
    else:
        index1 = find_point_index(grid_coordinates, [z1, r1])
        index2 = find_point_index(grid_coordinates, [z2, r2])
        index3 = find_point_index(grid_coordinates, [z3, r3])
        index4 = find_point_index(grid_coordinates, [z4, r4])
        v1, v2, v3, v4 = open_original_concentration[index1], open_original_concentration[index2], \
                         open_original_concentration[index3], \
                         open_original_concentration[index4]

    interpolated_value = (1 - alpha) * (1 - beta) * v1 + alpha * (1 - beta) * v2 + alpha * beta * v3 + (
            1 - alpha) * beta * v4
    return interpolated_value


# 一般用这个，最准确但是很慢60000步要跑6、7天
def process_concentration_v4(nano_processed, original_concentration, c_val, position_list_i, grids_rz, xy_index,
                             radius_list, grid_coordinates):
    print("process_concentration 开始")
    dis = position_list_i[0]  # 0,100,300,400
    x_arr, y_arr, z_arr, process_points = generate_interval(dis, half_length)
    d1, d2, d3 = len(x_arr), len(y_arr), len(z_arr)
    x_square, y_square = np.square(x_arr), np.square(y_arr)
    x_abs, y_abs, z_abs = np.abs(x_arr), np.abs(y_arr), np.abs(z_arr)
    concentration = np.full((d1, d2, d3), fill_value=c_val)
    for i in range(d1):
        x2 = x_square[i]
        a = x_abs[i]
        x_i = x_arr[i]
        for j in range(d2):
            y2 = y_square[j]
            b = y_abs[j]
            y_j = y_arr[j]
            for k in range(d3):
                height = z_abs[k]
                radius = np.sqrt(x2 + y2)
                c = z_abs[k]
                if height >= 7.5 or radius >= 300:
                    lower_r, upper_r, lower_z, upper_z = grids_rz[np.int32(a), np.int32(b), np.int32(c)]
                    concentration[i, j, k] = interpolation_calculation_v4(lower_r, upper_r, lower_z, upper_z,
                                                                          radius,
                                                                          height,
                                                                          original_concentration,
                                                                          grid_coordinates, nano_processed, xy_index,
                                                                          radius_list)
                else:
                    concentration[i, j, k] = nano_processed[np.int32(x_i) + 300, np.int32(y_j) + 300, 8]
    return concentration, len(process_points)


def find_point_index_v2(coordinates_dict, target_coordinates):
    if target_coordinates in coordinates_dict:
        index = coordinates_dict[target_coordinates]
        return index
    else:
        print("Target coordinates not found in the array.")


def interpolation_calculation_v5(lower_r, upper_r, lower_z, upper_z, radius, height, open_original_concentration,
                                 coordinates_dict, nano_processed, xy_index, radius_list):
    r1, z1 = lower_r, lower_z  # 左下
    r2, z2 = upper_r, lower_z  # 右下
    r3, z3 = upper_r, upper_z  # 右上
    r4, z4 = lower_r, upper_z  # 左上
    alpha = (radius - lower_r) / (upper_r - lower_r)
    beta = (height - lower_z) / (upper_z - lower_z)
    # 该网格底部为上侧交界处（开放网格中没有）
    if lower_z == 7.5 and radius < 300:
        # index1 = np.int32(lower_r)
        # index2 = np.int32(upper_r)
        index3 = find_point_index_v2(coordinates_dict, (z3, r3))
        index4 = find_point_index_v2(coordinates_dict, (z4, r4))
        v1, v2, v3, v4 = same_radius_avg(nano_processed, lower_r, xy_index, radius_list), same_radius_avg(
            nano_processed, upper_r, xy_index, radius_list), \
                         open_original_concentration[index3], open_original_concentration[index4]
    # 该网格左侧在纳米空间内，右侧在交界处
    elif radius == 300 and height <= 7.5:
        # index1 = np.int32(lower_r)
        index2 = find_point_index_v2(coordinates_dict, (z2, r2))
        index3 = find_point_index_v2(coordinates_dict, (z3, r3))
        # index4 = np.int32(lower_r)
        vv = same_radius_avg(nano_processed, lower_r, xy_index, radius_list)
        v1, v2, v3, v4 = vv, open_original_concentration[index2], open_original_concentration[index3], vv
    # 该网格左下在上侧交界处
    elif radius == 300 and height <= 13.0:
        # index1 = np.int32(lower_r)
        index2 = find_point_index_v2(coordinates_dict, (z2, r2))
        index3 = find_point_index_v2(coordinates_dict, (z3, r3))
        index4 = find_point_index_v2(coordinates_dict, (z4, r4))
        v1, v2, v3, v4 = same_radius_avg(nano_processed, lower_r, xy_index, radius_list), open_original_concentration[
            index2], \
                         open_original_concentration[index3], open_original_concentration[index4]
    # 四个点的浓度都在开放空间文件中
    else:
        index1 = find_point_index_v2(coordinates_dict, (z1, r1))
        index2 = find_point_index_v2(coordinates_dict, (z2, r2))
        index3 = find_point_index_v2(coordinates_dict, (z3, r3))
        index4 = find_point_index_v2(coordinates_dict, (z4, r4))
        v1, v2, v3, v4 = open_original_concentration[index1], open_original_concentration[index2], \
                         open_original_concentration[index3], \
                         open_original_concentration[index4]

    interpolated_value = (1 - alpha) * (1 - beta) * v1 + alpha * (1 - beta) * v2 + alpha * beta * v3 + (
            1 - alpha) * beta * v4
    return interpolated_value


def process_concentration_v5(nano_processed, original_concentration, c_val, position_list_i, grids_rz, xy_index,
                             radius_list, coordinates_dict):
    print("process_concentration 开始")
    dis = position_list_i[0]  # 0,100,300,400
    x_arr, y_arr, z_arr, process_points = generate_interval(dis, half_length)
    d1, d2, d3 = len(x_arr), len(y_arr), len(z_arr)
    x_square, y_square = np.square(x_arr), np.square(y_arr)
    x_abs, y_abs, z_abs = np.abs(x_arr), np.abs(y_arr), np.abs(z_arr)
    concentration = np.full((d1, d2, d3), fill_value=c_val)
    for i in range(d1):
        x2 = x_square[i]
        a = x_abs[i]
        x_i = x_arr[i]
        for j in range(d2):
            y2 = y_square[j]
            b = y_abs[j]
            y_j = y_arr[j]
            for k in range(d3):
                height = z_abs[k]
                radius = np.sqrt(x2 + y2)
                c = z_abs[k]
                if height >= 7.5 or radius >= 300:
                    lower_r, upper_r, lower_z, upper_z = grids_rz[np.int32(a), np.int32(b), np.int32(c)]
                    concentration[i, j, k] = interpolation_calculation_v5(lower_r, upper_r, lower_z, upper_z,
                                                                          radius,
                                                                          height,
                                                                          original_concentration,
                                                                          coordinates_dict, nano_processed, xy_index,
                                                                          radius_list)
                else:
                    concentration[i, j, k] = nano_processed[np.int32(x_i) + 300, np.int32(y_j) + 300, 8]
    return concentration, len(process_points)
