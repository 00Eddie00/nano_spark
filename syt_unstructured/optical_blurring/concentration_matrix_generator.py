import numpy as np
from nano_spark.nano_parameters import *
from tool.cal_bcnl import cal_elements
from scipy.ndimage import convolve1d


# 三维卷积
def convolve3d(original_matrix, kernel, xy_c_val, z_c_val):
    matrix = original_matrix.copy()
    c_val = xy_c_val
    i = 0
    for kernel_i in kernel:
        if i == 2:
            c_val = z_c_val
        # 在第一个维度上进行一维卷积操作
        result = convolve1d(matrix, kernel_i, axis=i, mode='constant', cval=c_val)
        matrix = result
        i = i + 1
    return matrix


def process_concentration(original_concentration, c_val):
    concentration = np.full((601, 601, 16), c_val, dtype=float)
    nods = np.loadtxt(nod_file_name) - 1
    relations = np.load("../optical_blurring/refined_relations.npy")
    for i in range(601):
        for j in range(601):
            x = i - 300
            y = j - 300
            relation = relations[i, j, 0]
            c_id = relations[i, j, 1]
            if relation == 1:
                concentration[i, j, :] = original_concentration[c_id]
            elif relation == 2:
                single_area, control_area, near_triangle, index_in_triangle, nix_multiply_l, niy_multiply_l, a_arr, b_arr, c_arr, nmax, total_area = cal_elements(
                    nano_grid_file_name, nod_file_name)
                approximate_concentration = 0
                for k in range(3):
                    nod_k = nods[c_id, k]
                    f_k = original_concentration[nod_k]
                    a_k = a_arr[c_id, k]
                    b_k = b_arr[c_id, k]
                    c_k = c_arr[c_id, k]
                    approximate_concentration = approximate_concentration + f_k * (a_k + b_k * x + c_k * y)
                concentration[i, j, :] = approximate_concentration
    return concentration
