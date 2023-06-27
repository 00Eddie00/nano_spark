from nano_parameters import *
import numpy as np
from tool.cal_bcnl import *


def cal_dye_and_buffers(f, caf, cag, cab1, cab2, cab3, cab4):
    j_fdye = np.zeros(NP)
    j_gdye = np.zeros(NP)
    j_1 = np.zeros(NP)
    j_2 = np.zeros(NP)
    j_3 = np.zeros(NP)
    j_4 = np.zeros(NP)
    new_cag = np.zeros(NP)
    new_cab1 = np.zeros(NP)
    new_cab2 = np.zeros(NP)
    new_cab3 = np.zeros(NP)
    new_cab4 = np.zeros(NP)
    for i in range(0, NP):
        j_fdye[i] = -K_F3_PLUS * f[i] * (F3_T - caf[i]) + K_F3_MINUS * caf[i]
        j_gdye[i] = -K_GCaMP6f_PLUS * f[i] * (GCaMP6f_T - cag[i]) + K_GCaMP6f_MINUS * cag[i]
        j_1[i] = -K_Calmodulin_PLUS * f[i] * (Calmodulin_T - cab1[i]) + K_Calmodulin_MINUS * cab1[i]
        j_2[i] = -K_TroponinC_PLUS * f[i] * (TroponinC_T - cab2[i]) + K_TroponinC_MINUS * cab2[i]
        j_3[i] = -K_SR_PLUS * f[i] * (SR_T - cab3[i]) + K_SR_MINUS * cab3[i]
        j_4[i] = -K_SL_PLUS * f[i] * (SL_T - cab4[i]) + K_SL_MINUS * cab4[i]

    for i in range(0, NP):
        new_cag[i] = (-j_gdye[i]) * DT + cag[i]
        new_cab1[i] = (-j_1[i]) * DT + cab1[i]
        new_cab2[i] = (-j_2[i]) * DT + cab2[i]
        new_cab3[i] = (-j_3[i]) * DT + cab3[i]
        new_cab4[i] = (-j_4[i]) * DT + cab4[i]
    return j_fdye, j_gdye, j_1, j_2, j_3, j_4, new_cag, new_cab1, new_cab2, new_cab3, new_cab4


def nano_calculation_f(k_ryr, f, caf, cag, cab1, cab2, cab3, cab4, ca_jsr, c_ca_out):
    # 计算染料和缓冲物
    j_fdye, j_gdye, j_1, j_2, j_3, j_4, new_cag, new_cab1, new_cab2, new_cab3, new_cab4 = cal_dye_and_buffers(f, caf,
                                                                                                              cag, cab1,
                                                                                                              cab2,
                                                                                                              cab3,
                                                                                                              cab4)
    # 计算三角形系数
    single_area, control_area, near_triangle, index_in_triangle, nix_multiply_l, niy_multiply_l, a_arr, b_arr, c_arr, nmax, total_area = cal_elements(
        nano_grid_file_name, nod_file_name)

    # 上一次迭代得到的值
    last = np.zeros(NP)
    # 此时刻计算最终结果
    temp = np.zeros(NP)
    for j in range(0, 10):
        for i in range(0, NP):
            item_1 = (j_fdye[i] + j_gdye[i] + j_1[i] + j_2[i] + j_3[i] + j_4[i]) * control_area[i]
            item_2, item_3, item_4, item_5, item_6, item_7 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            # nmax=9
            for k in range(0, nmax):
                # 三角形编号
                near_triangle_index = near_triangle[i, k]
                # 判断有无相邻三角形
                if near_triangle_index != -1:
                    # 该点在该三角形（near_triangle_index）中的编号
                    point_index = index_in_triangle[i, k]
                    point2 = (point_index + 1) % 3
                    point3 = (point_index + 2) % 3
                    # 该点所对应的边的两点
                    nod2 = nods[near_triangle_index, point2]
                    nod3 = nods[near_triangle_index, point3]
                    in_boundary, out_boundary, inner = judge_point(near_triangle_index, nods, npoch)

                    # 求导部分的项没有边界边
                    if npoch[nod2] == B_INNER or npoch[nod3] == B_INNER:
                        f2i, f3i = f[nod2], f[nod3]
                        if j != 0:
                            # f2i, f3i = last[nod2], last[nod3]
                            f2i, f3i = (f[nod2] + last[nod2]) / 2.0, (f[nod3] + last[nod3]) / 2.0

                        item_2 = item_2 + D_CA * (
                                (f2i * b_arr[near_triangle_index, point2] + f3i * b_arr[near_triangle_index, point3]) *
                                nix_multiply_l[near_triangle_index, point_index] + (
                                        f2i * c_arr[near_triangle_index, point2] + f3i * c_arr[
                                    near_triangle_index, point3]) * niy_multiply_l[
                                    near_triangle_index, point_index])
                        item_3 = item_3 + D_CA * (b_arr[near_triangle_index, point_index] * nix_multiply_l[
                            near_triangle_index, point_index] + c_arr[near_triangle_index, point_index] *
                                                  niy_multiply_l[
                                                      near_triangle_index, point_index])

                    f2jk, f3jk = f[nod2], f[nod3]
                    if j != 0:
                        f2jk, f3jk = last[nod2], last[nod3]
                    fjk = (f2jk + f3jk) / 3.0
                    if len(in_boundary) == 2 and k_ryr != 0:
                        length_j = cal_length(in_boundary, grids)
                        item_4 = item_4 + k_ryr * (ca_jsr - fjk) * length_j
                        item_5 = item_5 - k_ryr * length_j / 3.0
                    elif len(out_boundary) == 2:
                        length_k = cal_length(out_boundary, grids)
                        item_6 = item_6 + (D_CA / Delta_r) * (c_ca_out - fjk) * length_k
                        item_7 = item_7 - (D_CA / Delta_r) * length_k / 3.0

            temp[i] = ((item_2 + item_1 + item_4 + item_6) * DT + control_area[i] * f[i]) / (
                    control_area[i] - (item_3 + item_5 + item_7) * DT)
        # 此次迭代的值保存到last
        last = np.copy(temp)
    # 这一步的值保存到new_caf
    new_f = np.copy(last)
    return new_f, new_cag, new_cab1, new_cab2, new_cab3, new_cab4


def nano_calculation_g2(f, caf, c_caf_out):
    j_fdye = np.zeros(NP)
    for i in range(0, NP):
        j_fdye[i] = -K_F3_PLUS * f[i] * (F3_T - caf[i]) + K_F3_MINUS * caf[i]
    single_area, control_area, near_triangle, index_in_triangle, nix_multiply_l, niy_multiply_l, a_arr, b_arr, c_arr, nmax, total_area = cal_elements(
        nano_grid_file_name, nod_file_name)

    last = np.zeros(NP)
    temp = np.zeros(NP)
    for j in range(0, 10):
        for i in range(0, NP):
            item_1 = -j_fdye[i] * control_area[i]
            item_2, item_3, item_4, item_5 = 0.0, 0.0, 0.0, 0.0
            for k in range(0, nmax):
                # 三角形编号
                near_triangle_index = near_triangle[i, k]
                # 判断有无相邻三角形
                if near_triangle_index != -1:
                    # 该点在该三角形（near_triangle_index）中的编号
                    point_index = index_in_triangle[i, k]
                    point2 = (point_index + 1) % 3
                    point3 = (point_index + 2) % 3
                    nod2 = nods[near_triangle_index, point2]
                    nod3 = nods[near_triangle_index, point3]
                    in_boundary, out_boundary, inner = judge_point(near_triangle_index, nods, npoch)
                    if npoch[nod2] == B_INNER or npoch[nod3] == B_INNER:
                        g2i, g3i = caf[nod2], caf[nod3]
                        if j != 0:
                            # g2i, g3i = last[nod2], last[nod3]
                            g2i, g3i = (caf[nod2] + last[nod2]) / 2.0, (caf[nod3] + last[nod3]) / 2.0
                        item_2 = item_2 + D_CAF * ((g2i * b_arr[near_triangle_index, point2] + g3i * b_arr[
                            near_triangle_index, point3]) * nix_multiply_l[near_triangle_index, point_index] + (
                                                           g2i * c_arr[near_triangle_index, point2] + g3i * c_arr[
                                                       near_triangle_index, point3]) * niy_multiply_l[
                                                       near_triangle_index, point_index])
                        item_3 = item_3 + D_CAF * (b_arr[near_triangle_index, point_index] * nix_multiply_l[
                            near_triangle_index, point_index] + c_arr[near_triangle_index, point_index] *
                                                   niy_multiply_l[near_triangle_index, point_index])
                    g2k, g3k = caf[nod2], caf[nod3]
                    if j != 0:
                        g2k, g3k = last[nod2], last[nod3]
                    gk = (g2k + g3k) / 3.0
                    if len(out_boundary) == 2:
                        length_k = cal_length(out_boundary, grids)
                        item_4 = item_4 + (D_CAF / Delta_r) * (c_caf_out - gk) * length_k
                        item_5 = item_5 - (D_CAF / Delta_r) * length_k / 3.0
            temp[i] = ((item_2 + item_1 + item_4) * DT + control_area[i] * caf[i]) / (
                    control_area[i] - (item_3 + item_5) * DT)
        # 此次迭代的值保存到last
        last = np.copy(temp)
    # 这一步的值保存到new_caf
    new_caf = np.copy(last)
    return new_caf
