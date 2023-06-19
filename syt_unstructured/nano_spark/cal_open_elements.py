from open_parameters import *
from nano_parameters import *
import numpy as np


def neighbors_concentration(point_id, concentration, p32, p12, p23, p21):
    if p32 != -1:
        f32 = concentration[p32]
    else:
        f32 = concentration[point_id]
    if p12 != -1:
        f12 = concentration[p12]
    else:
        f12 = concentration[point_id]
    if p23 != -1:
        f23 = concentration[p23]
    else:
        f23 = concentration[point_id]
    if p21 != -1:
        f21 = concentration[p21]
    else:
        f21 = concentration[point_id]
    return f32, f12, f23, f21


def cal_dye_and_buffers(f, caf, cab1, cab2, cab3, cab4):
    grid_coordinates = np.loadtxt(open_grid_file_name, delimiter=",")  # 点坐标
    point_count = len(grid_coordinates)
    j_fdye = np.zeros(point_count)
    j_1 = np.zeros(point_count)
    j_2 = np.zeros(point_count)
    j_3 = np.zeros(point_count)
    j_4 = np.zeros(point_count)
    new_cab1 = np.zeros(point_count)
    new_cab2 = np.zeros(point_count)
    new_cab3 = np.zeros(point_count)
    new_cab4 = np.zeros(point_count)
    for i in range(0, point_count):
        j_fdye[i] = -K_F3_PLUS * f[i] * (F3_T - caf[i]) + K_F3_MINUS * caf[i]
        j_1[i] = -K_Calmodulin_PLUS * f[i] * (Calmodulin_T - cab1[i]) + K_Calmodulin_MINUS * cab1[i]
        j_2[i] = -K_TroponinC_PLUS * f[i] * (TroponinC_T - cab2[i]) + K_TroponinC_MINUS * cab2[i]
        j_3[i] = -K_SR_PLUS * f[i] * (SR_T - cab3[i]) + K_SR_MINUS * cab3[i]
        j_4[i] = -K_SL_PLUS * f[i] * (SL_T - cab4[i]) + K_SL_MINUS * cab4[i]

    for i in range(0, point_count):
        new_cab1[i] = (-j_1[i]) * DT + cab1[i]
        new_cab2[i] = (-j_2[i]) * DT + cab2[i]
        new_cab3[i] = (-j_3[i]) * DT + cab3[i]
        new_cab4[i] = (-j_4[i]) * DT + cab4[i]
    return j_fdye, j_1, j_2, j_3, j_4, new_cab1, new_cab2, new_cab3, new_cab4


def open_calculation_f(f, caf, cab1, cab2, cab3, cab4):
    #  计算缓冲物和染料
    j_fdye, j_1, j_2, j_3, j_4, new_cab1, new_cab2, new_cab3, new_cab4 = cal_dye_and_buffers(f, caf, cab1, cab2, cab3,
                                                                                             cab4)
    grid_coordinates = np.loadtxt(open_grid_file_name, delimiter=",")  # 点坐标
    neighbors = np.loadtxt(neighbors_file_name, int, delimiter=",")  # 邻点
    # coefficients = np.loadtxt(coefficient_file_name, delimiter=",")*2  # 系数
    # coefficients = np.load("../config/open/coefficient.npy") * 2
    coefficients = np.load("../config/open/coefficient.npy")
    point_count = len(grid_coordinates)
    last = np.copy(f)
    # last = np.zeros(point_count)
    temp = np.copy(f)
    for i in range(0, 10):
        # last[:3] = f[:3]
        # 不计算边界点
        for j in range(3, point_count):
            j_buffers = j_1[j] + j_2[j] + j_3[j] + j_4[j]
            # 此处没有Jout，因为交界处在纳米空间中计算过了，开放空间只是用它的值，无需再计算
            p32, p12, p23, p21 = neighbors[j][0], neighbors[j][1], neighbors[j][2], neighbors[j][3]
            coe1, coe2, coe3, coe4 = coefficients[j][0], coefficients[j][1], coefficients[j][2], coefficients[j][3]
            # 学长没有取平均值
            f32, f12, f23, f21 = neighbors_concentration(j, temp, p32, p12, p23, p21)
            f32_n, f12_n, f23_n, f21_n = neighbors_concentration(j, f, p32, p12, p23, p21)
            if i != 0:
                f32, f12, f23, f21 = (f32 + f32_n) / 2.0, (f12 + f12_n) / 2.0, (f23 + f23_n) / 2.0, (f21 + f21_n) / 2.0
            item1 = D_CA * (f32 * coe1 + f12 * coe2 + f23 * coe3 + f21 * coe4)
            item3 = D_CA * (coe1 + coe2 + coe3 + coe4)
            last[j] = ((item1 + j_fdye[j] + j_buffers) * DT + f[j]) / (1 + item3 * DT)
        temp = np.copy(last)
    return temp, new_cab1, new_cab2, new_cab3, new_cab4


def open_calculation_caf(f, caf):
    grid_coordinates = np.loadtxt(open_grid_file_name, delimiter=",")  # 点坐标
    neighbors = np.loadtxt(neighbors_file_name, int, delimiter=",")  # 邻点
    # coefficients = np.loadtxt(coefficient_file_name, delimiter=",")*2  # 系数
    # coefficients = np.load("../config/open/coefficient.npy") * 2
    coefficients = np.load("../config/open/coefficient.npy")
    point_count = len(grid_coordinates)
    last = np.copy(caf)
    # last = np.zeros(point_count)
    temp = np.copy(caf)
    j_fdye = np.zeros(point_count)
    for j in range(0, point_count):
        j_fdye[j] = -K_F3_PLUS * f[j] * (F3_T - caf[j]) + K_F3_MINUS * caf[j]
    for i in range(0, 10):
        # last[:3]=caf[:3]
        for j in range(3, point_count):
            p32, p12, p23, p21 = neighbors[j][0], neighbors[j][1], neighbors[j][2], neighbors[j][3]
            coe1, coe2, coe3, coe4 = coefficients[j][0], coefficients[j][1], coefficients[j][2], coefficients[j][3]
            # 学长没有取平均值
            g32, g12, g23, g21 = neighbors_concentration(j, temp, p32, p12, p23, p21)
            g32_n, g12_n, g23_n, g21_n = neighbors_concentration(j, caf, p32, p12, p23, p21)
            if i != 0:
                g32, g12, g23, g21 = (g32 + g32_n) / 2, (g12 + g12_n) / 2, (g23 + g23_n) / 2, (g21 + g21_n) / 2
            item4 = D_F * (g32 * coe1 + g12 * coe2 + g23 * coe3 + g21 * coe4)
            item5 = D_F * (coe1 + coe2 + coe3 + coe4)
            last[j] = ((item4 - j_fdye[j]) * DT + caf[j]) / (1 + item5 * DT)
        temp = np.copy(last)
    return temp
