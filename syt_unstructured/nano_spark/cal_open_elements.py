from open_parameters import *
from nano_parameters import *
import numpy as np


def cal_dye_and_buffers(f, caf, cab1, cab2, cab3, cab4):
    grid_coordinates = np.loadtxt(grid_file_name, delimiter=",")  # 点坐标
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
    for i in range(0, NP):
        j_fdye[i] = -K_F3_PLUS * f[i] * (F3_T - caf[i]) + K_F3_MINUS * caf[i]
        j_1[i] = -K_Calmodulin_PLUS * f[i] * (Calmodulin_T - cab1[i]) + K_Calmodulin_MINUS * cab1[i]
        j_2[i] = -K_TroponinC_PLUS * f[i] * (TroponinC_T - cab2[i]) + K_TroponinC_MINUS * cab2[i]
        j_3[i] = -K_SR_PLUS * f[i] * (SR_T - cab3[i]) + K_SR_MINUS * cab3[i]
        j_4[i] = -K_SL_PLUS * f[i] * (SL_T - cab4[i]) + K_SL_MINUS * cab4[i]

    for i in range(0, NP):
        new_cab1[i] = (-j_1[i]) * DT + cab1[i]
        new_cab2[i] = (-j_2[i]) * DT + cab2[i]
        new_cab3[i] = (-j_3[i]) * DT + cab3[i]
        new_cab4[i] = (-j_4[i]) * DT + cab4[i]
    return j_fdye, j_1, j_2, j_3, j_4, new_cab1, new_cab2, new_cab3, new_cab4


def open_calculation_f(f, caf, cab1, cab2, cab3, cab4):
    j_fdye, j_1, j_2, j_3, j_4, new_cab1, new_cab2, new_cab3, new_cab4 = cal_dye_and_buffers(f, caf, cab1, cab2, cab3,
                                                                                             cab4)
    grid_coordinates = np.loadtxt(grid_file_name, delimiter=",")  # 点坐标
    neighbors = np.loadtxt(neighbors_file_name, int, delimiter=",")  # 邻点
    coefficient = np.loadtxt(coefficient_file_name, delimiter=",")  # 系数
    point_count = len(grid_coordinates)
    last = np.zeros(point_count)
    temp = np.zeros(point_count)
    for i in range(0, point_count):
        j_buffers = j_1[i] + j_2[i] + j_3[i] + j_4[i]
        item1, item2, item3 = 0, 0, 0


