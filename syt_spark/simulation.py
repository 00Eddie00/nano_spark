import numpy as np
from parameters import *
from tool_mkdir import *


def cal_r(point_index, two_neighbor):
    r1, r2 = -1, -1
    for i in range(len(two_neighbor)):
        if point_index == two_neighbor[i][0]:
            r1 = two_neighbor[i][1]
            r2 = two_neighbor[i][2]
            return r1, r2
    return r1, r2


def cal_f(f, g, h_1, h_2, h_3, h_4):
    neighbor = np.loadtxt("config/grid/neighbor.csv", int, delimiter=",")  # 邻点
    two_neighbor = np.loadtxt("config/grid/two_neighbor.csv", int, delimiter=",")  # 交界处的特殊邻点
    coefficient = np.load("config/coefficient/coefficient.npy")
    grid_coordinates = np.loadtxt("config/grid/grid_coordinates.csv", delimiter=",")  # 点坐标
    point_nums = len(grid_coordinates)
    last = np.copy(f)
    temp = np.copy(last)
    for k in range(10):
        j = 23
        for i in range(point_nums):
            J_h1 = -K_Calmodulin_PLUS * f[i] * (Calmodulin_T - h_1[i]) + K_Calmodulin_MINUS * h_1[i]
            J_h2 = -K_TroponinC_PLUS * f[i] * (TroponinC_T - h_2[i]) + K_TroponinC_MINUS * h_2[i]
            J_h3 = -K_SR_PLUS * f[i] * (SR_T - h_3[i]) + K_SR_MINUS * h_3[i]
            J_h4 = -K_SL_PLUS * f[i] * (SL_T - h_4[i]) + K_SL_MINUS * h_4[i]
            J_buffers = J_h1 + J_h2 + J_h3 + J_h4
            J_dye = 0.0
            if 0 <= i < 3360 or 3425 <= i < 3450:
                J_dye = -K_F_PLUS * f[i] * (F_T - g[i]) + K_F_MINUS * g[i]
            k1, k2 = 0.0, 0.0
            #  在入流边界
            if i == j and i != 287:
                j = j + 24
                k1 = DT * CA_JSR * K_RYR / V
                k2 = DT * K_RYR / V
            p3, p7, p1, p5 = neighbor[i]  # 上下外内
            f32, f12, f23, f21 = 0.0, 0.0, 0.0, 0.0
            if p3 == -1:  # 上面没有
                f32 = last[i]
            else:
                f32 = last[p3]
            if p7 == -1:  # 下面没有
                f12 = last[i]
            else:
                f12 = last[p7]
            if p1 == -1:  # 外面没有
                f23 = last[i]
            elif p1 == -2:  # 外面没有对应的
                r1, r2 = cal_r(i, two_neighbor)
                dist_1 = abs(grid_coordinates[r2][0] - grid_coordinates[i][0])
                dist_2 = abs(grid_coordinates[r1][0] - grid_coordinates[i][0])
                dist = dist_1 + dist_2
                f23 = (last[r1] * dist_1 + last[r2] * dist_2) / dist
            else:
                f23 = last[p1]
            if p5 == -1:  # 里面没有
                f21 = last[i]
            elif p5 == -2:  # 里面没有对应的
                r1, r2 = cal_r(i, two_neighbor)
                dist_1 = abs(grid_coordinates[r2][0] - grid_coordinates[i][0])
                dist_2 = abs(grid_coordinates[r1][0] - grid_coordinates[i][0])
                dist = dist_1 + dist_2
                f21 = (last[r1] * dist_1 + last[r2] * dist_2) / dist
            else:
                f21 = last[p5]

            diffusion_numerator = D_CA * (
                        f32 * coefficient[i][0] - f12 * coefficient[i][1] + f23 * coefficient[i][2] - f21 *
                        coefficient[i][3]) * DT

            diffusion_denominator = D_CA * (
                        -coefficient[i][0] + coefficient[i][1] - coefficient[i][2] + coefficient[i][3]) * DT
            temp[i] = (diffusion_numerator + (J_dye + J_buffers) * DT + k1 + f[i]) / (1 - diffusion_denominator + k2)

        last = np.copy(temp)
    return last


def cal_g_h(f, g, h_1, h_2, h_3, h_4, point_nums):
    h_1_new = np.empty(point_nums)
    h_2_new = np.empty(point_nums)
    h_3_new = np.empty(point_nums)
    h_4_new = np.empty(point_nums)
    g_1_new = np.empty(point_nums)

    for i in range(point_nums):
        J_h1 = -K_Calmodulin_PLUS * f[i] * (Calmodulin_T - h_1[i]) + K_Calmodulin_MINUS * h_1[i]
        J_h2 = -K_TroponinC_PLUS * f[i] * (TroponinC_T - h_2[i]) + K_TroponinC_MINUS * h_2[i]
        J_h3 = -K_SR_PLUS * f[i] * (SR_T - h_3[i]) + K_SR_MINUS * h_3[i]
        J_h4 = -K_SL_PLUS * f[i] * (SL_T - h_4[i]) + K_SL_MINUS * h_4[i]
        J_dye = 0.0
        if 0 <= i < 3360 or 3425 <= i < 3450:
            J_dye = -K_F_PLUS * f[i] * (F_T - g[i]) + K_F_MINUS * g[i]

        g_1_new[i] = DT * (-J_dye) + g[i]
        h_1_new[i] = DT * (-J_h1) + h_1[i]
        h_2_new[i] = DT * (-J_h2) + h_2[i]
        h_3_new[i] = DT * (-J_h3) + h_3[i]
        h_4_new[i] = DT * (-J_h4) + h_4[i]
    return g_1_new, h_1_new, h_2_new, h_3_new, h_4_new


def simulation(dir_name, total_steps):
    # 创建用于存储结果的文件夹
    mkdir(dir_name + "Ca")
    mkdir(dir_name + "CaF")
    mkdir(dir_name + "CaB")
    #  初始化各物质浓度
    grid_coordinates = np.loadtxt("config/grid/grid_coordinates.csv", delimiter=",")
    point_nums = len(grid_coordinates)
    f = np.full(point_nums, INITIAL_C_CA)
    INITIAL_G = K_F_PLUS * INITIAL_C_CA * F_T / (K_F_PLUS * INITIAL_C_CA + K_F_MINUS)
    g = np.empty(point_nums)
    g[:3360] = INITIAL_G
    g[3425:3450] = INITIAL_G
    INITIAL_H1 = K_Calmodulin_PLUS * INITIAL_C_CA * Calmodulin_T / (
            K_Calmodulin_PLUS * INITIAL_C_CA + K_Calmodulin_MINUS)
    h1 = np.full(point_nums, INITIAL_H1)
    INITIAL_H2 = K_TroponinC_PLUS * INITIAL_C_CA * TroponinC_T / (K_TroponinC_PLUS * INITIAL_C_CA + K_TroponinC_MINUS)
    h2 = np.full(point_nums, INITIAL_H2)
    INITIAL_H3 = K_SR_PLUS * INITIAL_C_CA * SR_T / (K_SR_PLUS * INITIAL_C_CA + K_SR_MINUS)
    h3 = np.full(point_nums, INITIAL_H3)
    INITIAL_H4 = K_SL_PLUS * INITIAL_C_CA * SL_T / (K_SL_PLUS * INITIAL_C_CA + K_SL_MINUS)
    h4 = np.full(point_nums, INITIAL_H4)
    current_step = 1
    # 记录初始时刻钙离子浓度值
    ca_file_name = f"{dir_name}Ca/Ca00000000.csv"
    np.savetxt(ca_file_name, f)
    print(ca_file_name, "SAVED")
    # 记录初始时刻荧光钙离子浓度值
    caf_file_name = f"{dir_name}CaF/CaF00000000.csv"
    np.savetxt(caf_file_name, g)
    print(caf_file_name, "SAVED")

    for i in range(current_step, total_steps + 1):
        last = cal_f(f, g, h1, h2, h3, h4)
        g_1_new, h_1_new, h_2_new, h_3_new, h_4_new = cal_g_h(f, g, h1, h2, h3, h4, point_nums)
        length = len(str(i))
        rest_name = "0" * (8 - length) + str(i) + ".dat"
        ca_file_name = f"{dir_name}Ca/Ca{rest_name}.csv"
        np.savetxt(ca_file_name, last)
        f = np.copy(last)
        print(ca_file_name, "SAVED")
        caf_file_name = f"{dir_name}CaF/CaF{rest_name}.csv"
        np.savetxt(caf_file_name, g_1_new)
        g = np.copy(g_1_new)
        print(caf_file_name, "SAVED")

        # 记录各缓冲物浓度
        cab_file_name = f"{dir_name}CaB/CaB{rest_name}.csv"
        np.savetxt(cab_file_name,
                   np.concatenate(
                       (h_1_new[:, np.newaxis], h_2_new[:, np.newaxis], h_3_new[:, np.newaxis], h_4_new[:, np.newaxis]),
                       axis=1), delimiter=',')
        h1 = np.copy(h_1_new)
        h2 = np.copy(h_2_new)
        h3 = np.copy(h_3_new)
        h4 = np.copy(h_4_new)
        print(cab_file_name, "SAVED")
        # 记录已经计算的步数
        with open(dir_name + "STATE.csv", "w") as state_file:
            state_file.write(str(i))


if __name__ == '__main__':
    simulation("", 100)
