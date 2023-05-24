from cal_nano_element import *
from nano_parameters import *
from syt_unstructured.tool.tool_mkdir import mkdir


def nano_spark():
    # ryr通道开放时间对应的步数
    release_step = int(RELEASE_TIME / DT)

    # 初始化各点Ca浓度
    nano_f = np.full(NP, INITIAL_C_CA)
    open_f = np.full(NP, INITIAL_C_CA)

    # 初始化各点CaF浓度
    INITIAL_F3 = K_F3_PLUS * INITIAL_C_CA * F3_T / (K_F3_PLUS * INITIAL_C_CA + K_F3_MINUS)
    nano_caf = np.full(NP, INITIAL_F3)
    open_caf = np.full(NP, INITIAL_F3)

    # 初始化各点CaG浓度
    INITIAL_G = K_GCaMP6f_PLUS * INITIAL_C_CA * GCaMP6f_T / (K_GCaMP6f_PLUS * INITIAL_C_CA + K_GCaMP6f_MINUS)
    nano_cag = np.full(NP, INITIAL_G)

    # 初始化各点Calmodulin浓度
    INITIAL_H1 = K_Calmodulin_PLUS * INITIAL_C_CA * Calmodulin_T / (
            K_Calmodulin_PLUS * INITIAL_C_CA + K_Calmodulin_MINUS)
    nano_cab1 = np.full(NP, INITIAL_H1)
    open_cab1 = np.full(NP, INITIAL_H1)

    # 初始化各点Troponin C浓度
    INITIAL_H2 = K_TroponinC_PLUS * INITIAL_C_CA * TroponinC_T / (K_TroponinC_PLUS * INITIAL_C_CA + K_TroponinC_MINUS)
    nano_cab2 = np.full(NP, INITIAL_H2)
    open_cab2 = np.full(NP, INITIAL_H2)

    # 初始化各点SR浓度
    INITIAL_H3 = K_SR_PLUS * INITIAL_C_CA * SR_T / (K_SR_PLUS * INITIAL_C_CA + K_SR_MINUS)
    nano_cab3 = np.full(NP, INITIAL_H3)
    open_cab3 = np.full(NP, INITIAL_H3)

    # 初始化各点SL浓度
    INITIAL_H4 = K_SL_PLUS * INITIAL_C_CA * SL_T / (K_SL_PLUS * INITIAL_C_CA + K_SL_MINUS)
    nano_cab4 = np.full(NP, INITIAL_H4)
    open_cab4 = np.full(NP, INITIAL_H4)

    out_boundary_length = cal_out_boundary_length(grid)

    # 开始计算
    total_steps = 2000
    avg_arr = np.empty((total_steps + 1, 7))
    dir_name = "../result"
    mkdir(f"{dir_name}/NANO/Ca")
    mkdir(f"{dir_name}/NANO/CaG")
    mkdir(f"{dir_name}/NANO/CaF")
    mkdir(f"{dir_name}/NANO/CaB")
    mkdir(f"{dir_name}/NANO/nano_avg")
    nano_c_ca_prefix = f"{dir_name}/NANO/Ca/Ca"
    nano_c_cag_prefix = f"{dir_name}/NANO/CaG/CaG"
    nano_c_caf_prefix = f"{dir_name}/NANO/CaF/CaF"
    nano_c_cab_prefix = f"{dir_name}/NANO/CaB/CaB"

    # 记录初始时刻钙离子浓度值
    ca_file_name = f"{nano_c_ca_prefix}00000000.csv"
    np.savetxt(ca_file_name, nano_f)
    # 记录初始时刻钙离子平均值
    avg_arr[0, 0] = cal_avg(nano_f, grid_file_name, nod_file_name)
    print(ca_file_name, "SAVED")

    # 记录初始时刻荧光（GCaMP6f）钙离子浓度值
    cag_file_name = f"{nano_c_cag_prefix}00000000.csv"
    np.savetxt(cag_file_name, nano_cag)
    # 记录初始时刻荧光（GCaMP6f）钙离子平均值
    avg_arr[0, 1] = cal_avg(nano_cag, grid_file_name, nod_file_name)
    print(cag_file_name, "SAVED")

    # 记录初始时刻荧光（Fluo-3）钙离子浓度值
    caf_file_name = f"{nano_c_caf_prefix}00000000.csv"
    np.savetxt(caf_file_name, nano_caf)
    # 记录初始时刻荧光（Fluo-3）钙离子平均值
    avg_arr[0, 2] = cal_avg(nano_caf, grid_file_name, nod_file_name)
    print(caf_file_name, "SAVED")

    # 记录初始时刻缓冲物浓度值
    cab_file_name = f"{nano_c_cab_prefix}00000000.csv"
    np.savetxt(cab_file_name, np.stack((nano_cab1, nano_cab2, nano_cab3, nano_cab4), 1))
    # 记录初始时刻缓冲物平均值
    avg_arr[0, 3] = cal_avg(nano_cab1, grid_file_name, nod_file_name)
    avg_arr[0, 4] = cal_avg(nano_cab2, grid_file_name, nod_file_name)
    avg_arr[0, 5] = cal_avg(nano_cab3, grid_file_name, nod_file_name)
    avg_arr[0, 6] = cal_avg(nano_cab4, grid_file_name, nod_file_name)
    print(cab_file_name, "SAVED")

    c_ca_out = open_f[69]
    c_caf_out = open_caf[69]
    k_ryr = K_RYR
    ca_jsr = CA_JSR
    for i in range(1, total_steps + 1):
        new_nano_f, new_nano_cag, new_nano_cab1, new_nano_cab2, new_nano_cab3, new_nano_cab4 = nano_calculation_f(k_ryr,
                                                                                                                  nano_f,
                                                                                                                  nano_caf,
                                                                                                                  nano_cag,
                                                                                                                  nano_cab1,
                                                                                                                  nano_cab2,
                                                                                                                  nano_cab3,
                                                                                                                  nano_cab4,
                                                                                                                  ca_jsr,
                                                                                                                  c_ca_out)
        new_nano_caf = nano_calculation_g2(nano_f, nano_caf, c_caf_out)

        # temp_out_boundary_f = 0.0
        # temp_out_boundary_caf = 0.0
        # total_out_boundary_length = 0.0
        # for j in range(0, len(out_boundary_length)):
        #     temp_out_boundary_f = temp_out_boundary_f + new_f[j + 84] * out_boundary_length[j]
        #     temp_out_boundary_caf = temp_out_boundary_caf + new_caf[j + 84] * out_boundary_length[j]
        #     total_out_boundary_length = total_out_boundary_length + out_boundary_length[j]
        # out_boundary_c_ca_n1 = temp_out_boundary_f / total_out_boundary_length
        # out_boundary_c_caf_n1 = temp_out_boundary_caf / total_out_boundary_length
        #
        # for j in range(84,205):
        #     new_f[j] = out_boundary_c_ca_n1
        #     new_caf[j] = out_boundary_c_caf_n1
        # n + 1 步的浓度，计算出流边界平均值
        out_boundary_c_ca = np.sum(new_nano_f[84:205] * out_boundary_length) / np.sum(out_boundary_length)
        out_boundary_c_caf = np.sum(new_nano_caf[84:205] * out_boundary_length) / np.sum(out_boundary_length)
        new_nano_f[84:205] = out_boundary_c_ca
        new_nano_caf[84:205] = out_boundary_c_caf

        length = len(str(i))
        rest_name = "0" * (8 - length) + str(i) + ".dat"
        # 记录该时刻钙离子浓度值
        ca_file_name = f"{nano_c_ca_prefix}{rest_name}.csv"
        np.savetxt(ca_file_name, new_nano_f)
        avg_arr[i, 0] = cal_avg(new_nano_f, grid_file_name, nod_file_name)
        nano_f = np.copy(new_nano_f)
        print(ca_file_name, "SAVED")
        # 记录该时刻荧光（GCaMP6f）钙离子浓度值
        cag_file_name = f"{nano_c_cag_prefix}{rest_name}.csv"
        np.savetxt(cag_file_name, new_nano_cag)
        avg_arr[i, 1] = cal_avg(new_nano_cag, grid_file_name, nod_file_name)
        nano_cag = np.copy(new_nano_cag)
        print(cag_file_name, "SAVED")
        # 记录该时刻荧光（Fluo-3）钙离子浓度值
        caf_file_name = f"{nano_c_caf_prefix}{rest_name}.csv"
        np.savetxt(caf_file_name, new_nano_caf)
        avg_arr[i, 2] = cal_avg(new_nano_caf, grid_file_name, nod_file_name)
        nano_caf = np.copy(new_nano_caf)
        print(caf_file_name, "SAVED")
        # 记录该时刻缓冲物浓度值
        cab_file_name = f"{nano_c_cab_prefix}{rest_name}.csv"
        np.savetxt(cab_file_name, np.stack((new_nano_cab1, new_nano_cab2, new_nano_cab3, new_nano_cab4), 1))
        avg_arr[i, 3] = cal_avg(new_nano_cab1, grid_file_name, nod_file_name)
        avg_arr[i, 4] = cal_avg(new_nano_cab2, grid_file_name, nod_file_name)
        avg_arr[i, 5] = cal_avg(new_nano_cab3, grid_file_name, nod_file_name)
        avg_arr[i, 6] = cal_avg(new_nano_cab4, grid_file_name, nod_file_name)
        nano_cab1 = np.copy(new_nano_cab1)
        nano_cab2 = np.copy(new_nano_cab2)
        nano_cab3 = np.copy(new_nano_cab3)
        nano_cab4 = np.copy(new_nano_cab4)
        print(cab_file_name, "SAVED")
        # 控制ryr通道关闭
        if i == release_step:
            k_ryr = 0
    np.savetxt("../result/NANO/nano_avg/avg.dat", avg_arr)
    print("avg", "SAVED")


if __name__ == '__main__':
    nano_spark()
