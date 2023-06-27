from cal_nano_element import *
from cal_open_elements import *
from nano_parameters import *
from open_parameters import open_grid_file_name
from tool.tool_mkdir import *
import os


def nano_spark(is_continue, total_steps):
    # ryr通道开放时间对应的步数
    release_step = int(RELEASE_TIME / DT)
    grid_coordinates = np.loadtxt(open_grid_file_name, delimiter=",")  # 点坐标
    out_boundary_length = cal_out_boundary_length(grids)
    # 点个数
    point_count = len(grid_coordinates)
    # 纳米入流处仍然固定
    k_ryr = K_RYR
    ca_jsr = CA_JSR

    # 文件路径前缀
    dir_name = "../result"
    nano_c_ca_prefix = f"{dir_name}/NANO/Ca/Ca"
    nano_c_cag_prefix = f"{dir_name}/NANO/CaG/CaG"
    nano_c_caf_prefix = f"{dir_name}/NANO/CaF/CaF"
    nano_c_cab_prefix = f"{dir_name}/NANO/CaB/CaB"
    open_c_ca_prefix = f"{dir_name}/OPEN/Ca/Ca"
    open_c_caf_prefix = f"{dir_name}/OPEN/CaF/CaF"
    open_c_cab_prefix = f"{dir_name}/OPEN/CaB/CaB"
    # 从0开始
    current_step = 1
    if is_continue:
        path = "../result/"
        dir_nano = "NANO/"
        dir_open = "OPEN/"
        dirnames = ["Ca", "CaB", "CaF", "CaG"]
        # 初始化各点Ca浓度
        filenames = os.listdir(f"{path}{dir_nano}{dirnames[0]}/")  # 目前所有步数的浓度文件
        current_step = len(filenames)  # 当前生成文件数
        last_filename = filenames[-1]
        nano_f = np.loadtxt(f"{path}{dir_nano}{dirnames[0]}/{last_filename}")
        open_f = np.loadtxt(f"{path}{dir_open}{dirnames[0]}/{last_filename}")

        # 初始化各点缓冲物浓度
        filenames = os.listdir(f"{path}{dir_nano}{dirnames[1]}/")  # 目前所有步数的浓度文件
        last_filename = filenames[-1]
        nano_cab = np.loadtxt(f"{path}{dir_nano}{dirnames[1]}/{last_filename}")
        open_cab = np.loadtxt(f"{path}{dir_open}{dirnames[1]}/{last_filename}")
        nano_cab1 = nano_cab[:, 0]
        open_cab1 = open_cab[:, 0]
        nano_cab2 = nano_cab[:, 1]
        open_cab2 = open_cab[:, 1]
        nano_cab3 = nano_cab[:, 2]
        open_cab3 = open_cab[:, 2]
        nano_cab4 = nano_cab[:, 3]
        open_cab4 = open_cab[:, 3]

        # 初始化各点CaF浓度
        filenames = os.listdir(f"{path}{dir_nano}{dirnames[2]}/")  # 目前所有步数的浓度文件
        last_filename = filenames[-1]
        nano_caf = np.loadtxt(f"{path}{dir_nano}{dirnames[2]}/{last_filename}")
        open_caf = np.loadtxt(f"{path}{dir_open}{dirnames[2]}/{last_filename}")

        # 初始化各点CaG浓度
        filenames = os.listdir(f"{path}{dir_nano}{dirnames[3]}/")  # 目前所有步数的浓度文件
        last_filename = filenames[-1]
        nano_cag = np.loadtxt(f"{path}{dir_nano}{dirnames[3]}/{last_filename}")
        # 开放空间中无该染料
    else:
        # 初始化各点Ca浓度
        nano_f = np.full(NP, INITIAL_C_CA)
        open_f = np.full(point_count, INITIAL_C_CA)

        # 初始化各点CaF浓度
        INITIAL_F3 = K_F3_PLUS * INITIAL_C_CA * F3_T / (K_F3_PLUS * INITIAL_C_CA + K_F3_MINUS)
        nano_caf = np.full(NP, INITIAL_F3)
        open_caf = np.full(point_count, INITIAL_F3)

        # 初始化各点CaG浓度
        INITIAL_G = K_GCaMP6f_PLUS * INITIAL_C_CA * GCaMP6f_T / (K_GCaMP6f_PLUS * INITIAL_C_CA + K_GCaMP6f_MINUS)
        nano_cag = np.full(NP, INITIAL_G)
        # 开放空间中无该染料

        # 初始化各点Calmodulin浓度
        INITIAL_H1 = K_Calmodulin_PLUS * INITIAL_C_CA * Calmodulin_T / (
                K_Calmodulin_PLUS * INITIAL_C_CA + K_Calmodulin_MINUS)
        nano_cab1 = np.full(NP, INITIAL_H1)
        open_cab1 = np.full(point_count, INITIAL_H1)

        # 初始化各点Troponin C浓度
        INITIAL_H2 = K_TroponinC_PLUS * INITIAL_C_CA * TroponinC_T / (
                K_TroponinC_PLUS * INITIAL_C_CA + K_TroponinC_MINUS)
        nano_cab2 = np.full(NP, INITIAL_H2)
        open_cab2 = np.full(point_count, INITIAL_H2)

        # 初始化各点SR浓度
        INITIAL_H3 = K_SR_PLUS * INITIAL_C_CA * SR_T / (K_SR_PLUS * INITIAL_C_CA + K_SR_MINUS)
        nano_cab3 = np.full(NP, INITIAL_H3)
        open_cab3 = np.full(point_count, INITIAL_H3)

        # 初始化各点SL浓度
        INITIAL_H4 = K_SL_PLUS * INITIAL_C_CA * SL_T / (K_SL_PLUS * INITIAL_C_CA + K_SL_MINUS)
        nano_cab4 = np.full(NP, INITIAL_H4)
        open_cab4 = np.full(point_count, INITIAL_H4)

        # 生成保存纳米空间的文件
        mkdir(f"{dir_name}/NANO/Ca")
        mkdir(f"{dir_name}/NANO/CaG")
        mkdir(f"{dir_name}/NANO/CaF")
        mkdir(f"{dir_name}/NANO/CaB")
        mkdir(f"{dir_name}/NANO/nano_avg")

        # 生成保存开放空间的文件
        mkdir(f"{dir_name}/OPEN/Ca")
        mkdir(f"{dir_name}/OPEN/CaF")
        mkdir(f"{dir_name}/OPEN/CaB")

        # 纳米空间
        # 记录初始时刻纳米空间钙离子浓度值
        ca_file_name = f"{nano_c_ca_prefix}00000000.csv"
        np.savetxt(ca_file_name, nano_f)
        print(ca_file_name, "SAVED")

        # 记录初始时刻纳米空间荧光（GCaMP6f）钙离子浓度值
        cag_file_name = f"{nano_c_cag_prefix}00000000.csv"
        np.savetxt(cag_file_name, nano_cag)
        print(cag_file_name, "SAVED")

        # 记录初始时刻纳米空间荧光（Fluo-3）钙离子浓度值
        caf_file_name = f"{nano_c_caf_prefix}00000000.csv"
        np.savetxt(caf_file_name, nano_caf)
        print(caf_file_name, "SAVED")

        # 记录初始时刻纳米空间缓冲物浓度值
        cab_file_name = f"{nano_c_cab_prefix}00000000.csv"
        np.savetxt(cab_file_name, np.stack((nano_cab1, nano_cab2, nano_cab3, nano_cab4), 1))
        print(cab_file_name, "SAVED")

        # 开放空间
        # 记录初始时刻开放空间钙离子浓度值
        ca_file_name = f"{open_c_ca_prefix}00000000.csv"
        np.savetxt(ca_file_name, open_f)
        print(ca_file_name, "SAVED")

        # 记录初始时刻开放空间荧光（Fluo-3）钙离子浓度值
        caf_file_name = f"{open_c_caf_prefix}00000000.csv"
        np.savetxt(caf_file_name, open_caf)
        print(caf_file_name, "SAVED")

        # 记录初始时刻开放空间缓冲物浓度值
        cab_file_name = f"{open_c_cab_prefix}00000000.csv"
        np.savetxt(cab_file_name, np.stack((open_cab1, open_cab2, open_cab3, open_cab4), 1))
        print(cab_file_name, "SAVED")
    # 开始计算
    for i in range(current_step, total_steps + 1):
        # 计算纳米空间
        c_ca_out = open_f[69]
        c_caf_out = open_caf[69]
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

        # n + 1 步的浓度，计算出流边界平均值
        out_boundary_c_ca = np.sum(new_nano_f[84:205] * out_boundary_length) / np.sum(out_boundary_length)
        out_boundary_c_caf = np.sum(new_nano_caf[84:205] * out_boundary_length) / np.sum(out_boundary_length)
        new_nano_f[84:205] = out_boundary_c_ca
        new_nano_caf[84:205] = out_boundary_c_caf

        # 计算开放空间
        # 因此开放空间不计算交界处
        # n 步的浓度，计算出流边界平均值
        out_boundary_c_ca = np.sum(nano_f[84:205] * out_boundary_length) / np.sum(out_boundary_length)
        out_boundary_c_caf = np.sum(nano_caf[84:205] * out_boundary_length) / np.sum(out_boundary_length)
        open_f[:3] = out_boundary_c_ca
        open_caf[:3] = out_boundary_c_caf
        new_open_f, new_open_cab1, new_open_cab2, new_open_cab3, new_open_cab4 = open_calculation_f(open_f, open_caf,
                                                                                                    open_cab1,
                                                                                                    open_cab2,
                                                                                                    open_cab3,
                                                                                                    open_cab4)
        new_open_caf = open_calculation_caf(open_f, open_caf)

        # 保存文件
        length = len(str(i))
        rest_name = "0" * (8 - length) + str(i)

        # 纳米空间
        # 记录该时刻钙离子浓度值
        ca_file_name = f"{nano_c_ca_prefix}{rest_name}.csv"
        np.savetxt(ca_file_name, new_nano_f)
        nano_f = np.copy(new_nano_f)
        print(ca_file_name, "SAVED")
        # 记录该时刻荧光（GCaMP6f）钙离子浓度值
        cag_file_name = f"{nano_c_cag_prefix}{rest_name}.csv"
        np.savetxt(cag_file_name, new_nano_cag)
        nano_cag = np.copy(new_nano_cag)
        print(cag_file_name, "SAVED")
        # 记录该时刻荧光（Fluo-3）钙离子浓度值
        caf_file_name = f"{nano_c_caf_prefix}{rest_name}.csv"
        np.savetxt(caf_file_name, new_nano_caf)
        nano_caf = np.copy(new_nano_caf)
        print(caf_file_name, "SAVED")
        # 记录该时刻缓冲物浓度值
        cab_file_name = f"{nano_c_cab_prefix}{rest_name}.csv"
        np.savetxt(cab_file_name, np.stack((new_nano_cab1, new_nano_cab2, new_nano_cab3, new_nano_cab4), 1))
        nano_cab1 = np.copy(new_nano_cab1)
        nano_cab2 = np.copy(new_nano_cab2)
        nano_cab3 = np.copy(new_nano_cab3)
        nano_cab4 = np.copy(new_nano_cab4)
        print(cab_file_name, "SAVED")

        # 开放空间
        # 记录初始时刻开放空间钙离子浓度值
        open_ca_file_name = f"{open_c_ca_prefix}{rest_name}.csv"
        np.savetxt(open_ca_file_name, new_open_f)
        open_f = np.copy(new_open_f)
        print(open_ca_file_name, "SAVED")

        # 记录初始时刻开放空间荧光（Fluo-3）钙离子浓度值
        open_caf_file_name = f"{open_c_caf_prefix}{rest_name}.csv"
        np.savetxt(open_caf_file_name, new_open_caf)
        open_caf = np.copy(new_open_caf)
        print(open_caf_file_name, "SAVED")

        # 记录初始时刻开放空间缓冲物浓度值
        open_cab_file_name = f"{open_c_cab_prefix}{rest_name}.csv"
        np.savetxt(open_cab_file_name, np.stack((new_open_cab1, new_open_cab2, new_open_cab3, new_open_cab4), 1))
        open_cab1 = np.copy(new_open_cab1)
        open_cab2 = np.copy(new_open_cab2)
        open_cab3 = np.copy(new_open_cab3)
        open_cab4 = np.copy(new_open_cab4)
        print(open_cab_file_name, "SAVED")

        # 控制ryr通道关闭
        if i == release_step:
            k_ryr = 0


if __name__ == '__main__':
    is_continue = True
    total_steps = 20000
    nano_spark(is_continue, total_steps)
