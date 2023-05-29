import os
from cal_bcnl import *
from nano_spark.nano_parameters import *
import numpy as np


# 计算各个点平均值
def cal_avg():
    grid_file_name = "../config/nano/4RYRgridt.dat"
    nod_file_name = "../config/nano/4RYRnod.dat"
    single_area, control_area, near_triangle, index_in_triangle, nix_multiply_l, niy_multiply_l, b_arr, c_arr, nmax, total_area = cal_elements(
        grid_file_name, nod_file_name)
    path = "../result/"
    dir_nano = "NANO/"
    dir_avg = "AVG/"
    dirnames = ["Ca", "CaF", "CaG"]
    for dirname in dirnames:
        # 初始化各点Ca浓度
        filenames = os.listdir(f"{path}{dir_nano}{dirname}/")  # 目前所有步数的浓度文件
        current_step = len(filenames)  # 当前生成文件数
        c_avg = np.zeros(current_step)
        i = 0
        for filename in filenames:
            nano_c = np.loadtxt(f"{path}{dir_nano}{dirname}/{filename}")
            total = 0.0
            NP = len(nano_c)
            for j in range(0, NP):
                total = total + nano_c[j] * control_area[j]
            total = total / 3.0
            avg = total / total_area
            c_avg[i] = avg
            i = i + 1
        np.savetxt(f"{path}{dir_avg}{dirname}_avg.csv", c_avg)


def main():
    cal_avg()


if __name__ == "__main__":
    main()
