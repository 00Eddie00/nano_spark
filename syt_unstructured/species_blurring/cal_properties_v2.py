import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import interp1d

# 使用插值函数，使画出的FDHM更直
def cal_temporal_properties(t, ampl, time_interval):
    peak_value = np.max(ampl)  # 峰值，即某空间点在时间尺度上的最大荧光强度
    base_value = ampl[0]
    t_rise_index = np.argmax(ampl)  # 峰化时间对应ampl中的索引
    half_value = (peak_value + base_value) / 2
    # 使用线性插值来估计缺失点的 X 值
    interpolator1 = interp1d(ampl[:t_rise_index], t[:t_rise_index], kind='linear')
    # 计算缺失点的 X 值
    missing_x1 = interpolator1(half_value)
    # 使用线性插值来估计缺失点的 X 值
    interpolator2 = interp1d(ampl[t_rise_index:], t[t_rise_index:], kind='linear')
    # 计算缺失点的 X 值
    missing_x2 = interpolator2(half_value)
    t_rise = t_rise_index * time_interval  # 峰化时间，即荧光强度从基线达到峰值所需时间
    t50 = missing_x2  # 荧光强度从峰值衰减50%所需时间
    fdhm = missing_x2 - missing_x1  # full duration at half maximum，半峰全持续时间
    fdhm_start = missing_x1
    return peak_value, t_rise, t50, fdhm, fdhm_start, half_value


def temporal_plotter(temporal_path, time_interval, save=False):
    temporal_fluorescence = np.loadtxt(temporal_path)
    # ampl = temporal_fluorescence / temporal_fluorescence[0] - 1  # 振幅，即荧光强度：fluorescence intensity
    ampl = temporal_fluorescence / temporal_fluorescence[0]  # 振幅，即荧光强度（直接用浓度）：fluorescence intensity
    t = [i * time_interval for i in range(len(temporal_fluorescence))]  # 时间间隔，单位ms
    peak_value, t_rise, t50, fdhm, fdhm_start, half_value = cal_temporal_properties(t, ampl, time_interval)  # 计算钙火花的各属性
    plt.figure(dpi=300, figsize=(10, 10))
    plt.plot(t, ampl, ls='-', lw=2)
    # 绘制tRise时间点及峰值
    plt.plot(t_rise, peak_value, 'o', color='orange')
    plt.text(t_rise, peak_value, f'tRise={np.round(t_rise, 4)}ms peak={np.round(peak_value, 4)}', fontsize=14)
    # 绘制t50时间点
    t50_t_index = t50
    t50_ampl = half_value
    plt.plot(t50_t_index, t50_ampl, 'o', color='orange')
    plt.text(t50_t_index, t50_ampl, f't50={np.round(t50 - t_rise, 4)}ms', fontsize=14)
    # 绘制FDHM
    fdhm_start_ampl = half_value
    plt.plot(fdhm_start, fdhm_start_ampl, 'o', color='orange')
    plt.plot([fdhm_start, t50_t_index], [fdhm_start_ampl, t50_ampl], ls='--', color='orange')
    plt.text(t_rise, fdhm_start_ampl, f'FDHM={np.round(fdhm, 4)}ms', fontsize=12)
    plt.xlabel('t(ms)', fontsize=14)
    # plt.ylabel('ΔF/F0', fontsize=14)
    # 直接用浓度
    plt.ylabel('F/F0', fontsize=14)
    plt.suptitle('中心：时间分布')
    file_without_extension = os.path.splitext(temporal_path)[0]
    plt.title(f"数据来源：{file_without_extension}")
    plt.grid()
    if save:
        plt.savefig(f"{file_without_extension}.png")
    # plt.show()


def main():
    dir_path = "../result/"
    time_interval = 2 * 10 ** -6 * 100 * 1000

    # # caf卷积画图
    # dis_dirname = "NANO"
    # species = "CaF"
    # caf_position_list = [[0, 0], [100, 0], [300, 0], [400, 0]]
    # for i in range(len(caf_position_list)):
    #     dis_suffix = f"_psf_({caf_position_list[i][0]},{caf_position_list[i][1]}).csv"
    #     temporal_path = f"{dir_path}{dis_dirname}/{species}{dis_suffix}"
    #     # print(temporal_path)
    #     temporal_plotter(temporal_path, time_interval, save=True)

    # cag卷积画图
    dis_dirname = "NANO_basic_parameter"
    species = "CaG"
    cag_position_list = [[0, 0], [300, 0]]
    for i in range(len(cag_position_list)):
        dis_suffix = f"_psf_({cag_position_list[i][0]},{cag_position_list[i][1]}).csv"
        temporal_path = f"{dir_path}{dis_dirname}/{species}{dis_suffix}"
        temporal_plotter(temporal_path, time_interval, save=True)

    # 纳米空间平均值图
    # dis_dirname = "NANO_AVG"
    # avg_suffix = "_avg.csv"
    # species = "Ca"
    # temporal_path = f"{dir_path}{dis_dirname}/{species}{avg_suffix}"
    # temporal_plotter(temporal_path, time_interval, save=True)

    # 开放空间平均值图
    # dis_dirname = "OPEN_AVG"
    # avg_prefix = "avg_c_"
    # species = "CaF"
    # temporal_path = f"{dir_path}{dis_dirname}/{avg_prefix}{species}.csv"
    # temporal_plotter(temporal_path, time_interval, save=True)


if __name__ == "__main__":
    matplotlib.use('TkAgg')
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
    main()
