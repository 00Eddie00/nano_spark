import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mlp

plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号


#  画折线图
def plotter(xx, yy, name,dir_name, file_name):
    path = dir_name + file_name
    # plt.figure(dpi=300, figsize=(4, 4))  # 创建图形窗口
    plt.plot(xx, yy, ls='-', lw=2)

    plt.xlabel('steps', fontsize=14)  # 设置x轴标签
    plt.ylabel('concentration', fontsize=14)  # 设置y轴标签
    plt.tick_params(axis='x', labelrotation=20)  # 设置x轴刻度逆时针旋转20度
    plt.suptitle(name + '均值', fontsize=24)  # 设置折线图更上一级标题
    plt.title("数据来源：" + path)  # 设置折线图标题
    # 设置刻度间隔
    # x_major_locator = MultipleLocator(2500)  # 把x轴的刻度间隔设置为1000，并存在变量里
    # y_major_locator = MultipleLocator(y_locator)  # 把y轴的刻度间隔设置为0.02，并存在变量里
    # ax = plt.gca()  # ax为两条坐标轴的实例
    # ax.xaxis.set_major_locator(x_major_locator)  # 把x轴的主刻度设置为1000的倍数
    # ax.yaxis.set_major_locator(y_major_locator)  # 把y轴的主刻度设置为0.001的倍数
    # 设置刻度范围
    # plt.xlim(0, 50000)  # 把x轴的刻度范围设置为-500到15000，因为500不满一个刻度间隔，所以数字不会显示出来，但是能看到一点空白
    # plt.ylim(0, 0.00016)  # 把y轴的刻度范围设置为-0.0005到1，同理，-0.0005不会标出来，但是能看到一点空白

    # plt.legend(loc=0)

    plt.grid()  # 添加网格
    plt.savefig(dir_name + name + '均值' + '.png')
    plt.show()


def get_xy(dir_name, file_name, name_arr):
    path = dir_name + file_name
    nano_avg = np.loadtxt(path)
    l = len(name_arr)
    avg = nano_avg[:, :l]
    # x轴坐标
    xx = np.arange(len(avg))
    for i in range(0, l):
        # y轴坐标
        yy = avg[:, i]
        name = name_arr[i]
        plotter(xx, yy, name,dir_name, file_name)



def main():
    dir_name = "../result/NANO/nano_avg/"
    file_name = "avg.dat"
    name_arr = ["ca", "cag", "caf"]
    get_xy(dir_name, file_name, name_arr)
    print("plotterfn 执行结束")


if __name__ == '__main__':
    mlp.use('TkAgg')
    main()
