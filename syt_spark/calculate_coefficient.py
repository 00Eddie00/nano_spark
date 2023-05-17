import numpy as np


def cal_r(point_index, two_neighbor):
    r = -1
    for i in range(len(two_neighbor)):
        if point_index == two_neighbor[i][0]:
            r = two_neighbor[i][1]
            return r
    return r


def cal_coe():
    grid_coordinates = np.loadtxt("config/grid/grid_coordinates.csv", delimiter=",")  # 点坐标
    neighbor = np.loadtxt("config/grid/neighbor.csv", int, delimiter=",")  # 邻点
    two_neighbor = np.loadtxt("config/grid/two_neighbor.csv", int, delimiter=",")  # 交界处的特殊邻点
    coefficient = np.empty((len(grid_coordinates), 4), float)
    for i in range(len(grid_coordinates)):
        z2 = grid_coordinates[i][0]
        r2 = grid_coordinates[i][1]
        p3, p7, p1, p5 = neighbor[i]  # 上下外内
        #  z轴系数
        if p3 == -1:  # 没上
            z1 = grid_coordinates[p7][0]
            h = z2 - z1
            coefficient[i][0] = 0.0
            coefficient[i][1] = 1 / (2 * h * (z1 - z2))
        elif p7 == -1:  # 没下
            z3 = grid_coordinates[p3][0]
            h = z3 - z2
            coefficient[i][0] = 1 / (2 * h * (z3 - z2))
            coefficient[i][1] = 0.0
        else:
            z1 = grid_coordinates[p7][0]
            z3 = grid_coordinates[p3][0]
            h = z3 - z1
            coefficient[i][0] = 1 / (2 * h * (z3 - z2))
            coefficient[i][1] = 1 / (2 * h * (z1 - z2))
        #  r轴系数
        if p1 == -1:  # 没外
            r1 = grid_coordinates[p5][1]
            coefficient[i][2] = 0.0
            coefficient[i][3] = r1 / ((r2 * r2 - r1 * r1) * (r1 - r2))
        elif p1 == -2:  # 没对应的外
            r = cal_r(i, two_neighbor)
            r1 = grid_coordinates[p5][1]
            r3 = grid_coordinates[r][1]
            coefficient[i][2] = r3 / ((r3 * r3 - r1 * r1) * (r3 - r2))
            coefficient[i][3] = r1 / ((r3 * r3 - r1 * r1) * (r1 - r2))
        elif p5 == -1:  # 没内
            r3 = grid_coordinates[p1][1]
            coefficient[i][2] = r3 / ((r3 * r3 - r2 * r2) * (r3 - r2))
            coefficient[i][3] = 0.0
        elif p5 == -2:  # 没对应的内
            r = cal_r(i, two_neighbor)
            r1 = grid_coordinates[r][1]
            r3 = grid_coordinates[p1][1]
            coefficient[i][2] = r3 / ((r3 * r3 - r1 * r1) * (r3 - r2))
            coefficient[i][3] = r1 / ((r3 * r3 - r1 * r1) * (r1 - r2))
        else:
            r1 = grid_coordinates[p5][1]
            r3 = grid_coordinates[p1][1]
            coefficient[i][2] = r3 / ((r3 * r3 - r1 * r1) * (r3 - r2))
            coefficient[i][3] = r1 / ((r3 * r3 - r1 * r1) * (r1 - r2))
    np.save("config/coefficient/coefficient", coefficient)


if __name__ == '__main__':
    cal_coe()
