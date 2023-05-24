import numpy as np
from open_parameters import *


def open_coefficient():
    grid_coordinates = np.loadtxt(grid_file_name, delimiter=",")  # 点坐标
    neighbors = np.loadtxt(neighbors_file_name, int, delimiter=",")  # 邻点
    point_count = len(grid_coordinates)
    # 上下外内
    coefficient = np.empty((point_count, 4))
    for i in range(point_count):
        z2 = grid_coordinates[i, 0]
        r2 = grid_coordinates[i, 1]
        # 上下外内
        coe1, coe2, coe3, coe4 = 0, 0, 0, 0
        neighbor_up, neighbor_down, neighbor_right, neighbor_left = neighbors[i, 0], neighbors[i, 1], neighbors[i, 2], \
                                                                    neighbors[i, 3]
        # 没有上顶面
        if neighbor_up == -1:
            # 则一定有下顶面
            z1 = grid_coordinates[neighbor_down, 0]
            h = z2 - z1
            coe2 = 1 / (2 * h * (z1 - z2))
        # 没有下顶面
        elif neighbor_down == -1:
            # 则一定有上顶面
            z3 = grid_coordinates[neighbor_up, 0]
            h = z3 - z2
            coe1 = 1 / (2 * h * (z3 - z2))
        else:
            z3 = grid_coordinates[neighbor_up, 0]
            z1 = grid_coordinates[neighbor_down, 0]
            h = z3 - z1
            if z2 == big_z_min:
                z3 = grid_coordinates[neighbor_up, 0]
                z1 = -z3
                h = z3 - z1
            coe1 = 1 / (2 * h * (z3 - z2))
            coe2 = 1 / (2 * h * (z1 - z2))
        # 没有外侧面
        if neighbor_right == -1:
            # 则一定有内侧面
            r1 = grid_coordinates[neighbor_left, 1]
            coe4 = r1 / ((r2 ** 2 - r1 ** 2) * (r1 - r2))
        # 没有内侧面
        elif neighbor_left == -1:
            # 则一定有外侧面
            r3 = grid_coordinates[neighbor_right, 1]
            coe3 = r3 / ((r3 ** 2 - r2 ** 2) * (r3 - r2))
        else:
            r1 = grid_coordinates[neighbor_left, 1]
            r3 = grid_coordinates[neighbor_right, 1]
            coe3 = r3 / ((r3 ** 2 - r1 ** 2) * (r3 - r2))
            coe4 = r1 / ((r3 ** 2 - r1 ** 2) * (r1 - r2))
        coefficient[i, 0], coefficient[i, 1], coefficient[i, 2], coefficient[i, 3] = coe1, -coe2, coe3, -coe4
    return coefficient


def main():
    coefficient = open_coefficient()
    np.savetxt(coefficient_file_name, coefficient, delimiter=",")


if __name__ == "__main__":
    main()
