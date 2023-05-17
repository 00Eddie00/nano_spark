import numpy as np
from generate_grid import generate_coordinates


def generate_neighbor():
    grid_coordinates = np.loadtxt("config/grid/grid_coordinates.csv", delimiter=",")
    point_nums = len(grid_coordinates)
    # 建立一个二维数组存储各点的邻居；索引0：z轴正方向，索引1：z轴负方向，索引2：ρ轴正方向，索引3：ρ负方向；
    # -1 代表此方向上无邻居；-2代表此方向有两个邻居，需要求中点。
    neighbor = np.full((point_nums, 4), -1)
    two_neighbor = np.empty((23, 3))
    nano_coordinates, open_coordinates, boundary_coordinates, new_open_coordinates1, new_open_coordinates2, nano_length, boundary_length, open_length = generate_coordinates()
    grid_coordinates = nano_coordinates + boundary_coordinates + open_coordinates
    open_max = list(max(grid_coordinates))
    # 纳米空间
    for i in range(nano_length):
        # z轴邻居
        if 7.5 > grid_coordinates[i][0] > -7.5:
            neighbor[i][0] = i + 1
            neighbor[i][1] = i - 1
        elif grid_coordinates[i][0] == 7.5:
            neighbor[i][1] = i - 1
        else:
            neighbor[i][0] = i + 1
        # r轴邻居
        if 296 > grid_coordinates[i][1] > 0:
            neighbor[i][2] = i + 24
            neighbor[i][3] = i - 24
        elif grid_coordinates[i][1] == 0:
            neighbor[i][2] = i + 24
        else:
            neighbor[i][3] = i - 24
            for j in range(nano_length, nano_length + boundary_length):
                if grid_coordinates[i][0] == grid_coordinates[j][0]:
                    neighbor[i][2] = j
                    neighbor[j][3] = i
    # 开放空间
    for i in range(nano_length + boundary_length, point_nums):
        # z轴邻居
        if open_max[0] > grid_coordinates[i][0] > -open_max[0]:
            neighbor[i][0] = i + 1
            neighbor[i][1] = i - 1
        elif grid_coordinates[i][0] == open_max[0]:
            neighbor[i][1] = i - 1
        else:
            neighbor[i][0] = i + 1

        # r轴邻居
        if open_max[1] > grid_coordinates[i][1] > 305:
            neighbor[i][2] = i + 133
            neighbor[i][3] = i - 133
        elif grid_coordinates[i][1] == open_max[1]:
            neighbor[i][3] = i - 133
        else:
            neighbor[i][2] = i + 133
            for j in range(nano_length, nano_length + boundary_length):
                if grid_coordinates[i][0] == grid_coordinates[j][0]:
                    neighbor[i][3] = j
                    neighbor[j][2] = i

    k = 0
    # 边界
    for i in range(nano_length, nano_length + boundary_length):
        # z轴邻居
        if open_max[0] > grid_coordinates[i][0] > -open_max[0]:
            neighbor[i][0] = i + 1
            neighbor[i][1] = i - 1
        elif -open_max[0] == grid_coordinates[i][0]:
            neighbor[i][0] = i + 1
        else:
            neighbor[i][1] = i - 1

        # r轴邻居
        if grid_coordinates[i][0] < -7.5:
            neighbor[i][2] = i + 155
        elif grid_coordinates[i][0] > 7.5:
            neighbor[i][2] = i + 133
        elif grid_coordinates[i][0] <= -5.5:
            neighbor[i][2] = -2
            two_neighbor[k][0] = i
            two_neighbor[k][1] = 3579
            two_neighbor[k][2] = 3580
            k = k + 1
        elif grid_coordinates[i][0] == -5:
            neighbor[i][3] = -2
            two_neighbor[k][0] = i
            two_neighbor[k][1] = 3338
            two_neighbor[k][2] = 3339
            k = k + 1
        elif grid_coordinates[i][0] <= -0.5:
            neighbor[i][2] = -2
            two_neighbor[k][0] = i
            two_neighbor[k][1] = 3580
            two_neighbor[k][2] = 3581
            k = k + 1
        elif 4.5 >= grid_coordinates[i][0] >= 0.5:
            neighbor[i][2] = -2
            two_neighbor[k][0] = i
            two_neighbor[k][1] = 3581
            two_neighbor[k][2] = 3582
            k = k + 1
        elif 7.5 >= grid_coordinates[i][0] >= 5.5:
            neighbor[i][2] = -2
            two_neighbor[k][0] = i
            two_neighbor[k][1] = 3582
            two_neighbor[k][2] = 3583
            k = k + 1
    # print(k)
    np.savetxt("config/grid/neighbor.csv", neighbor, fmt='%d', delimiter=",")
    np.savetxt("config/grid/two_neighbor.csv", two_neighbor, fmt='%d', delimiter=",")


if __name__ == '__main__':
    generate_neighbor()
