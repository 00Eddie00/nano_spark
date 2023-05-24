import numpy as np
from open_parameters import *


def open_neighbor():
    grid_coordinates = np.loadtxt(grid_file_name, delimiter=",")
    point_count = len(grid_coordinates)
    # 上下外内
    neighbors = np.empty((point_count, 4))
    # 大空间
    for i in range(0, big_open_count):
        z, r = grid_coordinates[i][0], grid_coordinates[i][1]
        # 上下外内
        neighbor_up, neighbor_down, neighbor_right, neighbor_left = -1, -1, -1, -1
        if big_z_min < z < big_z_max:
            neighbor_up, neighbor_down = i + 1, i - 1
        elif z == big_z_min:
            neighbor_up, neighbor_down = i + 1, i + 1
        else:
            neighbor_down = i - 1

        if big_r_min < r < big_r_max:
            neighbor_right, neighbor_left = i + big_z_length, i - big_z_length
        elif r == big_r_min:
            neighbor_right = i + big_z_length
            if z not in boundary_z:
                neighbor_left = i + 5610
        else:
            neighbor_left = i - big_z_length
        neighbors[i, 0], neighbors[i, 1], neighbors[i, 2], neighbors[
            i, 3] = neighbor_up, neighbor_down, neighbor_right, neighbor_left

    # 小空间
    for i in range(big_open_count, point_count):
        z, r = grid_coordinates[i][0], grid_coordinates[i][1]
        neighbor_up, neighbor_down, neighbor_right, neighbor_left = -1, -1, -1, -1
        if small_z_min < z < small_z_max:
            neighbor_up, neighbor_down = i + 1, i - 1
        elif z == small_z_min:
            neighbor_up = i + 1
        else:
            neighbor_down = i - 1

        if small_r_min < r < small_r_max:
            neighbor_right, neighbor_left = i + small_z_length, i - small_z_length
        elif r == small_r_min:
            neighbor_right = i + small_z_length
        else:
            neighbor_right, neighbor_left = i - 5610, i - small_z_length
        neighbors[i, 0], neighbors[i, 1], neighbors[i, 2], neighbors[
            i, 3] = neighbor_up, neighbor_down, neighbor_right, neighbor_left
    return neighbors


def main():
    neighbors = open_neighbor()
    np.savetxt(neighbors_file_name, neighbors, fmt='%d', delimiter=",")


if __name__ == "__main__":
    main()
