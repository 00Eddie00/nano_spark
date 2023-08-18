import unittest
import numpy as np
import os
from scipy.spatial import Delaunay
from shapely.geometry import Point, LineString, Polygon, MultiPolygon
import matplotlib.pyplot as plt
import matplotlib
import math


def find_multiple_dimensions_index(array, point):
    r, z = point[0], point[1]

    # 在数组中寻找匹配的索引
    indices = np.where((array[:, :, :, 0] == r) & (array[:, :, :, 1] == z))

    i, j, k = indices[0][0], indices[1][0], indices[2][0]
    return i, j, k


def find_point_index(array, point):
    # 在数组中寻找匹配的索引
    index = np.where((array == point).all(axis=1))

    if index[0].size > 0:
        return index[0][0]
    # else:
    #     return None


# 插值计算
def interpolation_calculation(lower_r, upper_r, lower_z, upper_z, radius, height, grid_coordinates):
    r1, z1 = lower_r, lower_z  # 左下
    r2, z2 = upper_r, lower_z  # 右下
    r3, z3 = upper_r, upper_z  # 右上
    r4, z4 = lower_r, upper_z  # 左上
    if lower_z == 7.5 and radius < 300:
        index1 = np.int32(math.floor(radius))
        index2 = np.int32(math.floor(radius)) + 1
        index3 = find_point_index(grid_coordinates, [z3, r3])
        index4 = find_point_index(grid_coordinates, [z4, r4])
        # 有一丢丢问题
    elif radius == 300 and height <= 13.0:
        index1 = [np.int32(r1) + 300, 0, 8]
        index2 = find_point_index(grid_coordinates, [z2, r2])
        index3 = find_point_index(grid_coordinates, [z3, r3])
        if height <= 7.0:
            index4 = [np.int32(r1) + 300, 0, 8]
        else:
            index4 = find_point_index(grid_coordinates, [z4, r4])
    else:
        index1 = find_point_index(grid_coordinates, [z1, r1])
        index2 = find_point_index(grid_coordinates, [z2, r2])
        index3 = find_point_index(grid_coordinates, [z3, r3])
        index4 = find_point_index(grid_coordinates, [z4, r4])
    if index1 == None or index2 == None or index3 == None or index4 == None:
        print(f"radius={radius}height={height}")
        print(f"index1={index1}index2={index2}index3={index3}index4={index4}")
        print(f"z1={z1}r1={r1}")
        print(f"z4={z4}r4={r4}")


half_length = 500


def generate_interval(dis):
    # 每隔1nm取一个点
    interval_arr = [5, 10, 20, 25, 50]
    end_value = 0.0
    interval = half_length // len(interval_arr)
    process_points = np.empty(0)
    for i in range(len(interval_arr)):
        start_value = end_value
        end_value = end_value + interval
        num_points = interval // interval_arr[i] + 1
        data_array = np.linspace(start_value, end_value, np.int(num_points))
        if i != len(interval_arr) - 1:
            data_array = data_array[:-1]
        process_points = np.concatenate((process_points, data_array))

    x_half_left = dis - np.flip(process_points)
    x_half_right = process_points + dis
    x_arr = np.concatenate((x_half_left[:-1], x_half_right))

    y_half_left = - np.flip(process_points)
    y_half_right = process_points
    y_arr = np.concatenate((y_half_left[:-1], y_half_right))

    z_half_left = - np.flip(process_points)
    z_half_right = process_points
    z_arr = np.concatenate((z_half_left[:-1], z_half_right))
    print(x_arr)
    print(y_arr)
    print(z_arr)

def point_scatter(start, end, base, k=1, positive=True):
    result = [start]
    flag = 1 if positive else -1
    i = 0
    while True:
        step = (k * i + base) * flag
        value = result[i] + step
        if positive is True and value > end:
            result[-1] = end
            break
        elif positive is False and value < end:
            result[-1] = end
            break
        else:
            result.append(value)
            i = i + 1
    return result


class MyTestCase(unittest.TestCase):
    def test_something(self):
        open_r = point_scatter(300, 5000, 5, k=2)
        new_open_r = np.flipud(point_scatter(295, 0, 5, k=2, positive=False))
        r = np.concatenate((new_open_r, open_r))[:17]
        print(r)


if __name__ == '__main__':
    unittest.main()
