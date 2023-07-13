import unittest
import numpy as np
import os
from scipy.spatial import Delaunay
from shapely.geometry import Point, LineString, Polygon, MultiPolygon
from math import exp, sqrt, log
import matplotlib.pyplot as plt
import matplotlib


def find_point_index(array, point):
    # # 将数组转换为NumPy数组
    # np_array = np.array(array)
    #
    # # 在数组中寻找匹配的索引
    # index = np.where((np_array == point).all(axis=1))
    # print(index)
    #
    # if index[0].size > 0:
    #     return index[0][0]
    # else:
    #     return None
    # 生成一个示例的四维数组
    # arr = np.random.randint(0, 10, size=(3, 4, 2, 2))
    arr = np.array([
        [[1, 2], [5, 4], [7, 8]],
        [[10, 11], [4, 5], [16, 17]],
        [[19, 20], [5, 4], [4, 26]]
    ])
    # 要查找的值
    a = 4
    b = 5

    # 使用 np.where 找到值所在位置的下标
    indices = np.where((arr[:, :, 0] == a) & (arr[:, :, 1] == b))
    i, j = indices[0][0], indices[1][0]

    print("Value (a, b) = ({}, {}) is located at index (i, j) = ({}, {})".format(a, b, i, j))


class MyTestCase(unittest.TestCase):
    def test_something(self):
        # 定义二维数组
        array = [[1, 2],
                 [3, 4],
                 [5, 6],
                 [7, 8]]

        # 定义要查找的点
        point = [5, 6]

        # 查找点的下标
        find_point_index(array, point)

        # if index is not None:
        #     print(f"点所在位置的下标为: {index}")
        # else:
        #     print("找不到该点在数组中的位置")


if __name__ == '__main__':
    unittest.main()
