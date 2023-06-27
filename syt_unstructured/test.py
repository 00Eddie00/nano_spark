import unittest
import numpy as np
import os
from scipy.spatial import Delaunay
from shapely.geometry import Point, LineString, Polygon, MultiPolygon
from math import exp, sqrt, log


def compare_file():
    array1 = np.load('optical_blurring/refined_relations.npy')
    array2 = np.load('optical_blurring/refined_relations_xmh.npy')
    for i in range(601):
        for j in range(601):
            if array1[i, j, 0] == 3:
                print(f"relation[{i}, {j}, 1]={array1[i, j, 1]}")
    # for i in range(601):
    #     for j in range(601):
    #         if array1[i, j, 0] == 2:
    #             count = count + 1
    # print(count)
    # for i in range(601):
    #     for j in range(601):
    #         if array2[i, j, 0] == 3:
    #             array2[i, j, 0] = 2
    # for i in range(601):
    #     for j in range(601):
    #         # if array1[i, j, 0] == -1:
    #         #     print(f"array1[{i},{j},0]={array1[i, j, 0]},array1[{i},{j},1]={array1[i, j, 1]}")
    #         # if array2[i, j, 0] == -1:
    #         #     print(f"array2[{i},{j},0]={array2[i, j, 0]},array2[{i},{j},1]={array2[i, j, 1]}")
    #         # if array1[i][j][1] != array2[i][j][1]:
    #         #     print(f"array1[{i}][{j}][1]={array1[i][j][1]}!=array2[{i}][{j}][1]={array2[i][j][1]}")
    #         # if not np.isclose(array1[i][j][1], array2[i][j][1]) :
    #         #     print(f"array1[{i}][{j}][1]={array1[i][j][1]}!=array2[{i}][{j}][1]={array2[i][j][1]}")
    #         if array1[i][j][0] == 2 and array2[i][j][0] == 2 :
    #             continue
    #         elif array1[i][j][0] == 1 and array2[i][j][0] == 1:
    #             continue
    #         elif array1[i][j][0] == -1 and array2[i][j][0] == -1:
    #             continue
    #         else:
    #             print(f"array1[{i}][{j}][1]={array1[i][j][1]}\tarray2[{i}][{j}][1]={array2[i][j][1]}")
    # # 检查两个数组是否相同
    # are_equal = np.array_equal(array1, array2)
    #
    # if are_equal:
    #     print("两个.npy文件相同")
    # else:
    #     print("两个.npy文件不相同")


class MyTestCase(unittest.TestCase):
    def test_something(self):
        compare_file()


if __name__ == '__main__':
    unittest.main()
