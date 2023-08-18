import unittest
import numpy as np
from scipy.interpolate import griddata
import time


def cha1():
    # 记录程序开始时间
    start_time = time.time()
    # 假设四个顶点的坐标和对应的值
    vertices = np.array([(0, 0), (0, 1), (1, 0), (1, 1)])
    values = np.array([10, 20, 30, 40])
    # 要插入点的坐标
    insert_point = np.array([0.43, 0.86])
    # 使用griddata进行双线性插值
    inserted_value = griddata(vertices, values, insert_point, method='linear')
    # 记录程序结束时间
    end_time = time.time()
    # 计算程序运行时间
    elapsed_time = end_time - start_time
    # print(f"程序运行时间: {elapsed_time:.6f} 秒")
    elapsed_time_ms = elapsed_time * 1000  # 将秒转换为毫秒
    print(f"程序1运行时间: {elapsed_time_ms:.3f} 毫秒")
    print("插入点的值1:", inserted_value[0])


def cha2():
    # 记录程序开始时间
    start_time = time.time()
    radius, height = 0.43, 0.86
    lower_r, upper_r, lower_z, upper_z = 0, 1, 0, 1
    alpha = (radius - lower_r) / (upper_r - lower_r)
    beta = (height - lower_z) / (upper_z - lower_z)
    v1, v2, v3, v4 = 10, 30, 40, 20
    interpolated_value = (1 - alpha) * (1 - beta) * v1 + alpha * (1 - beta) * v2 + alpha * beta * v3 + (
            1 - alpha) * beta * v4
    # 记录程序结束时间
    end_time = time.time()
    # 计算程序运行时间
    elapsed_time = end_time - start_time
    # print(f"程序运行时间: {elapsed_time:.6f} 秒")
    elapsed_time_ms = elapsed_time * 1000  # 将秒转换为毫秒
    print(f"程序2运行时间: {elapsed_time_ms:.3f} 毫秒")
    print("插入点的值2:", interpolated_value)


class MyTestCase(unittest.TestCase):
    def test_something(self):
        # 假设 coordinates_array 是你的坐标数组
        coordinates_array = np.array([[1, 2, 3],
                                      [4, 5, 6],
                                      [7, 8, 9]])

        # 创建一个哈希表（字典）用于存储坐标值与索引的映射
        coordinates_dict = {(x, y, z): i for i, (x, y, z) in enumerate(coordinates_array)}

        # 要查找的坐标值
        target_coordinates = (4, 5, 6)

        # 使用哈希表直接查找索引
        if target_coordinates in coordinates_dict:
            index = coordinates_dict[target_coordinates]
            print("Target coordinates:", target_coordinates)
            print("Index of target coordinates:", index)
        else:
            print("Target coordinates not found in the array.")


if __name__ == '__main__':
    unittest.main()
