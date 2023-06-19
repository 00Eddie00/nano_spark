import unittest
import numpy as np
from scipy.signal import convolve2d
from scipy.ndimage import convolve1d


# 三维卷积
def convolve3d(original_matrix, kernel, xy_c_val, z_c_val):
    matrix = original_matrix.copy()
    # 第一维
    for i in np.arange(matrix.shape[1]):
        for j in np.arange(matrix.shape[2]):
            one_line = matrix[:, i, j]
            matrix[:, i, j] = convolve1d(one_line, kernel[0], mode='constant', cval=xy_c_val)
    # 第二维
    for i in np.arange(matrix.shape[0]):
        for j in np.arange(matrix.shape[2]):
            one_line = matrix[i, :, j]
            matrix[i, :, j] = convolve1d(one_line, kernel[1], mode='constant', cval=xy_c_val)
    # 第三维
    for i in np.arange(matrix.shape[0]):
        for j in np.arange(matrix.shape[1]):
            one_line = matrix[i, j, :]
            matrix[i, j, :] = convolve1d(one_line, kernel[2], mode='constant', cval=z_c_val)
    return matrix


# 三维卷积
def convolve3d2(original_matrix, kernel, xy_c_val, z_c_val):
    matrix = original_matrix.copy()
    c_val = xy_c_val
    i = 0
    for kernel_i in kernel:
        if i == 2:
            c_val = z_c_val
        # 在每个维度上进行一维卷积操作
        result = convolve1d(matrix, kernel_i, axis=i, mode='constant', cval=c_val)
        matrix = result
        i = i + 1
    return matrix


class MyTestCase(unittest.TestCase):
    def test_something(self):
        C_VAL = 4.081632653061224858e-03  # 用于钙火花
        # 创建两个三维数组
        concentration = np.random.rand(601, 601, 16)
        # concentration = np.full((601, 601, 16), 20)
        kernel = np.load("optical_blurring/kernel_v1.npy", allow_pickle=True)
        matrix1 = convolve3d(concentration, kernel, C_VAL, 0)
        print(concentration[456, 356, 6])
        print("******************************")
        print(matrix1[456, 356, 6])
        print("******************************")
        matrix2 = convolve3d2(concentration, kernel, C_VAL, 0)
        print(matrix2[456, 356, 6])


if __name__ == '__main__':
    unittest.main()
