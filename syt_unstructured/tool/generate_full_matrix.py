import numpy as np


def process_matrix(quarter_matrix):
    # 生成完整的三维矩阵
    full_matrix = np.full((2 * quarter_matrix.shape[0], 2 * quarter_matrix.shape[1], 2 * quarter_matrix.shape[2]),
                          dtype=float, fill_value=-1.0)

    # 复制八分之一矩阵的对称部分
    full_matrix[:quarter_matrix.shape[0], :quarter_matrix.shape[1], :quarter_matrix.shape[2]] = quarter_matrix

    # 反转复制的部分得到完整矩阵
    full_matrix[quarter_matrix.shape[0]:, :, :] = np.flip(full_matrix[:quarter_matrix.shape[0], :, :], axis=0)
    full_matrix[:, quarter_matrix.shape[1]:, :] = np.flip(full_matrix[:, :quarter_matrix.shape[1], :], axis=1)
    full_matrix[:, :, quarter_matrix.shape[2]:] = np.flip(full_matrix[:, :, :quarter_matrix.shape[2]], axis=2)

    print(full_matrix)


def main():
    # 在这里编写你的主要逻辑
    print("Hello, world!")
    # grids_rz = np.load("../optical_blurring/grids_rz.npy")
    # points_rz = np.load("../optical_blurring/points_rz.npy")
    # 已知的八分之一矩阵
    quarter_matrix = np.array([
        [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
        [[10, 11, 12], [13, 14, 15], [16, 17, 18]],
        [[19, 20, 21], [22, 23, 24], [25, 26, 27]]
    ])
    process_matrix(quarter_matrix)


if __name__ == "__main__":
    main()
