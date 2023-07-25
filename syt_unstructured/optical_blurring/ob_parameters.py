import numpy as np

sigma_xy = np.square(400) / (8 * np.log(2))
sigma_z = np.square(800) / (8 * np.log(2))
xy_count = 517  # (卷积核矩阵在x、y维度上点个数-1)/2
z_count = 15  # (卷积核矩阵在z维度上点个数-1)/2
length_side = 1000
half_length = 1000 / 2
