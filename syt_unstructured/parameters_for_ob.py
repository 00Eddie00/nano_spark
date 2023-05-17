from math import log

# 用于处理浓度矩阵
nano_point_totals = 3360  # 点总数

xy_length = 300  # x、y方向的长度的一半，单位nm
z_length = 15  # z方向的长度，单位nm

C_VAL = 0.0

# 用于产生卷积核
# xy_count = 517  # (卷积核矩阵在x、y维度上点个数-1)/2
xy_count = 697
# z_count = 1033  # (卷积核矩阵在z维度上点个数-1)/2
z_count = 15  # (卷积核矩阵在z维度上点个数-1)/2
# z_count = 1033
# sigma_xy = 400 ** 2 / (8 * log(2))
sigma_xy = 540 ** 2 / (8 * log(2))
# sigma_z = 800 ** 2 / (8 * log(2))
sigma_z = 540 ** 2 / (8 * log(2))
