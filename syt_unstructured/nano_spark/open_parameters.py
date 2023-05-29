import numpy as np

open_grid_file_name = "../config/open/open_grid_coordinates.csv"
neighbors_file_name = "../config/open/open_neighbor.csv"
coefficient_file_name = "../config/open/coefficient.csv"

# 亚空间与大开放空间交界点
boundary_r = 300
boundary_z = [0, 3, 7.5]
# 大开放区域
big_z_min = 0
big_z_max = 2500
big_r_min = 300
big_r_max = 5000
big_z_length = 69  # z轴点的个数
big_open_count = 4623
# 小开放区域
small_z_min = 13
small_z_max = 2500
small_r_min = 0
small_r_max = 295
small_z_length = 66
