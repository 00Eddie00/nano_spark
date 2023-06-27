import numpy as np
from ob_parameters import *


def gaussian_function(w, sigma):
    return np.power((2 * np.pi * sigma), -1 / 2) * np.exp(-np.square(w) / (2 * sigma))


# 生成卷积核
def generate_kernel():
    sigma = sigma_xy
    # x、y方向上
    xy_parameter = np.array([gaussian_function(i, sigma) for i in range(xy_count + 1)])
    sigma = sigma_z
    # z方向上
    z_parameter = np.array([gaussian_function(i, sigma) for i in range(z_count + 1)])
    # 拼接
    xy_parameter = np.concatenate((xy_parameter[xy_count + 1:0:-1], xy_parameter))
    z_parameter = np.concatenate((z_parameter[z_count + 1:0:-1], z_parameter))

    kernel = np.asarray([xy_parameter, xy_parameter, z_parameter], dtype='object')
    return kernel


def main():
    kernel = generate_kernel()
    np.save("kernel_v1", kernel)


if __name__ == "__main__":
    main()
