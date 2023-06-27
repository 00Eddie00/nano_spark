import numpy as np


# 计算每个三角形的系数b c n L
def cal_elements(grid_file_name, nod_file_name):
    grid = np.loadtxt(grid_file_name, dtype=np.float64)
    nod = np.loadtxt(nod_file_name, dtype=int) - 1
    NE = len(nod)
    NP = len(grid)
    # 每个三角形面积
    single_area = np.zeros(NE, float)

    # 每个点的控制面积
    control_area = np.zeros(NP, float)

    #  中间变量
    cal_max = np.zeros(NP, int)
    area_matrix = np.zeros([3, 3], float)
    # nix*l
    nix_multiply_l = np.zeros([NE, 3], float)

    # niy*l
    niy_multiply_l = np.zeros([NE, 3], float)

    # ai
    a_arr = np.zeros([NE, 3], float)

    # bi
    b_arr = np.zeros([NE, 3], float)

    # ci
    c_arr = np.zeros([NE, 3], float)

    # 非结构网格总面积
    total_area = 0.0

    for i in range(0, NE):
        for j in range(0, 3):
            p = nod[i, j]
            area_matrix[j, 0] = 1.0
            area_matrix[j, 1] = grid[p, 0]
            area_matrix[j, 2] = grid[p, 1]
        area = np.linalg.det(area_matrix) * 0.5
        single_area[i] = area
        total_area = total_area + area
        for v in range(0, 3):
            p = nod[i, v]
            control_area[p] = control_area[p] + area
            cal_max[p] = cal_max[p] + 1
    nmax = int(max(cal_max))
    # 该点的相邻三角形编号
    near_triangle = np.full([NP, nmax], -1, int)
    # 当前点在当前三角形中的编号
    index_in_triangle = np.full([NP, nmax], -1, int)
    # 记录某点在其相邻三角形中序列中的第几个三角形
    num_temp = np.zeros(NP, int)
    for j in range(0, NE):
        for k in range(0, 3):
            i = nod[j, k]
            t = num_temp[i]
            near_triangle[i, t] = j
            index_in_triangle[i, t] = k
            num_temp[i] = num_temp[i] + 1
            # 当前点对应的那条边的b c n*l
            nod2 = nod[j, (k + 1) % 3]
            nod3 = nod[j, (k + 2) % 3]
            # nix_multiply_l[j, k] = grid[nod2, 1] - grid[nod3, 1]
            # niy_multiply_l[j, k] = grid[nod3, 0] - grid[nod2, 0]
            # a_arr[j, k] = (grid[nod3, 0] * grid[nod2, 1] - grid[nod2, 0] * grid[nod3, 1]) / (2 * single_area[j])
            # b_arr[j, k] = (grid[nod3, 1] - grid[nod2, 1]) / (2 * single_area[j])
            # c_arr[j, k] = (grid[nod2, 0] - grid[nod3, 0]) / (2 * single_area[j])
            nix_multiply_l[j, k] = grid[nod3, 1] - grid[nod2, 1]
            niy_multiply_l[j, k] = grid[nod2, 0] - grid[nod3, 0]
            a_arr[j, k] = (grid[nod2, 0] * grid[nod3, 1] - grid[nod3, 0] * grid[nod2, 1]) / (2 * single_area[j])
            b_arr[j, k] = (grid[nod2, 1] - grid[nod3, 1]) / (2 * single_area[j])
            c_arr[j, k] = (grid[nod3, 0] - grid[nod2, 0]) / (2 * single_area[j])
    return single_area, control_area, near_triangle, index_in_triangle, nix_multiply_l, niy_multiply_l, a_arr, b_arr, c_arr, nmax, total_area


# 判断三角形类型
def judge_point(n, nod, npoch):
    in_boundary, out_boundary, inner = [], [], []
    for j in range(0, 3):
        k = nod[n, j]
        if npoch[k] == 2:
            in_boundary.append(k)
        elif npoch[k] == 4:
            out_boundary.append(k)
        elif npoch[k] == 0:
            inner.append(k)
    return in_boundary, out_boundary, inner


# 计算控制面某段边界长度
def cal_length(boundary, grid):
    a = [grid[boundary[0], 0], grid[boundary[0], 1]]
    b = [grid[boundary[1], 0], grid[boundary[1], 1]]
    return np.linalg.norm(np.array(a) - np.array(b))


# 计算各个点平均值
def cal_avg(concentration, grid_file_name, nod_file_name):
    single_area, control_area, near_triangle, index_in_triangle, nix_multiply_l, niy_multiply_l, a_arr, b_arr, c_arr, nmax, total_area = cal_elements(
        grid_file_name, nod_file_name)
    total = 0.0
    NP = len(concentration)
    for i in range(0, NP):
        total = total + concentration[i] * control_area[i]
    total = total / 3.0
    avg = total / total_area
    return avg


# 计算出流边界长度
def cal_out_boundary_length(grid):
    out_boundary = grid[84:205]
    out_boundary_num = len(out_boundary)
    out_boundary_length = np.empty(out_boundary_num)
    for i in range(out_boundary_num):
        a = [out_boundary[i, 0], out_boundary[i, 1]]
        b = [out_boundary[(i + 1) % out_boundary_num, 0], out_boundary[(i + 1) % out_boundary_num, 1]]
        out_boundary_length[i] = np.linalg.norm(np.array(a) - np.array(b))
    return out_boundary_length


def main():
    grid_file_name = "../config/nano/4RYRgridt.dat"
    nod_file_name = "../config/nano/4RYRnod.dat"
    npoch_file_name = "../config/nano/4RYRnpoch.dat"
    npoch = np.loadtxt(npoch_file_name, dtype=int)
    grid = np.loadtxt(grid_file_name, dtype=np.float64)
    nod = np.loadtxt(nod_file_name, dtype=int) - 1
    single_area, control_area, near_triangle, index_in_triangle, nix_multiply_l, niy_multiply_l, a_arr, b_arr, c_arr, nmax, total_area = cal_elements(
        grid_file_name, nod_file_name)
    in_boundary, out_boundary, inner = judge_point(5, nod, npoch)
    print(in_boundary)
    print(out_boundary)
    print(inner)


if __name__ == '__main__':
    main()
