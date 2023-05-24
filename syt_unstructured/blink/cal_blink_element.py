import numpy as np
from syt_unstructured.blink.blink_parameters import *
from syt_unstructured.tool.cal_bcnl import *


def fn(f, g):
    global NEWCCA
    tempcca = np.empty(NP, float)  # 保存每一次迭代后的结果
    single_area, control_area, near_triangle, index_in_triangle, nix_multiply_l, niy_multiply_l, b_arr, c_arr, nmax, total_area = cal_elements(
        grid_file_name, nod_file_name)
    for it in range(0, 10):
        for p in range(0, NP):
            m1, m5, m6 = 0.0, 0.0, 0.0  # 分别对应(1)(5)(6)式
            m2 = DT * (K2 * g[p] - K1 * f[p] * (F - g[p])) * control_area[p]  # 计算（2）式
            n7, n8, n9, n10 = 0.0, 0.0, 0.0, 0.0  # 分别对应(7)(8)(9)(10)式
            in_boundary_list = []  # 保存有入流边界的三角形编号
            out_boundary_list = []  # 保存有出流边界的三角形编号
            for o in range(0, nmax):  # 计算相邻三角形
                n = near_triangle[p, o]  # n为三角形编号
                if n == -1:
                    break
                i = index_in_triangle[p, o]
                nod_2, nod_3 = nod[n, (i + 1) % 3], nod[n, (i + 2) % 3]  # 其余两点的编号及在三角形中的序号
                f1i, f2i, f3i = f[p], f[nod_2], f[nod_3]  # f1i, f2i, f3i取n时刻的值

                if it > 0:  # 判断是否为第一次迭代
                    f1i, f2i, f3i = (f[p] + NEWCCA[p]) / 2.0, (f[nod_2] + NEWCCA[nod_2]) / 2.0, (
                            f[nod_3] + NEWCCA[nod_3]) / 2.0
                avg = (1 / 3) * (f1i + f2i + f3i)  # fni为第 i 个三角形单元的三个顶点 n 时刻的平均值
                m5 = m5 + (1.0 + ((BCSQ * KDCSQ) / (KDCSQ + avg) ** 2)) * single_area[n]  # 计算（5）式
                in_boundary, out_boundary, inner = judge_point(n, nod, npoch)
                if len(in_boundary) == 2:
                    in_boundary_list.append(o)  # 记录该点的triangleNumber的下标
                elif len(out_boundary) == 2:
                    out_boundary_list.append(o)
                if npoch[nod_2] > B_INNER and npoch[nod_3] > B_INNER:
                    continue
                m1 = m1 + (nix_multiply_l[n, i] * (f2i * b_arr[n, (i + 1) % 3] + f3i * b_arr[n, (i + 2) % 3]) +
                           niy_multiply_l[n, i] * (f2i * c_arr[n, (i + 1) % 3] + f3i * c_arr[n, (i + 2) % 3]))  # 计算（1）式
                m6 = m6 + (nix_multiply_l[n, i] * b_arr[n, i] + niy_multiply_l[n, i] * c_arr[n, i])  # 计算（6）式
            tempcca[p] = (DT * DCAJSR * m1 + m2 + f[p] * m5 + 0.5 * f[p] * DT * DCAJSR * m6) / (
                    m5 - 0.5 * DT * DCAJSR * m6)  # 存入临时数组保存
            if len(in_boundary_list) > 0:  # 入流边界数目大于0
                for q in range(0, len(in_boundary_list)):
                    o1 = in_boundary_list[q]
                    in_n = near_triangle[p, o1]  # 取出三角形编号
                    i = index_in_triangle[p, o1]
                    in_boundary, out_boundary, inner = judge_point(in_n, nod, npoch)
                    nod_2, nod_3 = nod[in_n, (i + 1) % 3], nod[in_n, (i + 2) % 3]
                    f2j, f3j = f[nod_2], f[nod_3]
                    if it > 0:  # 判断是否为第一次迭代
                        f2j, f3j = NEWCCA[nod_2], NEWCCA[nod_3]
                    Lj = cal_length(in_boundary, grid)  # 入流边界长度
                    n7 = n7 + DCAFSR * (CCAFSR - (f2j + f3j) / 3.0) * Lj  # 计算（7）式
                    n8 = n8 + DCAFSR * Lj  # 计算（8）式
                tempcca[p] = (DT * DCAJSR * m1 + m2 + f[p] * m5 + 0.5 * f[p] * DT * DCAJSR * m6 + n7 * DT) / (
                        m5 - 0.5 * DT * DCAJSR * m6 + (n8 * DT) / 3.0)  # 若有入流边界则修改临时数组
            elif len(out_boundary_list) > 0:  # 出流边界数目大于0
                for r in range(0, len(out_boundary_list)):
                    o2 = out_boundary_list[r]
                    out_n = near_triangle[p, o2]  # 三角形编号
                    i = index_in_triangle[p, o2]
                    in_boundary, out_boundary, inner = judge_point(out_n, nod, npoch)
                    nod_2, nod_3 = nod[out_n, (i + 1) % 3], nod[out_n, (i + 2) % 3]
                    f2k, f3k = f[nod_2], f[nod_3]
                    if it > 0:
                        f2k, f3k = NEWCCA[nod_2], NEWCCA[nod_3]
                    Lk = cal_length(out_boundary, grid)  # 出流边界长度
                    n9 = n9 + DCARYR * (CCAMYO - (f2k + f3k) / 3.0) * Lk  # 计算（9）式
                    n10 = n10 + DCARYR * Lk  # 计算（10）式
                tempcca[p] = (DT * DCAJSR * m1 + m2 + f[p] * m5 + 0.5 * f[p] * DT * DCAJSR * m6 + n9 * DT) / (
                        m5 - 0.5 * DT * DCAJSR * m6 + (n10 * DT) / 3.0)  # 若有出流边界则修改临时数组
        NEWCCA = np.copy(tempcca)
    return NEWCCA


def gn(f, g):
    global NEWF
    tempf = np.empty(NP, float)  # 保存每一次迭代后的结果
    single_area, control_area, near_triangle, index_in_triangle, nix_multiply_l, niy_multiply_l, b_arr, c_arr, nmax, total_area = cal_elements(
        grid_file_name, nod_file_name)
    for it in range(0, 10):
        for p in range(0, NP):
            m1, m2, m6 = 0.0, 0.0, 0.0,  # 分别对应(1)(2)(6)式
            for o in range(0, nmax):
                n = near_triangle[p, o]  # 三角形编号
                if n == -1:
                    break
                i = index_in_triangle[p, o]
                nod_2, nod_3 = nod[n, (i + 1) % 3], nod[n, (i + 2) % 3]
                if npoch[nod_2] > B_INNER and npoch[nod_3] > B_INNER:
                    continue
                g2i, g3i = g[nod_2], g[nod_3]
                if it > 0:  # 判断是否为第一次迭代
                    g2i, g3i = (g[nod_2] + NEWF[nod_2]) / 2.0, (g[nod_3] + NEWF[nod_3]) / 2.0
                m1 = m1 + DT * DF * (
                            nix_multiply_l[n, i] * (g2i * b_arr[n, (i + 1) % 3] + g3i * b_arr[n, (i + 2) % 3]) +
                            niy_multiply_l[n, i] * (g2i * c_arr[n, (i + 1) % 3] + g3i * c_arr[n, (i + 2) % 3]))
                m6 = m6 + DT * DF * 0.5 * (nix_multiply_l[n, i] * b_arr[n, i] + niy_multiply_l[n, i] * c_arr[n, i])
            m2 = DT * (K1 * f[p] * (F - g[p]) - K2 * g[p]) * control_area[p]
            tempf[p] = (m1 + m2 + g[p] * control_area[p] + m6 * g[p]) / (
                    control_area[p] - m6)
        NEWF = np.copy(tempf)
    return NEWF
