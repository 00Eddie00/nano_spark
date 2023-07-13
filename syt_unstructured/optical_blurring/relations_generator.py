import numpy as np
from nano_spark.open.open_parameters import *
from shapely.geometry import Point, Polygon, LineString, MultiPolygon
from tool.generate_open_grid import point_scatter


def nano_judge_relation():
    relations = np.full((601, 601, 2), dtype=int, fill_value=-1)
    # 创建点坐标数组
    nods = np.loadtxt("../config/nano/4RYRnod.dat", dtype=int) - 1
    points = np.loadtxt("../config/nano/4RYRgridt.dat", dtype=np.float64)
    triangle_items = []
    for i in range(len(nods)):
        triangle1_vertices = []
        for j in range(3):
            nod_j = nods[i, j]
            x, y = points[nod_j, 0], points[nod_j, 1]
            triangle1_vertices.append((x, y))
        triangle_i = Polygon(triangle1_vertices)
        triangle_items.append(triangle_i)
    # 创建三角形集合对象
    print("三角形集合对象创建完成")
    for i in range(601):
        print(f"i={i}")
        for j in range(601):
            x = i - 300
            y = j - 300
            if np.square(x) + np.square(y) <= 90000:
                # 创建要插入的新点对象
                new_point = Point(x, y)
                # 检查新点是否在三角形集合中的任何一个三角形内部
                for p, triangle in enumerate(triangle_items):
                    # 检查新点是否在三角形内部，包括边界（顶点）
                    if triangle.contains(new_point) or triangle.touches(new_point):
                        triangle_index = p  # 记录包含新点的三角形索引
                        relations[i, j, 0] = 2
                        relations[i, j, 1] = triangle_index
                        # 检查该点是否在顶点上
                        for k in range(3):
                            nod_k = nods[triangle_index, k]
                            xx, yy = points[nod_k, 0], points[nod_k, 1]
                            if x == xx and y == yy:
                                relations[i, j, 0] = 1
                                relations[i, j, 1] = nod_k
                                break
                        break
    return relations


def nano_relations_refine():
    relations = np.load("relations.npy")
    refined_relations = np.copy(relations)
    print("修改一些点")
    for i in range(601):
        for j in range(601):
            x = i - 300
            y = j - 300
            r_square = np.square(x) + np.square(y)
            if r_square <= 90000 and relations[i, j, 0] == -1:
                r = np.sqrt(r_square)
                if r > 299:
                    refined_relations[i, j, 0] = 1
                    refined_relations[i, j, 1] = 84
                elif r < 299:
                    refined_relations[i, j, 0] = 1
                    if x == 15 and y == 15:
                        refined_relations[i, j, 1] = 0
                    elif x == 15 and y == -15:
                        refined_relations[i, j, 1] = 21
                    elif x == -15 and y == -15:
                        refined_relations[i, j, 1] = 42
                    elif x == -15 and y == 15:
                        refined_relations[i, j, 1] = 63
    np.save("refined_relations", refined_relations)


def save_rectangle_list():
    grid_coordinates = np.loadtxt("../config/open/open_grid_coordinates.csv", delimiter=",")
    neighbors = np.loadtxt("../config/open/open_neighbor.csv", int, delimiter=",")  # 邻点
    grid_list = []
    for i in range(len(grid_coordinates)):
        neighbor_down, neighbor_right = neighbors[i][1], neighbors[i][2]  # 下外
        if neighbor_down == -1 or neighbor_right == -1:
            continue
        low_right = neighbors[neighbor_right][1]
        zz, rr = grid_coordinates[i][0], grid_coordinates[i][1]  # 左上
        z1, r1 = grid_coordinates[neighbor_down][0], grid_coordinates[neighbor_down][1]  # 左下
        z2, r2 = grid_coordinates[neighbor_right][0], grid_coordinates[neighbor_right][1]  # 右上
        z3, r3 = grid_coordinates[low_right][0], grid_coordinates[low_right][1]  # 右下
        rectangle_i = Polygon([(r1, z1), (r3, z3), (r2, z2), (rr, zz)])  # 逆时针
        grid_list.append(rectangle_i)
    np.save("rectangle_list", grid_list)
    print("rectangle_list saved")


def cal_pre_next(arr, coordinate):
    index = np.searchsorted(arr, coordinate)
    lower, upper = -1.0, -1.0
    if index == 0:
        lower, upper = arr[index], arr[index + 1]
    elif index == len(arr):
        lower, upper = arr[index - 1], arr[index]
    else:
        lower, upper = arr[index - 1], arr[index]
    return lower, upper


def open_judge_relation():
    open_z = point_scatter(13, 2500, 6)
    open_z.insert(0, 7.5)
    open_z.insert(0, 3)
    open_z.insert(0, 0)
    open_r = point_scatter(300, 5000, 5, k=2)
    new_open_r = np.flipud(point_scatter(295, 0, 5, k=2, positive=False))
    r = np.concatenate((new_open_r, open_r))
    z = np.copy(open_z)
    # 因为只取到距离中心点400nm（0，400）处，且浓度体积为1立方微米，所以r:-500~500，z:-500~500
    xy_start = 0
    z_start = 0
    xy_end = 500
    z_end = 500
    # 每隔1nm取一个点
    interval = 1
    x_squared = np.square(np.arange(xy_start, xy_end + 1, interval))
    y_squared = np.square(np.arange(xy_start, xy_end + 1, interval))
    z_arange = np.arange(z_start, z_end + 1, interval)
    d1, d2, d3 = (xy_end - xy_start) // interval + 1, (xy_end - xy_start) // interval + 1, (
            z_end - z_start) // interval + 1  # 501 501 501
    grids_rz = np.full((d1, d2, d3, 4), dtype=float, fill_value=-1.0)
    points_rz = np.full((d1, d2, d3, 2), dtype=float, fill_value=-1.0)
    for i in range(d1):
        x2 = x_squared[i]
        print(f"{i}点开始")
        for j in range(d2):
            y2 = y_squared[j]
            radius = np.sqrt(x2 + y2)
            lower_r, upper_r = cal_pre_next(r, radius)
            for k in range(d3):
                height = z_arange[k]
                if height >= 7.5 and radius >= 300:
                    lower_z, upper_z = cal_pre_next(z, height)
                    grids_rz[i, j, k] = lower_r, upper_r, lower_z, upper_z
                    points_rz[i, j, k] = radius, height
    np.save("grids_rz_v2", grids_rz)
    np.save("points_rz_v2", points_rz)
    print("grid_rz saved")
    print("points_rz saved")


def main():
    # relations = nano_judge_relation()
    # np.save("relations", relations)
    # nano_relations_refine()
    open_judge_relation()


if __name__ == "__main__":
    main()
