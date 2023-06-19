import numpy as np
from tool.cal_bcnl import *
from nano_spark.nano_parameters import *
from shapely.geometry import Point, Polygon, LineString, MultiPolygon
from decimal import Decimal


# 计算点到线段距离
def point_to_edge_distance(point, edge_start, edge_end):
    """
    计算一个点到边的距离

    参数:
        point: 待计算距离的点的坐标，格式为 (x, y)
        edge_start: 边的起点坐标，格式为 (x1, y1)
        edge_end: 边的终点坐标，格式为 (x2, y2)

    返回值:
        distance: 点到边的距离
    """
    p = Point(point)
    line = LineString([edge_start, edge_end])

    distance = p.distance(line)
    return distance


#  计算三角形面积
def triangle_area(x1, y1, x2, y2, x3, y3):
    # 构造三角形的三个顶点
    a = Point(x1, y1)
    b = Point(x2, y2)
    c = Point(x3, y3)
    # 构造三角形对象
    triangle = Polygon([(a.x, a.y), (b.x, b.y), (c.x, c.y)])
    # 计算三角形的面积
    area = triangle.area
    return area


def judge_relation():
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
        print(f"x={i}")
        for j in range(601):
            x = i - 300
            y = j - 300
            if x ** 2 + y ** 2 <= 90000:
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


# def judge_relation(x, y, triangle_id):
#     grid = np.loadtxt(nano_grid_file_name, dtype=np.float64)
#     nod = np.loadtxt(nod_file_name, dtype=int) - 1
#     exclusive_relation = np.full(2, -1)  # 第一个值表示关系类型：1在顶点，2在边上，3在内部
#     single_area, control_area, near_triangle, index_in_triangle, nix_multiply_l, niy_multiply_l, b_arr, c_arr, nmax, total_area = cal_elements(
#         nano_grid_file_name, nod_file_name)
#     for i in nod[triangle_id]:
#         if x == grid[i, 0] and y == grid[i, 1]:
#             # 在顶点上
#             exclusive_relation = [1, i]
#             return exclusive_relation
#     nod1, nod2, nod3 = nod[triangle_id, 0], nod[triangle_id, 1], nod[triangle_id, 2]
#     x1, y1 = grid[nod1, 0], grid[nod1, 1]
#     x2, y2 = grid[nod2, 0], grid[nod2, 1]
#     x3, y3 = grid[nod3, 0], grid[nod3, 1]
#     area = single_area[triangle_id]
#     area1 = triangle_area(x, y, x2, y2, x3, y3)
#     area2 = triangle_area(x1, y1, x, y, x3, y3)
#     area3 = triangle_area(x1, y1, x2, y2, x, y)
#     distance1 = point_to_edge_distance((x, y), (x1, y1), (x2, y2))
#     distance2 = point_to_edge_distance((x, y), (x1, y1), (x3, y3))
#     distance3 = point_to_edge_distance((x, y), (x2, y2), (x3, y3))
#     if area - (area1 + area2 + area3) < 0:
#         #  在外部
#         return exclusive_relation
#     elif distance1 == 0.0 or distance2 == 0.0 or distance3 == 0.0:
#         # 在边界
#         exclusive_relation = [2, triangle_id]
#         return exclusive_relation
#     else:
#         # 在内部
#         exclusive_relation = [3, triangle_id]
#         return exclusive_relation


def get_relations():
    relations = judge_relation()
    print("修改一些点")
    refined_relations = np.copy(relations)
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


def main():
    get_relations()


if __name__ == "__main__":
    main()
