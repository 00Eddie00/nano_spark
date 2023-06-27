import unittest
import numpy as np
from tool.cal_bcnl import *
from shapely.geometry import Point, Polygon, LineString, MultiPolygon
from decimal import Decimal


def gene_triangle():
    # 创建点坐标数组
    nods = np.loadtxt("config/nano/4RYRnod.dat", dtype=int) - 1
    points = np.loadtxt("config/nano/4RYRgridt.dat", dtype=np.float64)
    triangle_items = []
    for i in range(len(nods)):
        triangle1_vertices = []
        for j in range(3):
            nod_j = nods[i, j]
            x, y = points[nod_j, 0], points[nod_j, 1]
            triangle1_vertices.append((x, y))
        triangle_i = Polygon(triangle1_vertices)
        triangle_items.append(triangle_i)
    return triangle_items


def judge_relation():
    relations = np.full((601, 601, 2), dtype=int, fill_value=-1)
    # 创建点坐标数组
    nods = np.loadtxt("config/nano/4RYRnod.dat", dtype=int) - 1
    points = np.loadtxt("config/nano/4RYRgridt.dat", dtype=np.float64)
    triangle_items = gene_triangle()
    # 创建三角形集合对象
    print("三角形集合对象创建完成")
    for i in range(601):
        print(f"x={i}")
        for j in range(601):
            x = i - 300
            y = j - 300
            print(f"x={x}\ty={y}")
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
                        print(f"relations[{i}, {j}]={relations[i, j]}")
                        break
    return relations


class MyTestCase(unittest.TestCase):
    def test_something(self):
        x,y=3,4
        a=np.square(x)
        b = np.square(y)
        print(np.sqrt(a+b))


if __name__ == '__main__':
    unittest.main()
