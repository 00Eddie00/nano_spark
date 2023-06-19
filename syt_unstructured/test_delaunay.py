import unittest
import numpy as np
from shapely.geometry import Point, Polygon, MultiPolygon


class MyTestCase(unittest.TestCase):
    def test_something(self):
        relations = np.full((601, 601, 2), dtype=int, fill_value=-1)
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
        print("三角形集合对象创建完成")
        for i in range(601):
            for j in range(601):
                x = i - 300
                y = j - 300
                if x ** 2 + y ** 2 <= 90000:
                    print(f"x={x},y={y}")
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
                            # print("在三角形内部，包括边界（顶点）")
                            # print(f"relations[{i}, {j}, 0]={relations[i, j, 0]},relations[{i}, {j}, 1] ={relations[i, j, 1]}")
                            break


if __name__ == '__main__':
    unittest.main()
