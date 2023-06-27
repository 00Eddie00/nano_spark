import unittest
from decimal import Decimal as De
import numpy as np

grid = np.loadtxt("config/nano/4RYRgridt.dat", dtype=np.float64)
triangle_matrix = np.loadtxt("config/nano/4RYRnod.dat", dtype=int) - 1


def cal_cross_product(x, y, x1, y1, x2, y2):
    return (De(str(x1)) - De(str(x))) * (De(str(y2)) - De(str(y))) - (De(str(x2)) - De(str(x))) * (
            De(str(y1)) - De(str(y)))


def judge_relation(x, y, triangle_id):
    exclusive_relation = np.full(2, -1)  # 第一个值表示关系类型：1在顶点，2在边上，3在内部

    vertex1 = triangle_matrix[triangle_id, 0]
    vertex2 = triangle_matrix[triangle_id, 1]
    vertex3 = triangle_matrix[triangle_id, 2]
    x1 = grid[vertex1, 0]
    y1 = grid[vertex1, 1]
    x2 = grid[vertex2, 0]
    y2 = grid[vertex2, 1]
    x3 = grid[vertex3, 0]
    y3 = grid[vertex3, 1]
    # 判断是否为三角形单元的某顶点
    if x == x1 and y == y1:
        exclusive_relation[0] = 1
        exclusive_relation[1] = vertex1
        return exclusive_relation
    if x == x2 and y == y2:
        exclusive_relation[0] = 1
        exclusive_relation[1] = vertex2
        return exclusive_relation
    if x == x3 and y == y3:
        exclusive_relation[0] = 1
        exclusive_relation[1] = vertex3
        return exclusive_relation
    # 计算叉积
    c1 = cal_cross_product(x, y, x1, y1, x2, y2)
    c2 = cal_cross_product(x, y, x2, y2, x3, y3)
    c3 = cal_cross_product(x, y, x3, y3, x1, y1)
    # 判断是否在内部
    if (c1 > 0 and c2 > 0 and c3 > 0) or (c1 < 0 and c2 < 0 and c3 < 0):
        exclusive_relation[0] = 3
        exclusive_relation[1] = triangle_id
        return exclusive_relation
    # 判断是否在边界上
    if (c1 == 0 and c2 * c3 > 0) or (c2 == 0 and c3 * c1 > 0) or (c3 == 0 and c1 * c2 > 0):
        exclusive_relation[0] = 2
        exclusive_relation[1] = triangle_id
        return exclusive_relation
    # 在外部
    return exclusive_relation


def get_relations():
    relations = np.full((601, 601, 2), dtype=int, fill_value=-1)
    count = 1
    for i in range(601):
        print(f"x={i}")
        for j in range(601):
            x = i - 300
            y = j - 300
            if x ** 2 + y ** 2 <= 90000:
                # print("开始")
                # len(triangle_matrix)三角形个数
                for triangle_id in range(len(triangle_matrix)):
                    exclusive_relation = judge_relation(x, y, triangle_id)
                    if exclusive_relation[0] != -1:
                        relations[i, j] = exclusive_relation
                        break
            # print(count)
            count += 1
    return relations


class MyTestCase(unittest.TestCase):
    def test_something(self):
        x, y = 50, 50
        for triangle_id in range(len(triangle_matrix)):
            exclusive_relation = judge_relation(x, y, triangle_id)
            if exclusive_relation[0]!=-1:
                print(exclusive_relation)
                break



if __name__ == '__main__':
    unittest.main()
