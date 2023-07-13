from nano_spark.open.open_parameters import *

def point_scatter(start, end, base, k=1, positive=True):
    result = [start]
    flag = 1 if positive else -1
    i = 0
    while True:
        step = (k * i + base) * flag
        value = result[i] + step
        if positive is True and value > end:
            result[-1] = end
            break
        elif positive is False and value < end:
            result[-1] = end
            break
        else:
            result.append(value)
            i = i + 1
    return result


def meshing():
    # 二维结构化网格划分 第一列为z轴 第二列为ρ轴，或称为r轴
    # 大开放区域
    open_z = point_scatter(13, 2500, 6)
    open_z.insert(0, 7.5)
    open_z.insert(0, 3)
    open_z.insert(0, 0)
    open_r = point_scatter(300, 5000, 5, k=2)
    open_coordinates = np.concatenate(
        (np.tile(open_z, len(open_r))[:, np.newaxis], np.repeat(open_r, len(open_z))[:, np.newaxis]), axis=1)
    # 小开放区域
    new_open_z = point_scatter(13, 2500, 6)
    new_open_r = np.flipud(point_scatter(295, 0, 5, k=2, positive=False))
    new_open_coordinates = np.concatenate(
        (np.tile(new_open_z, len(new_open_r))[:, np.newaxis], np.repeat(new_open_r, len(new_open_z))[:, np.newaxis]),
        axis=1)
    # 得到整个网格的点坐标
    grid_coordinates = np.concatenate((open_coordinates, new_open_coordinates))
    np.savetxt(open_grid_file_name, grid_coordinates, fmt='%.1f', delimiter=",")


if __name__ == "__main__":
    meshing()
