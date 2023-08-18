import numpy as np


def generate_interval(dis, half_length):
    # 每隔1nm取一个点
    interval_arr = [5, 10, 20, 25, 50]
    end_value = 0.0
    interval = half_length // len(interval_arr)
    process_points = np.empty(0)
    for i in range(len(interval_arr)):
        start_value = end_value
        end_value = end_value + interval
        num_points = interval // interval_arr[i] + 1
        data_array = np.linspace(start_value, end_value, int(num_points))
        if i != len(interval_arr) - 1:
            data_array = data_array[:-1]
        process_points = np.concatenate((process_points, data_array))

    x_half_left = dis - np.flip(process_points)
    x_half_right = process_points + dis
    x_arr = np.concatenate((x_half_left[:-1], x_half_right))

    y_half_left = - np.flip(process_points)
    y_half_right = process_points
    y_arr = np.concatenate((y_half_left[:-1], y_half_right))

    z_half_left = - np.flip(process_points)
    z_half_right = process_points
    z_arr = np.concatenate((z_half_left[:-1], z_half_right))
    return x_arr, y_arr, z_arr, process_points


def cal_nano_radius(position_list):
    position_list_len = len(position_list)  # 需要检测的点个数
    for position_index in range(position_list_len):
        radius_list = []
        dis = position_list[position_index][0]
        x_arr, y_arr, z_arr, process_points = generate_interval(dis, 1000 / 2)
        x_square, y_square = np.square(x_arr), np.square(y_arr)
        for i in range(len(x_arr)):
            x2 = x_square[i]
            for j in range(len(y_arr)):
                y2 = y_square[j]
                radius = np.sqrt(x2 + y2)
                if radius <= 300:
                    radius_list.append(radius)
        list = np.array(radius_list)
        list2 = np.unique(list)
        np.save(
            f"../optical_blurring/xy_list/({position_list[position_index][0]},{position_list[position_index][1]})_radius_list",
            list2)


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


def cal_nano_points():
    # open_r = point_scatter(300, 5000, 5, k=2)
    new_open_r = np.flipud(point_scatter(295, 0, 5, k=2, positive=False))
    radius_list = np.concatenate((new_open_r, (300,)))
    xy_index = np.full((len(radius_list), 50, 2), fill_value=-1)
    xy_len_index = np.zeros(len(radius_list), dtype=int)
    for i in range(601):
        x = i - 300
        x2 = np.square(x)
        for j in range(601):
            y = j - 300
            y2 = np.square(y)
            r = np.sqrt(x2 + y2)
            indexes = np.where(radius_list == r)
            if len(indexes[0]) > 0:
                a = indexes[0][0]
                b = xy_len_index[a]
                xy_index[a][b][0], xy_index[a][b][1] = i, j
                xy_len_index[a] = xy_len_index[a] + 1
    np.save(
        f"../optical_blurring/xy_list/xy_index", xy_index)


def main():
    # my_position_list = [[0, 0], [100, 0], [300, 0], [400, 0]]
    # cal_nano_radius(my_position_list)
    new_open_r = np.flipud(point_scatter(295, 0, 5, k=2, positive=False))
    radius_list = np.concatenate((new_open_r, (300,)))
    print(radius_list)
    # radius_list = np.concatenate((new_open_r, open_r))[:17]
    # indexes = np.where(radius_list == 295.0)
    # xy_index = np.load("../optical_blurring/xy_list/xy_index.npy")
    # a=indexes[0][0]
    # print(xy_index[a])
    # cal_nano_points()


if __name__ == "__main__":
    main()
