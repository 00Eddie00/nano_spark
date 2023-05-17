import numpy as np


def generate_nano():
    nano_r = []
    i = 0.0
    while i < 300:
        nano_r.append(i)
        if i < 15:
            i = i + 0.5
        elif i < 50:
            i = i + 1.0
        elif i < 100:
            i = i + 2.0
        else:
            i = i + 4.0
    j = 7.5
    nano_z = []
    while j >= -7.5:
        nano_z.append(j)
        if j > 0.5:
            j = j - 0.5
        else:
            if j == 0.5:
                nano_z.append(0.0)
            j = j - 1.0
    nano_z.reverse()
    nano_coordinates = [(x, y) for y in nano_r for x in nano_z]
    return nano_z, nano_r, nano_coordinates


def generate_open():
    n = 0
    open_z1 = []
    for i in range(5, 72):
        open_z1.append(n)
        n = n + i
    open_z2 = [-x for x in open_z1]
    open_z2.remove(0)
    open_z2.reverse()
    open_z = open_z2 + open_z1
    n = 0
    open_r = []
    for i in range(5, 97):
        n = n + i
        open_r.append(n + 300)
    open_coordinates = [(x, y) for y in open_r for x in open_z]
    return open_z, open_r, open_coordinates


# 新开放区域
def generate_new_open():
    n = 0
    new_open_z1 = []
    for i in range(5, 71):
        n = n + i
        new_open_z1.append(n)
    new_open_z1.remove(5)
    new_open_z2 = [-x for x in new_open_z1]
    new_open_z2.reverse()

    i = 0.0
    new_open_r = []
    while i < 300:
        new_open_r.append(i)
        if i < 200:
            i = i + 20
        elif i < 250:
            i = i + 10
        else:
            i = i + 5
    new_open_coordinates1 = [(x, y) for y in new_open_r for x in new_open_z1]
    new_open_coordinates2 = [(x, y) for y in new_open_r for x in new_open_z2]
    return new_open_coordinates1, new_open_coordinates2


def generate_coordinates():
    nano_z, nano_r, nano_coordinates = generate_nano()
    open_z, open_r, open_coordinates = generate_open()
    boundary_z = list(set(open_z + nano_z))
    # print(len(boundary_z))
    boundary_z.sort()
    boundary_coordinates = [(x, y) for y in [300] for x in boundary_z]
    new_open_coordinates1, new_open_coordinates2 = generate_new_open()

    nano_length = len(nano_coordinates)
    boundary_length = len(boundary_coordinates)
    open_length = len(open_coordinates)
    return nano_coordinates, open_coordinates, boundary_coordinates, new_open_coordinates1, new_open_coordinates2, nano_length, boundary_length, open_length


if __name__ == '__main__':
    nano_coordinates, open_coordinates, boundary_coordinates, new_open_coordinates1, new_open_coordinates2, nano_length, boundary_length, open_length = generate_coordinates()
    grid_coordinates = nano_coordinates + boundary_coordinates + open_coordinates
    np.savetxt("config/grid/grid_coordinates.csv", grid_coordinates, fmt='%.1f', delimiter=",")
