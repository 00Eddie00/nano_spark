import unittest
import numpy as np
import os
from scipy.spatial import Delaunay
from shapely.geometry import Point, LineString, Polygon, MultiPolygon
from math import exp, sqrt, log


# 高斯函数
def gaussian_function(w, sigma):
    result = exp(-(w ** 2) / (2 * sigma)) / sqrt(2 * np.pi * sigma)
    return result


def gaussian_function2(w, sigma):
    return np.power((2 * np.pi * sigma), -1 / 2) * np.exp(-np.square(w) / (2 * sigma))


class MyTestCase(unittest.TestCase):
    def test_something(self):
       for i in range(0,30,10):
           print(i)


if __name__ == '__main__':
    unittest.main()
