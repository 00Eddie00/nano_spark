import unittest
import numpy as np
from generate_grid import generate_coordinates


class MyTestCase(unittest.TestCase):
    def test_something(self):
        a=np.zeros(5)
        a[0]=1
        a[1] = 1
        a[2] = 1
        a[3] = 1
        a[4] = 1
        print(a * 2 +5)


if __name__ == '__main__':
    unittest.main()
