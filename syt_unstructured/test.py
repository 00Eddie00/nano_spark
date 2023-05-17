import unittest
import numpy as np
from tool.cal_bcnl import cal_out_boundary_length


class MyTestCase(unittest.TestCase):
    def test_something(self):
        nano_avg=np.loadtxt("result/NANO/nano_avg/avg.dat")
        avg=nano_avg[:,:3]
        xx=np.arange(2001)
        print(avg[:,2])


if __name__ == '__main__':
    unittest.main()
