import unittest
import numpy as np
import os

class MyTestCase(unittest.TestCase):
    def test_something(self):
        dirnames = ["Ca", "CaB", "CaF","CaB"]
        path = "result/"
        nano = "NANO/"
        open = "OPEN/"
        filenames1 = os.listdir("result/NANO/CaB/")  # 目前所有步数的浓度文件
        print(filenames1[-1])
        print(len(filenames1))



if __name__ == '__main__':
    unittest.main()
