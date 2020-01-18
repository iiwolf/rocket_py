import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import time
import unittest

from constants import *
from rocket_calc import run_rocket_config
# runs various test cases given by professor frederick to check correctness

class TestRocketConfigs(unittest.TestCase):
    def test1(self):
        self.assertAlmostEqual( run_rocket_config(1.0 * units["in"], 2.375 * units["in"], 8.0 * units["in"],
        4, 4, 1.0 * units["in**2"], 1.0 * units["lb"]),57570.44, places=2)

    def test2(self):
        self.assertAlmostEqual( run_rocket_config(1.0 * units["in"], 2.000 * units["in"], 8.0 * units["in"],
        4.0, 4.0, 1.0 * units["in**2"], 1.0 * units["lb"]),33278.44, places=2)

    def test3(self):
        self.assertAlmostEqual( run_rocket_config(1.0 * units["in"], 2.000 * units["in"], 1.0 * units["in"],
        10.0, 4.0, 1.0 * units["in**2"], 1.0 * units["lb"]),10948.62, places=2)

    def test4(self):
        self.assertAlmostEqual( run_rocket_config(1.0 * units["in"], 2.375 * units["in"], 1.0 * units["in"],
        2.0, 4.0, 1.0 * units["in**2"], 0.0 * units["lb"]),426.61, places=2)

    def test5(self):
        self.assertAlmostEqual( run_rocket_config(1.0 * units["in"], 1.5 * units["in"], 4.0 * units["in"],
        4.0, 4.0, 1.0 * units["in**2"], 0.887 * units["lb"]), 5000.00, places=2)

    def test_5000(self):
        self.assertAlmostEqual( run_rocket_config(0.001 * units["in"], 1.1968 * units["in"], 1.575 * units["in"],
        11.0, 2.666666, 4.0 * units["in**2"], 0.189105 * units["lb"]), 5000.00, places=2)

    def test_10000(self):
        self.assertAlmostEqual( run_rocket_config(1.0 * units["in"], 1.5 * units["in"], 4.0 * units["in"],
        4.0, 4.0, 1.0 * units["in**2"], 0.887 * units["lb"]), 5000.00, places=2)

    def test_15000(self):
        self.assertAlmostEqual( run_rocket_config(1.0 * units["in"], 1.5 * units["in"], 4.0 * units["in"],
        4.0, 4.0, 1.0 * units["in**2"], 0.887 * units["lb"]), 5000.00, places=2)



unittest.main() # OK
# def run_rocket_tests():
#     alt = run_rocket_config(1.0 * units["in"], 2.375 * units["in"], 8.0 * units["in"],
#         4 * units[""], 4, 1.0 * units["in**2"], 1.0 * units["lb"])
    
#     print(alt)
#     assertAlmostEqual(alt, 57570.44,2)

# run_rocket_tests