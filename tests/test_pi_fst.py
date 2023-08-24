# -*- coding: utf-8 -*-
# @Time    : 2023/7/19 11:15
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : test_pi_fst.py

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from haplot.chart import PiWithFstPlot


# generate random data
np.random.seed(20230719)
x = np.random.randn(1000)
y = np.random.randn(1000)
df = pd.DataFrame({'x': x, 'y': y})
fig = plt.figure(figsize=(8, 8))
PiWithFstPlot(df)
plt.show()
