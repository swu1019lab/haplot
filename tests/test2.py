# -*- coding: utf-8 -*-
# @Time    : 2023/7/4 10:30
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : test2.py


import numpy as np
import matplotlib.pyplot as plt
from haplot.data import generate_data
from haplot.chart import QQPlot
from haplot.theme import Theme


np.random.seed(123)
df = generate_data(sort=True)
fig, ax = plt.subplots(figsize=(5, 5))
QQPlot(df, 2)
Theme.apply()
plt.show()
