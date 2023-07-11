# -*- coding: utf-8 -*-
# @Time    : 2023/7/3 16:44
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : test1.py

import numpy as np
import matplotlib.pyplot as plt
from haplot.chart import ManhattanPlot
from haplot.theme import ManhattanTheme
from haplot.data import generate_data


# generate random data for manhattan plot
df = generate_data(sort=True)

np.random.seed(123)
fig, ax = plt.subplots(figsize=(5, 5))
ManhattanPlot(df, value_col=[2, 3])
ManhattanTheme.apply(ax)
plt.show()
