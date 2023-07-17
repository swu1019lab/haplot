# -*- coding: utf-8 -*-
# @Time    : 2023/7/14 20:55
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : test5.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from haplot.chart import LDHeatmapPlot

np.random.seed(123)
# generate random R^2 value for ld heatmap plot
df = pd.DataFrame(
    data=np.random.uniform(0, 1, (10, 10)),
    index=['rs{}'.format(i + 1) for i in range(10)],
    columns=['rs{}'.format(i + 1) for i in range(10)]
)

fig, ax = plt.subplots(figsize=(5, 5))
LDHeatmapPlot(df, plot_value=True, cmap='Reds')
plt.show()
