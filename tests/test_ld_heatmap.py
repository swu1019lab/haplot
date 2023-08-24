# -*- coding: utf-8 -*-
# @Time    : 2023/7/14 20:55
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : test_ld_heatmap.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from haplot.chart import LDHeatmapPlot

np.random.seed(123)
# generate random R^2 value for ld heatmap plot
# snp id and snp position
arrays = [['rs{}'.format(i + 1) for i in range(10)],
          np.sort(np.random.randint(1, 1000000, 10))]
index = pd.MultiIndex.from_arrays(arrays, names=('rs_id', 'pos'))

df = pd.DataFrame(
    data=np.random.uniform(0, 1, (10, 10)),
    index=index,
    columns=index
)

fig, ax = plt.subplots(figsize=(5, 5))
LDHeatmapPlot(df, plot_value=True, plot_snp=True, plot_diag=False, cmap='Reds')
plt.show()
