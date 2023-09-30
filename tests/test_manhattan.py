# -*- coding: utf-8 -*-
# @Time    : 2023/7/3 16:44
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : test_manhattan.py

import numpy as np
import matplotlib.pyplot as plt
from haplot.chart import ManhattanPlot
from haplot.data import generate_data
from haplot.theme import ManhattanTheme

np.random.seed(123)
# generate random data for manhattan plot
df = generate_data(sort=True, snp_num=1000, gene_num=2)
# print(df)
fig, ax = plt.subplots(figsize=(12, 4))
axs = ManhattanPlot(df, value_col=[2, 3], value_layout='stack', style='scatter', log=True)
for a in axs:
    ManhattanTheme.apply(a)
plt.tight_layout()
plt.show()
