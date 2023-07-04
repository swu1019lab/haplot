# -*- coding: utf-8 -*-
# @Time    : 2023/6/28 11:44
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : test.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from haplot.theme import BoxTheme
from haplot.stats import boxplot


np.random.seed(19680801)
all_data = [np.random.normal(0, std, 100) for std in range(1, 6)]
labels = ['x1', 'x2', 'x3', 'x4', 'x5']
df = pd.DataFrame(all_data, index=labels).T

fig, ax = plt.subplots()
boxplot(df)
BoxTheme.apply(ax)
plt.show()
fig.savefig('test.png')
