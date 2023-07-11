# -*- coding: utf-8 -*-
# @Time    : 2023/7/5 22:23
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : test4.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from haplot.chart import BarPlot
from haplot.theme import Theme


np.random.seed(123)
df = pd.DataFrame(
    {
        'x': ['A', 'B', 'C', 'D', 'E'],
        'y1': np.random.randint(1, 10, 5),
        'y2': np.random.randint(1, 10, 5),
        'y3': np.random.randint(1, 10, 5),
    }
)
fig, ax = plt.subplots(figsize=(5, 5))
BarPlot(df, x_col=0, y_col=[1, 2, 3])
Theme.apply()
plt.show()
fig.savefig('test4.png')
