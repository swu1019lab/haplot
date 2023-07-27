# -*- coding: utf-8 -*-
# @Time    : 2023/7/5 21:44
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : test3.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from haplot.data import generate_table_data


np.random.seed(123)
df = generate_table_data(row=10)
fig, ax = plt.subplots(figsize=(10, 5))
# TablePlot(df)
df.table()
plt.show()
fig.savefig('test3.png')
