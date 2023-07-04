# -*- coding: utf-8 -*-
# @Time    : 2023/7/3 16:44
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : test1.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from haplot.stats import ManhattanPlot
from haplot.theme import Theme


# generate random data for manhattan plot
def generate_data(sort=False):
    data = pd.DataFrame(
        {
            'chr': np.random.choice(['chr1', 'chr2', 'chr3'], 1000),
            'pos': np.random.randint(1, 100000000, 1000),
            'p_value': np.random.uniform(0, 1, 1000)
        }
    )
    if sort:
        data.sort_values(by=['chr', 'pos'], inplace=True)
    return data


np.random.seed(123)
fig, ax = plt.subplots(figsize=(10, 5))
df = generate_data(sort=True)
ManhattanPlot(df)
Theme.apply(ax)
plt.show()
