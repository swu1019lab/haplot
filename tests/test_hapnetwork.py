# -*- coding: utf-8 -*-
# @Time    : 2023/9/2 10:30
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : test_hapnetwork.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from haplot.chart import HapNetworkPlot


# generate data for network plot
np.random.seed(20230831)
node_num = 10
edge_num = 20
node_name = ['Hap{}'.format(i + 1) for i in range(node_num)]
node_data = pd.DataFrame(
    {
        'node': node_name,
        'value1': np.random.randint(1, 400, node_num),
        'value2': np.random.randint(1, 400, node_num),
        'value3': np.random.randint(1, 400, node_num)
    }
)
edge_data = pd.DataFrame(
    {
        'source': np.random.choice(node_name, edge_num),
        'target': np.random.choice(node_name, edge_num),
        'weight': np.random.randint(1, 10, edge_num),
        'color': ['C2'] * edge_num
    }
)

# plot
fig, ax = plt.subplots(figsize=(5, 5))
HapNetworkPlot(edge_data=edge_data, node_data=node_data, layout='spring', node_font_size=8)
plt.show()
