# -*- coding: utf-8 -*-
# @Time    : 2023/7/17 15:41
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : test_geomap.py

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from haplot.chart import GeoMapPlot


# generate data
np.random.seed(0)
x = np.random.randint(100, 115, 5)
y = np.random.randint(20, 40, 5)
value1 = np.random.randint(1, 100, 5)
value2 = np.random.randint(1, 100, 5)
value3 = np.random.randint(1, 100, 5)
data = pd.DataFrame({
    'x': x, 'y': y,
    'value1': value1,
    'value2': value2,
    'value3': value3
})

map_file = "china_data.json"
bound_file = "bound_data.geojson"

fig, ax = plt.subplots(figsize=(5, 5))
GeoMapPlot(data, value_col=[2, 3, 4])
plt.show()
