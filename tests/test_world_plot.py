# -*- coding: utf-8 -*-
# @Time    : 2023/8/23 20:33
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : test_world_plot.py

import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from haplot.chart import GeoMapPlot

# 查询世界地图上的城市坐标
# new york
ny_lon, ny_lat = -75, 43
# london
ld_lon, ld_lat = 0, 51.5
# beijing
bj_lon, bj_lat = 116.4, 39.9
# tokyo
tk_lon, tk_lat = 139.7, 35.7
# moscow
ms_lon, ms_lat = 37.6, 55.7

np.random.seed(20230823)
data = pd.DataFrame({
    'x': [ny_lon, ld_lon, bj_lon, tk_lon, ms_lon],
    'y': [ny_lat, ld_lat, bj_lat, tk_lat, ms_lat],
    'value1': np.random.randint(1, 100, 5),
    'value2': np.random.randint(1, 100, 5),
    'value3': np.random.randint(1, 100, 5)
})

fig, ax = plt.subplots(figsize=(10, 5), subplot_kw={'projection': ccrs.PlateCarree()})
ax.stock_img()
GeoMapPlot(data, value_col=[2, 3, 4])
plt.show()
