# -*- coding: utf-8 -*-
# @Time    : 2023/7/22 20:49
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : test_gene_struc.py

import numpy as np
import matplotlib.pyplot as plt
from haplot.data import generate_gene_data
from haplot.chart import GeneStrucPlot

# generate random data of gene structure
np.random.seed(20230722)

df = generate_gene_data(gene_num=1)
fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(111)
GeneStrucPlot(df)
plt.show()
