# -*- coding: utf-8 -*-
# @Time    : 2023/7/22 20:49
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : test8.py

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from haplot.chart import GeneStrucPlot

# generate random data of gene structure
np.random.seed(20230722)
gene = np.repeat(['gene1', 'gene2', 'gene3'], 5)
start = np.cumsum(np.random.randint(1, 1000, 15))
end = start + np.random.randint(1, 10, 15)
strand = np.repeat(np.random.choice(['+', '-'], 3), 5)
feature = np.tile(['exon', 'intron', '5UTR'], 5)

df = pd.DataFrame({'gene': gene, 'start': start, 'end': end, 'strand': strand, 'feature': feature})

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
GeneStrucPlot(df, ax=ax)
plt.show()
