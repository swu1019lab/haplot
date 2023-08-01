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


def generate_gene_data(
        gene_num=3,
        exon_min_num=1,
        exon_max_num=10,
        exon_min_len=100,
        exon_max_len=1000,
        intron_min_len=100,
        intron_max_len=1000,
        utr_min_len=0,
        utr_max_len=1000,
):
    features_data_list = []
    for i in range(gene_num):
        offset = np.random.randint(1, 100000)
        exon_num = np.random.randint(exon_min_num, exon_max_num)
        exon_len = np.random.randint(exon_min_len, exon_max_len, exon_num)
        intron_len = np.random.randint(intron_min_len, intron_max_len, exon_num - 1)
        utr_len = np.random.randint(utr_min_len, utr_max_len, 2)

        features = np.zeros(exon_num * 2 - 1, dtype=int)
        features[::2] = exon_len
        features[1::2] = intron_len
        features = np.insert(features, 0, utr_len[0])
        features = np.append(features, utr_len[1])
        features_end = np.cumsum(features)
        features_start = features_end - features

        features_type = ['exon' if i % 2 == 0 else 'intron' for i in range(exon_num * 2 - 1)]
        features_type = np.insert(features_type, 0, '5UTR')
        features_type = np.append(features_type, '3UTR')

        features_strand = np.repeat(np.random.choice(['+', '-'], 1), exon_num * 2 + 1)
        features_name = np.repeat('gene{}'.format(i + 1), exon_num * 2 + 1)

        features_data_list.append(
            pd.DataFrame(
                data={
                    'name': features_name,
                    'start': features_start + offset,
                    'end': features_end + offset,
                    'strand': features_strand,
                    'type': features_type,
                }
            )
        )
    return pd.concat(features_data_list, ignore_index=True)


df = generate_gene_data()

fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(111)
GeneStrucPlot(df)
plt.show()
