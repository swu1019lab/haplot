# -*- coding: utf-8 -*-
# @Time    : 2023/7/2 22:03
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : data.py

import numpy as np
import pandas as pd


def generate_data(sort=False):
    data = pd.DataFrame(
        {
            'chr': np.random.choice(['chr1', 'chr2', 'chr3'], 1000),
            'pos': np.random.randint(1, 100000000, 1000),
            'p_value1': np.random.uniform(0, 1, 1000),
            'p_value2': np.random.uniform(0, 1, 1000),
        }
    )
    if sort:
        data.sort_values(by=['chr', 'pos'], inplace=True)
    return data


def generate_table_data(row=10):
    data = pd.DataFrame(
        {
            'chr': np.random.choice(['chr1', 'chr2', 'chr3'], row),
            'pos': np.random.randint(1, 100000000, row),
            'p_value': np.random.uniform(0, 1, row),
            'rs_id': np.random.choice(['rs1', 'rs2', 'rs3'], row),
            'gene': np.random.choice(['gene1', 'gene2', 'gene3'], row),
            'beta': np.random.uniform(-1, 1, row),
            'se': np.random.uniform(0, 1, row),
            'maf': np.random.uniform(0, 1, row),
        }
    )
    return data
