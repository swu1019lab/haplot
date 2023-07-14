# -*- coding: utf-8 -*-
# @Time    : 2023/7/2 22:03
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : data.py

import numpy as np
import pandas as pd


def generate_data(sort=False, chr_num=3, snp_num=1000, gene_num=2):
    p_value1 = np.random.uniform(0, 1, snp_num)
    p_value2 = np.random.uniform(0, 1, snp_num)
    data = pd.DataFrame(
        data={
            'chr': np.random.choice(['chr{}'.format(i + 1) for i in range(chr_num)], snp_num),
            'pos': np.random.randint(1, 100000000, snp_num),
            'p_value1': p_value1,
            'p_value2': p_value2,
            # np.nan should be with numeric type
            'gene1': pd.Series(['gene{}'.format(i + 1) for i in range(gene_num)],
                               index=np.argsort(p_value1)[:gene_num]),
            'gene2': pd.Series(['gene{}'.format(i + 1) for i in range(gene_num)],
                               index=np.argsort(p_value2)[:gene_num])
        },
        index=np.arange(snp_num)
    )
    if sort:
        data.sort_values(by=['chr', 'pos'], inplace=True, ignore_index=True)
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
