# -*- coding: utf-8 -*-
# @Time    : 2023/7/2 22:03
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : data.py

import numpy as np
import pandas as pd


def generate_associated_data(sort=False, chr_num=3, snp_num=1000, gene_num=2):
    p_value1 = np.random.uniform(0, 1, snp_num)
    p_value2 = np.random.uniform(0, 1, snp_num)
    r2_value1 = np.random.uniform(0, 1, snp_num)
    r2_value2 = np.random.uniform(0, 1, snp_num)
    data = pd.DataFrame(
        data={
            'chr': np.random.choice(['chr{}'.format(i + 1) for i in range(chr_num)], snp_num),
            'pos': np.random.randint(1, 100000000, snp_num),
            'p_value1': p_value1,
            'r2_value1': r2_value1,
            'p_value2': p_value2,
            'r2_value2': r2_value2,
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
