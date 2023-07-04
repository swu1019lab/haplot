# -*- coding: utf-8 -*-
# @Time    : 2023/6/27 16:38
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : stats.py

from itertools import cycle
import numpy as np
import pandas as pd
from matplotlib import axes
import matplotlib.pyplot as plt


def boxplot(df: pd.DataFrame, by: str = 'column', ax: axes.Axes = None):
    if ax is None:
        ax = plt.gca()

    if by == "column":
        data = df.to_numpy()
        labels = df.columns.to_numpy()
    elif by == "index":
        data = df.T.to_numpy()
        labels = df.index.to_numpy()
    else:
        raise ValueError('by must be "column" or "index"')

    bp = ax.boxplot(data,
                    labels=labels,
                    widths=0.5,
                    showfliers=False,
                    # fill patch with color
                    patch_artist=True
                    )

    # Add in points to show each observation
    jitter = np.random.uniform(low=-0.25, high=0.25, size=data.shape)
    for i in range(len(bp.get('boxes'))):
        pc = ax.scatter([i + 1] * len(data[:, i]), data[:, i], s=2)
        pc.set_facecolor('gray')
        # set jitter for scatter
        pc.set_offsets(np.c_[jitter[:, i] + i + 1, data[:, i]])
        # move the scatters on top of the line.
        pc.set_zorder(2)

    # Set the axes ranges
    ax.set(
        xlim=(0.5, len(bp.get('boxes')) + 0.5),
        ylim=(None, None),
    )

    return ax


def ManhattanPlot(
        df: pd.DataFrame,
        chr_col: int = 0,
        pos_col: int = 1,
        p_value_col: int = 2,
        threshold0: float = None,
        threshold1: float = None,
        style: str = 'line',
        ax: axes.Axes = None
):
    """
    Manhattan plot

    The Manhattan plot is a type of scatter plot, usually used to display data
    with a large number of data-points

    Parameters
    ----------
    :param df: dataframe with columns of chromosome, position and p-value at least
    :param chr_col: column index of chromosome
    :param pos_col: column index of position
    :param p_value_col: column index of p-value
    :param threshold0: threshold for significant p-value
    :param threshold1: threshold for extremely significant p-value
    :param style: {'line', 'scatter', 'both'}, default 'line', style of plot
    :param ax: matplotlib axes object
    :return: matplotlib axes object
    """
    # transform xdata
    xdata = df.iloc[:, [chr_col, pos_col]]
    xdata.columns = ['chr', 'ps']
    xdata.set_index('chr', inplace=True)
    xdata = xdata.transform(lambda x: x * 1e-6)

    # transform ydata
    ydata = df.iloc[:, p_value_col]
    ydata = -np.log10(ydata.to_numpy())

    # chromosome colors, length and center loci
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = cycle(prop_cycle.by_key()['color'])
    chr_colors = xdata.groupby(level=0, group_keys=False).apply(lambda x: x.assign(colors=next(colors))).iloc[:, -1]
    chr_len = xdata.groupby(level=0).max()
    # add gap between chromosomes
    chr_len_new = chr_len.transform(lambda x: x + x.sum() / 100)
    chr_start = chr_len_new.cumsum() - chr_len_new
    chr_center = chr_start + chr_len_new / 2

    # transform x loci based on cum sum of chromosome length
    xdata_new = pd.concat([group + chr_start.loc[name] for name, group in xdata.groupby(level=0)])

    # plot
    if ax is None:
        ax = plt.gca()

    if style == 'line':
        ax.vlines(xdata_new.to_numpy().flatten(), 0, ydata, color=chr_colors, linewidth=0.5)
    elif style == 'scatter':
        ax.scatter(xdata_new.to_numpy().flatten(), ydata, s=1, c=chr_colors)
    elif style == 'both':
        ax.vlines(xdata_new.to_numpy().flatten(), 0, ydata, color=chr_colors, linewidth=0.5)
        ax.scatter(xdata_new.to_numpy().flatten(), ydata, s=1, c=chr_colors)
    else:
        raise ValueError('style must be "line", "scatter" or "both"')

    if threshold0 is not None:
        ax.axhline(-np.log10(threshold0), color='gray', linestyle='--', linewidth=0.5)
    if threshold1 is not None:
        ax.axhline(-np.log10(threshold1), color='gray', linewidth=0.5)

    ax.set_xlabel("Chromosome")
    ax.set_ylabel(r"$\mathrm{-log}_{10}(\mathit{p})$")

    ax.set_xlim(left=-1, right=np.max(xdata_new.values))
    ax.set_ylim(bottom=0, top=None)

    ax.set_xticks(chr_center.iloc[:, 0].to_numpy(), chr_center.index.to_list())
    ax.tick_params(axis='x', rotation=45)

    return ax
