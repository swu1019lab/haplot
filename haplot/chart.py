# -*- coding: utf-8 -*-
# @Time    : 2023/6/27 16:38
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : chart.py

from itertools import cycle
import numpy as np
import pandas as pd
from scipy import stats
from matplotlib import axes
import matplotlib.pyplot as plt


def boxplot(df: pd.DataFrame, by: str = 'column', ax: axes.Axes = None):
    """
    Boxplot

    a helpful function to plot boxplot

    Parameters
    ----------
    :param df: a DataFrame object to plot
    :param by: 'column' or 'index'
    :param ax: a matplotlib.axes.Axes object
    :return: a matplotlib.axes.Axes object
    """
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
        value_col: int | list = 2,
        threshold0: float = None,
        threshold1: float = None,
        chr_names: list | str = None,
        style: str = 'line',
        ax: axes.Axes = None):
    """
    Manhattan plot

    The Manhattan plot is a type of scatter plot, usually used to display data
    with a large number of data-points

    Parameters
    ----------
    :param df: dataframe with columns of chromosome, position and p-value at least
    :param chr_col: column index of chromosome
    :param pos_col: column index of position
    :param value_col: column index of value, or list of column index of value
    :param threshold0: threshold for significant p-value
    :param threshold1: threshold for extremely significant p-value
    :param chr_names: list of chromosome names selected to plot
    :param style: {'line', 'scatter'}, default 'line', style of plot
    :param ax: matplotlib axes object
    :return: matplotlib axes object
    """
    # filter data
    if chr_names is not None:
        if isinstance(chr_names, str):
            chr_names = [chr_names]
        df = df.loc[df.iloc[:, chr_col].isin(chr_names), :]

    # transform xdata
    xdata = df.iloc[:, [chr_col, pos_col]]
    xdata.columns = ['chr', 'ps']
    xdata.set_index('chr', inplace=True)
    xdata = xdata.transform(lambda x: x * 1e-6)

    # transform ydata
    ydata = df.iloc[:, value_col]
    ydata = -np.log10(ydata.to_numpy())

    # chromosome colors, length and center loci
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = cycle(prop_cycle.by_key()['color'])
    cdata = xdata.groupby(level=0, group_keys=False).apply(lambda x: x.assign(colors=next(colors))).iloc[:, -1]
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
        ax.vlines(xdata_new.to_numpy().flatten(), 0, ydata, color=cdata, linewidth=0.5)
    elif style == 'scatter':
        ax.scatter(xdata_new.to_numpy().flatten(), ydata, s=1, c=cdata)
    else:
        raise ValueError('style must be "line", "scatter"')

    if threshold0 is not None:
        ax.axhline(-np.log10(threshold0), color='gray', linestyle='--', linewidth=0.5)
    if threshold1 is not None:
        ax.axhline(-np.log10(threshold1), color='gray', linewidth=0.5)

    ax.set_xlim(left=-1, right=np.max(xdata_new.values))
    ax.set_xticks(chr_center.iloc[:, 0].to_numpy(), chr_center.index.to_list())

    return ax


def QQPlot(
        df: pd.DataFrame,
        p_value_col: int = 2,
        ax: axes.Axes = None):
    if ax is None:
        ax = plt.gca()

    p = df.iloc[:, p_value_col]
    # return two tuple: (osm, osr), (slope, intercept, r)
    # 1st tuple contains two arrays:
    # first is theoretical quantiles, second is sample quantiles
    res = stats.probplot(p, dist=stats.uniform)
    osm, osr = res[0]
    ax.scatter(-np.log10(osm), -np.log10(osr))
    # slope, intercept, r = res[1]
    reg = stats.linregress(-np.log10(osm), -np.log10(osr))
    ax.plot(-np.log10(osm), -np.log10(osm) * reg[0] + reg[1], 'r--')
    ax.text(0.05, 0.95, f"R^2 = {reg[2] ** 2:.3f}", transform=ax.transAxes)
    return ax


def TablePlot(
        df: pd.DataFrame,
        ax: axes.Axes = None):
    if ax is None:
        ax = plt.gca()

    ax.axis('off')
    ax.axis('tight')
    ax.table(cellText=df.values, colLabels=df.columns, loc='center')

    return ax


def BarPlot(
        df: pd.DataFrame,
        x_col: int = 0,
        y_col: int | list = 1,
        ax: axes.Axes = None):
    if ax is None:
        ax = plt.gca()

    x_data = df.iloc[:, x_col].to_numpy()
    y_data = df.iloc[:, y_col].to_numpy()

    if np.unique(x_data).size != x_data.size:
        raise ValueError('The value of x_col must be unique')

    if isinstance(y_col, int):
        ax.bar(x_data, y_data)
        return ax

    for i in range(len(y_col)):
        ax.bar(x_data, y_data[:, i], bottom=y_data[:, :i].sum(axis=1))

    return ax
