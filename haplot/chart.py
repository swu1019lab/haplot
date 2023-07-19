# -*- coding: utf-8 -*-
# @Time    : 2023/6/27 16:38
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : chart.py

from itertools import cycle
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib as mpl
from matplotlib import figure, axes, ticker
from matplotlib.patches import RegularPolygon, Patch
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.text import OffsetFrom
import matplotlib.pyplot as plt
import geopandas as gpd


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
        ann_col: int = None,
        threshold0: float = None,
        threshold1: float = None,
        chr_names: list | str = None,
        log: bool = False,
        colors: list = None,
        style: str = 'line',
        fig: figure.Figure = None):
    """
    Manhattan plot

    Features
    ----------
    Support line and scatter style \n
    Support multiple value columns \n
    Support multiple or selected chromosomes to plot \n
    Support annotation \n

    Parameters
    ----------
    :param df: dataframe with columns of chromosome, position and value at least
    :param chr_col: column index of chromosome
    :param pos_col: column index of position
    :param value_col: column index of value or list of column index of value
    :param ann_col: column index of annotation
    :param threshold0: threshold for significant p-value
    :param threshold1: threshold for extremely significant p-value
    :param chr_names: list of chromosome names selected to plot
    :param log: whether to plot -log10(value)
    :param colors: list of colors for each chromosome
    :param style: {'line', 'scatter'}, default 'line', style of plot
    :param fig: matplotlib figure object
    :return: matplotlib figure object
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
    if isinstance(value_col, int):
        value_col = [value_col]
    ydata = df.iloc[:, value_col].to_numpy()

    if log:
        ydata = -np.log10(ydata)

    # transform annotation
    if ann_col is not None:
        # drop empty rows if values are NA in all annotation columns
        ann_data = df.iloc[:, ann_col].dropna(how='all')
    else:
        ann_data = None

    # set colors for each chromosome
    if colors is None:
        colors = cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    else:
        colors = cycle(colors)
    cdata = xdata.groupby(level=0, group_keys=False).apply(lambda x: x.assign(colors=next(colors))).iloc[:, -1]
    # calculate chromosome length
    chr_len = xdata.groupby(level=0).max()
    # add gap between chromosomes
    chr_len_new = chr_len.transform(lambda x: x + x.sum() / 100)
    # calculate chromosome center location
    chr_start = chr_len_new.cumsum() - chr_len_new
    chr_center = chr_start + chr_len_new / 2

    # transform x loci based on cum sum of chromosome length
    xdata_new = pd.concat([group + chr_start.loc[name] for name, group in xdata.groupby(level=0)])

    # get figure object
    if fig is None:
        fig = plt.gcf()

    # create GridSpec object
    gs = fig.add_gridspec(len(value_col), 1)
    # add plot
    for i in range(len(value_col)):
        ax = fig.add_subplot(gs[i, 0])
        ax.set_title(df.columns[value_col[i]])

        if style == 'scatter':
            ax.scatter(xdata_new.to_numpy().flatten(), ydata[:, i], color=cdata, s=0.5)
        elif style == 'line':
            ax.vlines(xdata_new.to_numpy().flatten(), 0, ydata[:, i], color=cdata, linewidth=0.5)
        else:
            raise ValueError('style must be "line" or "scatter"')

        if i != len(value_col) - 1:
            ax.set_xticks([])
        else:
            ax.set_xticks(chr_center.iloc[:, 0].to_numpy(), chr_center.index.to_list())

        if threshold0 is not None:
            ax.axhline(-np.log10(threshold0), color='gray', linestyle='--', linewidth=0.5)
        if threshold1 is not None:
            ax.axhline(-np.log10(threshold1), color='gray', linewidth=0.5)

        # add annotation
        if ann_data is not None:
            ann_index = ann_data.index.to_numpy()
            ann_x = xdata_new.iloc[ann_index, :].to_numpy().flatten()
            ann_y = ydata[ann_index, i]
            # add an annotation for each point
            for xi, yi, text in zip(ann_x, ann_y, ann_data.to_numpy()):
                ax.annotate(
                    text,
                    xy=(xi, yi), xycoords='data',
                    xytext=(30, 10), textcoords='offset points',
                    ha='left', va='top',
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->", fc="C0", ec="C0", lw=0.5)
                )
    return fig


def QQPlot(
        df: pd.DataFrame,
        value_col: int = 2,
        ax: axes.Axes = None):
    """
    QQ plot

    Parameters
    ----------
    :param df: dataframe with p-value column at least
    :param value_col: column index of value
    :param ax: matplotlib axes object
    :return: matplotlib axes object
    """
    if ax is None:
        ax = plt.gca()

    p = df.iloc[:, value_col]
    # return two tuple: (osm, osr), (slope, intercept, r)
    # 1st tuple contains two arrays:
    # first is theoretical quantiles, second is sample quantiles
    res = stats.probplot(p, dist=stats.uniform)
    osm, osr = res[0]
    ax.scatter(-np.log10(osm), -np.log10(osr), c='k', s=5)
    # slope, intercept, r = res[1]
    reg = stats.linregress(-np.log10(osm), -np.log10(osr))
    ax.plot(-np.log10(osm), -np.log10(osm) * reg[0] + reg[1], 'b--')
    # calculate confidence interval for linear regression
    x = -np.log10(osm)
    y = -np.log10(osm) * reg[0] + reg[1]
    std = x.std()
    mean = x.mean()
    n = len(x)
    y_err = std * np.sqrt(1 / n + (x - mean) ** 2 / np.sum((x - mean) ** 2))
    ax.fill_between(x, y - y_err, y + y_err, alpha=0.2)
    return ax


def LDHeatmapPlot(
        df: pd.DataFrame,
        plot_diag: bool = True,
        plot_value: bool = False,
        cmap: any([str, 'Colormap']) = None,
        ax: axes.Axes = None):
    """
    LD heatmap plot

    Parameters
    ----------
    :param df: a dataframe with MxN R^2 values only
    :param plot_diag: whether to plot diagonal
    :param plot_value: whether to plot value
    :param cmap: colormap, can be a string or a matplotlib colormap object
    :param ax: matplotlib axes object
    :return: matplotlib axes object
    """
    if ax is None:
        ax = plt.gca()

    # get R^2 values
    data = df.to_numpy()

    # whether to plot diagonal
    start = 0
    stop = data.shape[0]
    n = data.shape[0]
    if not plot_diag:
        start = 1

    # set colormap
    if cmap is None:
        cmap = 'Reds'
    cmap = mpl.colormaps.get_cmap(cmap)
    norm = mpl.colors.Normalize(vmin=0, vmax=1)

    # get patch collection to plot
    patches = []
    values = []
    # for loop along the lower triangle of the matrix, then get the diagonal values
    for i in np.arange(start, stop):
        diag_values = np.diag(data, -i)
        values.extend(diag_values)
        for j in np.arange(0.5, len(diag_values) + 0.5):
            patches.append(RegularPolygon((j + i * 0.5, (n - i) / 2), numVertices=4, radius=0.5))

    patch_collection = PatchCollection(patches)
    patch_collection.set_array(values)
    patch_collection.set_cmap(cmap)
    patch_collection.set_norm(norm)

    ax.add_collection(patch_collection)
    ax.set_aspect('equal')
    ax.set_xlim(-0.1 + start - 0.5, n - start + 0.5 + 0.1)
    ax.set_ylim(-0.1, (n - start) / 2 + 0.5 + 0.1)

    # Add color bar
    cax = ax.inset_axes([0.8, 0.01, 0.03, 0.5])
    ax.figure.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, shrink=.5, label=r"$R^2$")

    # Add text
    # Loop over data dimensions and create text annotations.
    # Change the text's color depending on the data.
    # Change the font's size depending on the patch.
    if plot_value:
        text_colors = ("black", "white")
        color_array = patch_collection.get_array()
        threshold = patch_collection.norm(color_array.max()) / 2
        for i, p in enumerate(patches):
            text = ax.text(p.xy[0], p.xy[1], "{:.2f}".format(values[i]),
                           ha="center", va="center",
                           color=text_colors[int(values[i] > threshold)])
            patch_bbox = p.get_window_extent()
            text_width = text.get_window_extent().transformed(ax.transData.inverted()).width
            font_size = text.get_fontsize()
            while font_size > 1 and text_width > patch_bbox.width / 2:
                font_size -= 1
                text.set_fontsize(font_size)
                text_width = text.get_window_extent().transformed(ax.transData.inverted()).width

    # Turn axis off
    ax.set_axis_off()

    return ax


def GeoMapPlot(df: pd.DataFrame,
               lon_col: int = 0,
               lat_col: int = 1,
               value_col: int | list = 2,
               map_file: str = None,
               bound_file: str = None,
               ax: axes.Axes = None):
    """
    Geographical map plot

    Parameters
    ----------
    :param df: a dataframe with at least three columns: 'latitude', 'longitude' and 'value'
    :param lon_col: column index of longitude
    :param lat_col: column index of latitude
    :param value_col: column index of value
    :param map_file: map file path, should be a shapefile or a geojson file or json file
    :param bound_file: boundary file path, should be a shapefile or a geojson file or json file
    :param ax: matplotlib axes object
    :return: matplotlib axes object
    """
    xy_data = df.iloc[:, [lon_col, lat_col]].to_numpy()
    if isinstance(value_col, int):
        value_col = [value_col]
    value_data = df.iloc[:, value_col].to_numpy()
    value_sum = np.sum(value_data, axis=1)
    value_percent = np.cumsum(value_data, axis=1) / value_sum[:, None]

    if ax is None:
        ax = plt.gca()

    # plot map
    if map_file is None:
        raise ValueError("map_file is not specified")
    map_data = gpd.read_file(map_file)
    map_data.plot(ax=ax, color='white', edgecolor='lightgray')

    # plot boundary
    if bound_file is not None:
        bound_data = gpd.read_file(bound_file)
        bound_data.plot(ax=ax, color='mediumpurple', linewidth=3, alpha=0.4)

    ax.set_aspect('equal')

    # scatter plot for the size of each point
    scatter = ax.scatter(xy_data[:, 0], xy_data[:, 1], s=value_sum, c='white', edgecolors='m')

    # scatter plot using pie shape marker to show the percentage of each value
    # https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Wedge.html#matplotlib.patches.Wedge
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for i in range(value_percent.shape[0]):
        for j in range(value_percent.shape[1]):
            theta = 2 * np.pi * np.linspace(value_percent[i, j - 1] if j > 0 else 0, value_percent[i, j], 100)
            vertices = np.column_stack([np.cos(theta), np.sin(theta)])
            ax.scatter(xy_data[i, 0], xy_data[i, 1],
                       marker=np.append(vertices, [[0, 0]], axis=0),
                       s=value_sum[i], c=colors[j], linewidths=0)

    # add the legend for different values
    legend = ax.legend(
        handles=[Patch(color=colors[i], label=l) for i, l in enumerate(df.iloc[:, value_col].columns)],
        loc="lower left", title="Values", frameon=False
    )
    ax.add_artist(legend)

    # produce a legend with a cross-section of sizes from the scatter
    handles, labels = scatter.legend_elements(prop="sizes", num=5, markeredgecolor='m', markerfacecolor='w')
    ax.legend(handles, labels, loc="lower right", title="Sizes", frameon=False)

    return ax


def HistPlot(df: pd.DataFrame,
             value_col: int | list = 0,
             bins: int | list = 10,
             orientation: str = 'vertical',
             plot_cdf: bool = False,
             plot_ci: bool = False,
             plot_pdf: bool = False,
             ax: axes.Axes = None):
    """
    Histogram plot

    Parameters
    ----------
    :param df: a dataframe with at least one column: 'value'
    :param value_col: column index of value
    :param bins: number of bins, can be int or list
    :param orientation: 'vertical' or 'horizontal'
    :param plot_cdf: whether to plot the cumulative distribution function
    :param plot_ci: whether to plot the confidence interval
    :param plot_pdf: whether to plot the probability density function
    :param ax: matplotlib axes object
    :return: matplotlib axes object
    """
    if isinstance(value_col, int):
        value_col = [value_col]
    data = df.iloc[:, value_col].to_numpy()
    # calculate the mean and standard deviation of data
    mean = np.mean(data, axis=0)
    std = np.std(data, axis=0)

    if ax is None:
        ax = plt.gca()

    # plot histogram
    n, bins, patches = ax.hist(data, bins=bins, orientation=orientation,
                               edgecolor='white', linewidth=.1, density=True)
    if orientation == 'vertical':
        ax.spines[['top', 'bottom']].set_visible(False)
        ax.set_ylabel("Density")
        ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    else:
        ax.spines[['left', 'right']].set_visible(False)
        ax.set_xlabel("Density")
        ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
    # add legend
    # plot the confidence interval for each value
    if plot_ci:
        for m, s in zip(mean, std):
            ci = stats.norm.interval(0.95, loc=m, scale=s)
            # add the confidence interval to the plot
            if orientation == 'vertical':
                ax.axvline(x=ci[0], c='C7', linewidth=1, linestyle='--')
                ax.axvline(x=ci[1], c='C7', linewidth=1, linestyle='--')
            else:
                ax.axhline(y=ci[0], c='C7', linewidth=1, linestyle='--')
                ax.axhline(y=ci[1], c='C7', linewidth=1, linestyle='--')
            # change the color of the patch according to the confidence interval
            for i, p in enumerate(patches):
                if bins[i] < ci[0]:
                    p.set_facecolor('C2')
                elif bins[i] > ci[1]:
                    p.set_facecolor('C3')
                else:
                    p.set_facecolor('C7')
    # plot the probability density function
    if plot_pdf:
        for m, s in zip(mean, std):
            ax.plot(bins, stats.norm.pdf(bins, loc=m, scale=s), c='black', linewidth=1)
    # plot the cumulative distribution function
    if plot_cdf:
        for m, s in zip(mean, std):
            if orientation == 'vertical':
                ax_twin = ax.twinx()
                ax_twin.plot(bins, stats.norm.cdf(bins, loc=m, scale=s), c='black', linewidth=1)
                ax_twin.spines[['top', 'bottom']].set_visible(False)
                ax_twin.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0, symbol=None))
                ax_twin.set_ylabel('Cumulative (%)')
                ax_twin.set_ylim(0, 1)
            else:
                ax_twin = ax.twiny()
                ax_twin.plot(stats.norm.cdf(bins, loc=m, scale=s), bins, c='black', linewidth=1)
                ax_twin.spines[['left', 'right']].set_visible(False)
                ax_twin.xaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0, symbol=None))
                ax_twin.set_xlabel('Cumulative (%)')
                ax_twin.set_xlim(0, 1)

    return ax


def PiWithFstPlot(df: pd.DataFrame,
                  pi_ratio_col: int = 0,
                  fst_col: int = 1,
                  fig: figure.Figure = None):
    """
    Pi ratio with Fst plot

    Parameters
    ----------
    :param df: a dataframe with at least two columns: 'fst' and 'pi ratio'
    :param pi_ratio_col: column index of pi ratio value
    :param fst_col: column index of fst value
    :param fig: matplotlib figure object
    :return: matplotlib figure object
    """

    # data preparation
    xdata = df.iloc[:, pi_ratio_col].to_numpy()
    # get confidence interval of pi ratio using normal distribution with 95% confidence level
    x_ci = stats.norm.interval(0.95, loc=np.mean(xdata), scale=np.std(xdata))

    # calculate the mean and standard deviation of fst
    ydata = df.iloc[:, fst_col].to_numpy()
    # get confidence interval of fst using normal distribution with 95% confidence level
    y_ci = stats.norm.interval(0.95, loc=np.mean(ydata), scale=np.std(ydata))

    # subset data by confidence interval
    data = np.column_stack((xdata, ydata))
    data_top0 = data[np.all([data[:, 0] < x_ci[0], data[:, 1] > y_ci[1]], axis=0)]
    data_top1 = data[np.all([data[:, 0] > x_ci[1], data[:, 1] > y_ci[1]], axis=0)]

    # get figure object
    if fig is None:
        fig = plt.gcf()

    # create GridSpec object
    gs = fig.add_gridspec(2, 2, width_ratios=(3, 1), height_ratios=(1, 3),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.02, hspace=0.02)

    # create axes object
    ax = fig.add_subplot(gs[1, 0])
    ax_x = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_y = fig.add_subplot(gs[1, 1], sharey=ax)

    # plot scatter
    ax.scatter(xdata, ydata, s=10, c='black')
    if data_top0.shape[0] > 0:
        ax.scatter(data_top0[:, 0], data_top0[:, 1], s=10, c='C2')
    if data_top1.shape[0] > 0:
        ax.scatter(data_top1[:, 0], data_top1[:, 1], s=10, c='C3')
    legend_lines = [Line2D([0], [0], color='black', lw=4),
                    Line2D([0], [0], color='C2', lw=4),
                    Line2D([0], [0], color='C3', lw=4)]
    ax.legend(
        legend_lines,
        [
            'Whole genome',
            'Selected region (<{:.2f},>{:.2f})'.format(x_ci[0], y_ci[1]),
            'Selected region (>{:.2f},>{:.2f})'.format(x_ci[1], y_ci[1])
        ],
        loc='upper left', frameon=False
    )
    ax.set_xlabel('Log10(Pi ratio)')
    ax.set_ylabel('Fst')
    ax.axvline(x=x_ci[0], c='C7', linewidth=1, linestyle='--')
    ax.axvline(x=x_ci[1], c='C7', linewidth=1, linestyle='--')
    ax.axhline(y=y_ci[0], c='C7', linewidth=1, linestyle='--')
    ax.axhline(y=y_ci[1], c='C7', linewidth=1, linestyle='--')
    # change the color of the point combine x and y confidence interval

    # plot hist for x
    HistPlot(df, 0, bins=50, orientation='vertical', ax=ax_x, plot_pdf=False, plot_cdf=True, plot_ci=True)

    # plot hist for y
    HistPlot(df, 1, bins=50, orientation='horizontal', ax=ax_y, plot_pdf=False, plot_cdf=True, plot_ci=True)

    return fig
