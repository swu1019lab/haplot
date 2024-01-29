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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import figure, axes, ticker, transforms
from matplotlib.patches import RegularPolygon, Patch, Rectangle, FancyArrowPatch, ConnectionPatch
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.transforms import Affine2D
import matplotlib.pyplot as plt

import geopandas as gpd
import networkx as nx
from haplot.utils import AnchoredSizeLegend


def BoxPlot(df: pd.DataFrame, by: str = 'column', ax: axes.Axes = None):
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


def ScatterPlot(df: pd.DataFrame,
                x_col: int = 0,
                y_col: int = 1,
                color_col: int = None,
                size_col: int = None,
                color: str = None,
                size: float = None,
                title: str = None,
                x_label: str = None,
                x_unit: float = 1e-6,
                y_label: str = None,
                log_y: bool = False,
                plot_threshold: bool = False,
                threshold0: float = None,
                threshold1: float = None,
                plot_cmap: bool = False,
                cmap: any([str, 'Colormap']) = None,
                ax: axes.Axes = None):
    """
    Scatter plot

    Parameters
    ----------
    :param df: a dataframe with at least two columns: 'x' and 'y'
    :param x_col: column index of x
    :param y_col: column index of y
    :param color_col: column index of color
    :param size_col: column index of size
    :param color: color of scatter
    :param size: size of scatter
    :param title: title of plot
    :param x_label: label of x-axis
    :param x_unit: unit of x-axis, default 1e-6, unit is Mb
    :param y_label: label of y-axis
    :param log_y: whether to plot log10(y)
    :param plot_threshold: whether to plot threshold
    :param threshold0:  for significant p-value
    :param threshold1: threshold for extremely significant p-value
    :param plot_cmap: whether to plot colormap for color_col
    :param cmap: can be a string or a matplotlib colormap object
    :param ax: matplotlib axes object
    :return: matplotlib axes object
    """
    if ax is None:
        ax = plt.gca()

    # get x and y data
    x = df.iloc[:, x_col].to_numpy()
    y = df.iloc[:, y_col].to_numpy()

    # plot log10(y)
    if log_y:
        y = -np.log10(y)

    # get color data
    if color_col is not None:
        c = df.iloc[:, color_col].to_numpy()
    else:
        c = color

    # get size data
    if size_col is not None:
        s = df.iloc[:, size_col].to_numpy()
    else:
        s = size

    # plot scatter
    if plot_cmap:
        if cmap is None:
            cmap = 'Reds'
        cmap = mpl.colormaps.get_cmap(cmap)
        norm = mpl.colors.Normalize()
        ax.scatter(x, y, c=c, s=s, cmap=cmap, norm=norm)
        cax = ax.inset_axes([1.04, 0.2, 0.03, 0.6])
        ax.figure.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, shrink=.5, label=r"$R^2$")
    else:
        ax.scatter(x, y, c=c, s=s)

    # plot threshold
    if plot_threshold:
        if threshold0 is not None:
            ax.axhline(threshold0, color='C0', linestyle='--', linewidth=.5)
        if threshold1 is not None:
            ax.axhline(threshold1, color='C1', linewidth=.5)

    # set x ticks format
    formatter = ticker.FuncFormatter(lambda i, pos: '{:.1f}'.format(i * x_unit))
    ax.xaxis.set_major_formatter(formatter)

    # set x and y label
    if x_label is not None:
        ax.set_xlabel(x_label)

    if y_label is not None:
        ax.set_ylabel(y_label)

    # set title
    if title is not None:
        ax.set_title(title)

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
        chr_used: list = None,
        chr_units: float = 1e-6,
        value_layout: str = 'stack',
        log: bool = False,
        colors: list = None,
        style: str = 'line',
        ax: axes.Axes = None):
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
    :param chr_names: list of chromosome names to replace the original names
    :param chr_used: list of chromosome names used to plot
    :param chr_units: unit of chromosome length, default 1e-6, unit is Mb
    :param value_layout: {'stack', 'overlay'}, default 'stack', layout of value
    :param log: whether to plot -log10(value)
    :param colors: list of colors for each chromosome
    :param style: {'line', 'scatter'}, default 'line', style of plot
    :param ax: a matplotlib.axes.Axes object
    :return: list of Axes objects
    """
    # sort data by chromosome and position
    new_df = df.sort_values(by=[df.columns[chr_col], df.columns[pos_col]])
    new_df.dropna(inplace=True, how='any', subset=[df.columns[chr_col], df.columns[pos_col]])
    new_df.reset_index(drop=True, inplace=True)

    # filter data
    if chr_used is not None:
        new_df = new_df.loc[new_df.iloc[:, chr_col].isin(chr_used), :]

    # transform xdata
    xdata = new_df.iloc[:, [chr_col, pos_col]]
    xdata.columns = ['chr', 'ps']
    xdata.set_index('chr', inplace=True)
    xdata = xdata.transform(lambda x: x * chr_units)

    # transform ydata
    if isinstance(value_col, int):
        value_col = [value_col]
    ydata = new_df.iloc[:, value_col].to_numpy()

    if log:
        ydata = -np.log10(ydata)

    # transform annotation
    if ann_col is not None:
        # drop empty rows if values are NA in all annotation columns
        ann_data = new_df.iloc[:, ann_col].dropna(how='all')
    else:
        ann_data = None

    # set colors for each chromosome
    if colors is None:
        # colors = cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])
        colors = cycle(['skyblue', 'steelblue'])
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
    chr_labels = chr_center.index.to_numpy()
    if chr_names is not None and len(chr_names) == chr_center.shape[0]:
        chr_labels = chr_names

    # transform x loci based on cum sum of chromosome length
    xdata_new = pd.concat([group + chr_start.loc[name] for name, group in xdata.groupby(level=0)])

    # get axes object
    if ax is None:
        ax = plt.gca()

    # create GridSpec object
    if value_layout == 'stack':
        ax_divider = make_axes_locatable(ax)
        axs = [ax] + [ax_divider.append_axes("bottom", size="100%", pad=0.8, sharex=ax) for _ in
                      range(len(value_col) - 1)]
    elif value_layout == 'overlay':
        axs = [ax] * len(value_col)
    else:
        raise ValueError('value_layout must be "stack" or "overlay"')
    # add plot
    for i in range(len(value_col)):
        axs[i].set_title(new_df.columns[value_col[i]])

        if style == 'scatter' and value_layout == 'stack':
            axs[i].scatter(xdata_new.to_numpy().flatten(), ydata[:, i], color=cdata, s=0.5)
        elif style == 'scatter' and value_layout == 'overlay':
            axs[i].scatter(xdata_new.to_numpy().flatten(), ydata[:, i], s=0.5)
        elif style == 'line' and value_layout == 'stack':
            axs[i].vlines(xdata_new.to_numpy().flatten(), 0, ydata[:, i], color=cdata, linewidth=0.5)
        elif style == 'line' and value_layout == 'overlay':
            axs[i].vlines(xdata_new.to_numpy().flatten(), 0, ydata[:, i], color=next(colors), linewidth=0.5)
        else:
            raise ValueError('style must be "line" or "scatter"')

        axs[i].set_xticks(chr_center.iloc[:, 0].to_numpy(), chr_labels)

        if threshold0 is not None:
            axs[i].axhline(-np.log10(threshold0), color='gray', linestyle='--', linewidth=0.5)
        if threshold1 is not None:
            axs[i].axhline(-np.log10(threshold1), color='gray', linewidth=0.5)

        # add annotation
        if ann_data is not None:
            ann_index = ann_data.index.to_numpy()
            ann_x = xdata_new.iloc[ann_index, :].to_numpy().flatten()
            ann_y = ydata[ann_index, i]
            # add an annotation for each point
            for xi, yi, text in zip(ann_x, ann_y, ann_data.to_numpy()):
                axs[i].annotate(
                    text,
                    xy=(xi, yi), xycoords='data',
                    xytext=(30, 10), textcoords='offset points',
                    ha='left', va='top',
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->", fc="C0", ec="C0", lw=0.5)
                )
    return axs


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

    p = df.iloc[:, value_col].dropna().to_numpy()
    # return two tuple: (osm, osr), (slope, intercept, r)
    # 1st tuple contains two arrays:
    # first is theoretical quantiles, second is sample quantiles
    res = stats.probplot(p, dist=stats.uniform)
    osm, osr = res[0]
    ax.scatter(-np.log10(osm), -np.log10(osr), c='skyblue', s=5)
    ax.plot(-np.log10(osm), -np.log10(osm), 'k--')
    return ax


def LDHeatmapPlot(
        df: pd.DataFrame,
        plot_diag: bool = True,
        plot_value: bool = False,
        plot_snp: bool = False,
        cmap: any([str, 'Colormap']) = None,
        ax: axes.Axes = None):
    """
    LD heatmap plot

    Parameters
    ----------
    :param df: a dataframe with MxN R^2 values only
    :param plot_diag: whether to plot diagonal
    :param plot_value: whether to plot value
    :param plot_snp: whether to plot SNP location
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
    ax.set_xlim(start * 0.5, stop - start * 0.5)
    ax.set_ylim(-0.1, (n - start) / 2 + 0.5 + 0.1)

    # Turn axis off
    ax.set_axis_off()

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

    # Add SNP location
    # get SNP location from dataframe index name 'pos'
    if plot_snp and 'pos' in df.index.names:
        snp_loc = df.index.get_level_values('pos').to_numpy()
        sx = (stop - start) / (np.max(snp_loc) - np.min(snp_loc))
        scale_loc = Affine2D(). \
            translate(-np.min(snp_loc), ty=0). \
            scale(sx=sx, sy=1). \
            translate(tx=start * 0.5, ty=0). \
            transform(np.column_stack([snp_loc, [1] * n]))
        line_collection = LineCollection([[[a[0], 1], [i + 0.5, 0]] for i, a in enumerate(scale_loc)], linewidths=.5)
        line_collection.set_color('black')

        ax_divider = make_axes_locatable(ax)
        ax_snp = ax_divider.append_axes("top", size="10%", pad="0%", sharex=ax)
        ax_snp.add_collection(line_collection)
        ax_snp.set_xlim(start * 0.5, stop - start * 0.5)
        ax_snp.set_ylim(0, 1)
        ax_snp.set_xticks([])
        ax_snp.set_yticks([])
        ax_snp.set_yticklabels([])
        ax_snp.spines.bottom.set_visible(False)

        return ax_snp

    return ax


def GeoMapPlot(df: pd.DataFrame,
               lon_col: int = 0,
               lat_col: int = 1,
               value_col: int | list = 2,
               colors: list = None,
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
    :param colors: a list of colors
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
    if map_file is not None:
        map_data = gpd.read_file(map_file)
        map_data.plot(ax=ax, color='white', edgecolor='lightgray')

    # plot boundary
    if bound_file is not None:
        bound_data = gpd.read_file(bound_file)
        bound_data.plot(ax=ax, color='mediumpurple', linewidth=3, alpha=0.4)

    ax.set_aspect('equal')

    # scatter plot for the size of each point
    scatter = ax.scatter(xy_data[:, 0], xy_data[:, 1], s=value_sum, edgecolors='m', facecolors=None)
    scatter.set_visible(False)

    # add annotation for each point

    # scatter plot using pie shape marker to show the percentage of each value
    # https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Wedge.html#matplotlib.patches.Wedge
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color'] if colors is None else colors
    for i in range(value_percent.shape[0]):
        for j in range(value_percent.shape[1]):
            theta = 2 * np.pi * np.linspace(value_percent[i, j - 1] if j > 0 else 0, value_percent[i, j], 100)
            vertices = np.column_stack([np.cos(theta), np.sin(theta)])
            ax.scatter(xy_data[i, 0], xy_data[i, 1],
                       marker=np.append(vertices, np.asarray([[0, 0]]), axis=0),
                       s=value_sum[i], c=colors[j], linewidths=0, zorder=3)

    # add the legend for different values
    legend = ax.legend(
        handles=[Patch(color=colors[i], label=l) for i, l in enumerate(df.iloc[:, value_col].columns)],
        loc='lower left', bbox_to_anchor=(0, 1.01, 1, 0.1),
        borderaxespad=0, ncols=3, mode='expand', frameon=False
    )
    ax.add_artist(legend)

    # produce a legend with a cross-section of sizes from the scatter
    # add legend for different node size
    loc = mpl.ticker.MaxNLocator(nbins=3)
    size_label = loc.tick_values(min(value_sum), max(value_sum))
    asl = AnchoredSizeLegend(
        size_label[1:],
        size_label[1:],
        label_size=5,
        loc='center',
        bbox_to_anchor=(0, 0, 0.2, 0.5),
        bbox_transform=ax.transAxes,
        pad=0.1, borderpad=0.5,
        frameon=False
    )
    ax.add_artist(asl)
    # handles, labels = scatter.legend_elements(prop="sizes", num=5, markeredgecolor='m', markerfacecolor='w')
    # ax.legend(handles, labels,
    #           loc="center", bbox_to_anchor=(0, 0, 0.1, 0.5), borderaxespad=0,
    #           title="Sizes", frameon=False)

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
            'Selected region (pi<{:.2f} and Fst>{:.2f})'.format(x_ci[0], y_ci[1]),
            'Selected region (pi>{:.2f} and Fst>{:.2f})'.format(x_ci[1], y_ci[1])
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


def GeneStrucPlot(
        df: pd.DataFrame,
        gene_col: int = 0,
        start_col: int = 1,
        end_col: int = 2,
        strand_col: int = 3,
        feature_col: int = 4,
        feature_config: dict = None,
        center: int = 0.5,
        jitter: bool = False,
        ax: axes.Axes = None):
    """
    Gene structure plot

    For each gene, the plot will show the gene structure, including CDS, intron and exon.\n
    Features
    --------
    The plot will show the gene structure, including CDS, intron and exon.\n

    Parameters
    ----------
    :param df: a dataframe with at least five columns: 'gene', 'start', 'end', 'strand', 'feature'
    :param gene_col: column index of gene name
    :param start_col: column index of start position
    :param end_col: column index of end position
    :param strand_col: column index of strand
    :param feature_col: column index of feature
    :param feature_config: a dictionary of feature configuration
    :param center: the center of the feature, should be a float number between 0 and 1
    :param jitter: whether to jitter the feature to avoid overlap
    :param ax: matplotlib axes object
    :return: matplotlib axes object
    """
    if feature_config is None:
        feature_config = {
            'CDS': {
                'color': 'C0',
                # height should be limited within [0, 1]
                'height': 0.15,
                'patch': 'rectangle',
                'zorder': 2,
                'alpha': 1
            },
            'intron': {
                'color': 'C5',
                # height should be limited within [0, 1]
                'height': 0.05,
                'patch': 'arrow',
                'zorder': 1,
                'alpha': 1
            },
            'exon': {
                'color': 'C0',
                # height should be limited within [0, 1]
                'height': 0.13,
                'patch': 'rectangle',
                'zorder': 2,
                'alpha': 1
            },
            '5UTR': {
                'color': 'C7',
                # height should be limited within [0, 1]
                'height': 0.1,
                'patch': 'rectangle',
                'zorder': 1,
                'alpha': 1
            },
            '3UTR': {
                'color': 'C7',
                # height should be limited within [0, 1]
                'height': 0.1,
                'patch': 'rectangle',
                'zorder': 1,
                'alpha': 1
            }
        }

    # get axes object
    if ax is None:
        ax = plt.gca()

    # scale data
    x_min = df.iloc[:, start_col].min() - 2000
    x_max = df.iloc[:, end_col].max() + 2000
    new_df = df.copy()
    new_df.iloc[:, [start_col, end_col]] = (df.iloc[:, [start_col, end_col]] - x_min) / (x_max - x_min)

    # plot gene structure
    for name, group in new_df.groupby(new_df.columns[gene_col]):
        # get gene strand
        gene_strand = group.iloc[0, strand_col]
        # get gene start and end
        gene_start = group.iloc[:, start_col].min()
        gene_end = group.iloc[:, end_col].max()
        # jitter feature
        if jitter:
            offset = np.random.uniform(-0.25, 0.25)
        else:
            offset = 0
        # get gene feature
        for feature in group.itertuples(name="Feature", index=False):
            # get feature location
            feature_start = feature[start_col]
            feature_end = feature[end_col]
            # get feature type
            feature_type = feature[feature_col]
            # get feature color
            feature_color = feature_config[feature_type]['color']
            # get feature height
            feature_height = feature_config[feature_type]['height']
            # get feature patch type
            feature_patch = feature_config[feature_type]['patch']
            # get feature z-order
            feature_zorder = feature_config[feature_type]['zorder']
            # get feature alpha
            feature_alpha = feature_config[feature_type]['alpha']
            # plot feature
            if feature_patch == 'rectangle':
                ax.add_patch(Rectangle((feature_start, center - feature_height / 2 + offset),
                                       feature_end - feature_start, feature_height,
                                       fc=feature_color, ec=feature_color,
                                       zorder=feature_zorder, alpha=feature_alpha))
            elif feature_patch == 'arrow' and gene_strand == '+':
                ax.add_patch(FancyArrowPatch((feature_start, center + offset),
                                             (feature_end, center + offset),
                                             arrowstyle="->", fc=feature_color,
                                             ec=feature_color, shrinkA=0, shrinkB=0,
                                             mutation_scale=10, zorder=feature_zorder,
                                             alpha=feature_alpha))
            elif feature_patch == 'arrow' and gene_strand == '-':
                ax.add_patch(FancyArrowPatch((feature_start, center + offset),
                                             (feature_end, center + offset),
                                             arrowstyle="<-", fc=feature_color,
                                             ec=feature_color, shrinkA=0, shrinkB=0,
                                             mutation_scale=10, zorder=feature_zorder,
                                             alpha=feature_alpha))
            else:
                raise ValueError('Unknown feature patch type: {}'.format(feature_patch))
        # add gene name
        ax.annotate(name,
                    xy=((gene_start + gene_end) / 2, center + offset),
                    xycoords='data',
                    xytext=(0, 30),
                    textcoords='offset points',
                    va='center', ha='center')
    # add legend for feature
    legend_elements = []
    for feature_type, feature_config in feature_config.items():
        legend_elements.append(Patch(facecolor=feature_config['color'],
                                     edgecolor=feature_config['color'],
                                     label=feature_type))
    ax.legend(handles=legend_elements,
              loc='lower left', bbox_to_anchor=(0, 1.3, 1, 0.1),
              borderaxespad=0, ncols=5, mode='expand', frameon=False)
    # set x-axis
    ax_twin = ax.twiny()
    ax_twin.set_xlim(x_min, x_max)
    # set x-axis ticks with kb unit
    ax_twin.set_xlabel('Position (kb)')
    ax_twin.xaxis.set_major_formatter(
        ticker.FuncFormatter(lambda x, pos: '{:.2f}'.format(x / 1000))
    )
    # hide ticks and spines of ax and ax_twin
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[["left", "bottom", "right"]].set_visible(False)
    ax_twin.set_yticks([])
    ax_twin.spines[["left", "bottom", "right"]].set_visible(False)
    # set grid
    ax_twin.xaxis.grid(
        True,
        linestyle='--',
        which='major',
        color='lightgrey',
        alpha=.5
    )
    # Hide these grid behind plot objects
    ax_twin.set_axisbelow(True)

    return ax_twin


def HapNetworkPlot(
        edge_data: pd.DataFrame,
        node_data: pd.DataFrame = None,
        layout: str = 'spring',
        colors: list = None,
        node_font_size: int = 12,
        weight_show: bool = True,
        weight_show_style: int = 0,
        ax: axes.Axes = None):
    """
    Plot haplotype network using minimum spanning tree.

    Parameters
    ----------
    :param edge_data: a dataFrame with three columns: source, target, weight
    :param node_data: a dataFrame with node data, at least one column named 'node' and other columns for data statistics
    :param layout: layout algorithm for networkx
    :param colors: a list of colors for scatter-pie plot
    :param node_font_size: font size for node label
    :param weight_show: whether to show weight of edges
    :param weight_show_style: 0 for number, 1 for symbol with "|"
    :param ax: axes object to plot on (default: None)
    :return: axes object
    """
    # get axes object
    if ax is None:
        ax = plt.gca()

    # create networkx graph
    G = nx.from_pandas_edgelist(edge_data,
                                edge_attr=['weight', 'color'],
                                create_using=nx.Graph())
    # Find the minimum spanning tree
    T = nx.minimum_spanning_tree(G)

    # get node layout and plot nodes with label
    if layout == 'spring':
        pos = nx.spring_layout(T)
    elif layout == 'spectral':
        pos = nx.spectral_layout(T)
    elif layout == 'random':
        pos = nx.random_layout(T)
    elif layout == 'shell':
        pos = nx.shell_layout(T)
    elif layout == 'circular':
        pos = nx.circular_layout(T)
    else:
        raise ValueError('Unknown layout: {}'.format(layout))
    nx.draw_networkx_labels(T, pos, font_size=node_font_size)

    # show node data with scatter-pie plot
    node_size = []
    prop_cycle = plt.rcParams['axes.prop_cycle']
    if colors is None:
        colors = prop_cycle.by_key()['color']
    if node_data is not None:
        # loop node name in graph, get position and stat data of each node
        for node in T.nodes:
            x, y = pos[node]
            stat_data = node_data[node_data['node'] == node].values[0][1:]
            node_size.append(np.sum(stat_data))
            # stat percentage of each node and make the sum of percentage to 1
            stat_percentage = np.insert(np.cumsum(stat_data) / np.sum(stat_data), 0, 0)
            # plot scatter-pie
            for i in range(len(stat_percentage) - 1):
                theta = 2 * np.pi * np.linspace(stat_percentage[i], stat_percentage[i + 1], num=100)
                vertices = np.column_stack((np.cos(theta), np.sin(theta)))
                ax.scatter(x, y,
                           np.sum(stat_data),
                           colors[i],
                           marker=np.append(vertices, np.asarray([[0, 0]]), axis=0),
                           linewidths=0,
                           zorder=2)

        # add legend for different node size
        loc = mpl.ticker.MaxNLocator(nbins=2)
        node_size_label = loc.tick_values(min(node_size), max(node_size))
        asl = AnchoredSizeLegend(
            node_size_label[1:],
            node_size_label[1:],
            label_size=6,
            loc='lower left',
            bbox_to_anchor=(1.02, 0., 0.1, 1),
            bbox_transform=ax.transAxes,
            pad=0.1, borderpad=0.5,
            frameon=False
        )
        ax.add_artist(asl)

        # add legend for different types within same node
        legend_elements = []
        for i, column in enumerate(node_data.columns[1:]):
            legend_elements.append(
                Patch(facecolor=colors[i], edgecolor=colors[i], label=column)
            )
        ax.legend(handles=legend_elements,
                  loc='lower left', bbox_to_anchor=(0, 1.02, 1, 0.1),
                  borderaxespad=0, ncols=5, mode='expand', frameon=False)

    # get edge weights and plot edges with label
    if weight_show:
        if weight_show_style == 0:
            nx.draw_networkx_edge_labels(
                T, pos, edge_labels={(u, v): d["weight"] for u, v, d in T.edges(data=True)})
        elif weight_show_style == 1:
            nx.draw_networkx_edge_labels(
                T, pos, edge_labels={(u, v): d["weight"] * "|" for u, v, d in T.edges(data=True)}
            )
        else:
            raise ValueError('Unknown weight show style: {}'.format(weight_show_style))
    nx.draw_networkx_edges(T, pos, edge_color=[d["color"] for u, v, d in T.edges(data=True)])

    # hide axis
    ax.set_axis_off()

    return ax


def VennNetworkPlot(
        edge_data: pd.DataFrame,
        edge_width: int | float = 1.0,
        edge_style: int = 1,
        source_node_size: int | float = 100,
        source_font_size: int = 10,
        target_node_size: int | float = 5,
        target_font_size: int = 10,
        show_node_margin: bool = True,
        show_node_color: bool = False,
        show_source_label: bool = True,
        show_target_label: bool = False,
        r: float = 0.5,
        k: float = None,
        ax: axes.Axes = None):
    """
    Plot venn network.

    Parameters
    ----------
    :param edge_data: a dataFrame with three columns: source, target, color
    :param edge_width: width of edges
    :param edge_style: 1 for line, 2 for curve
    :param source_node_size: size of source nodes
    :param source_font_size: font size of source labels
    :param target_node_size: size of target nodes
    :param target_font_size: font size of target labels
    :param show_node_margin: whether to show margin of nodes when style is curve
    :param show_node_color: whether to show node color according to different source
    :param show_source_label: whether to show source labels
    :param show_target_label: whether to show target labels
    :param r: radius of the fixed nodes based on circle center (0, 0)
    :param k: optimal distance between nodes for spring layout
    :param ax: axes object to plot on (default: None)
    :return: axes object
    """
    # get axes object
    if ax is None:
        ax = plt.gca()

    # create networkx graph
    G = nx.from_pandas_edgelist(edge_data, edge_attr=['color'], create_using=nx.Graph())

    # get fixed nodes from source column
    fixed_nodes = edge_data['source'].unique()
    print(fixed_nodes)
    # set position of fixed nodes
    fixed_nodes_pos = {node: (r * np.cos(2 * np.pi * i / len(fixed_nodes)),
                              r * np.sin(2 * np.pi * i / len(fixed_nodes)))
                       for i, node in enumerate(fixed_nodes)}

    node_size = [target_node_size if u not in fixed_nodes else source_node_size for u, d in G.nodes(data=True)]
    edge_color = [d["color"] for u, v, d in G.edges(data=True)]

    node_options = {
        "node_color": "lightgray",
        "node_size": node_size,
    }

    edge_style1 = {
        'edge_color': edge_color,
        'width': edge_width,
    }
    edge_style2 = {
        'edge_color': edge_color,
        'width': edge_width,
        'arrows': True,
        'arrowsize': 10,
        'arrowstyle': None,
        'connectionstyle': 'arc3,rad=0.2'
    }
    edge_options = {
        1: edge_style1,
        2: edge_style2
    }
    # get node layout
    pos = nx.spring_layout(G, pos=fixed_nodes_pos, fixed=fixed_nodes, seed=42, k=k)

    # plot network
    if show_node_color:
        node_color = {}
        node_color.update(dict(zip(edge_data['source'], edge_data['color'])))
        node_color.update(dict(zip(edge_data['target'], edge_data['color'])))
        node_options['node_color'] = [node_color[node] for node in G.nodes]
    nx.draw_networkx_nodes(G, pos, **node_options)
    nx.draw_networkx_edges(G, pos, **edge_options[edge_style])

    if not show_node_margin:
        for p in ax.patches:
            p.shrinkA = p.shrinkB = 0.0

    # show fixed nodes or source nodes with label
    if show_source_label:
        nx.draw_networkx_labels(G, pos,
                                labels=dict(zip(fixed_nodes, fixed_nodes)),
                                font_size=source_font_size,
                                font_color='black',
                                font_weight='bold',
                                ax=ax)
    # show target nodes with label
    if show_target_label:
        nx.draw_networkx_labels(G, pos,
                                labels={node: node for node in G.nodes if node not in fixed_nodes},
                                font_size=target_font_size,
                                font_color='black',
                                font_weight='bold',
                                ax=ax)
    ax.set_axis_off()

    return ax


def HeatmapPlot(num_data: np.ndarray,
                txt_data: np.ndarray = None,
                cmap: str = 'PuBu',
                xlabel: str = None,
                ylabel: str = None,
                row_labels: list = None,
                col_labels: list = None,
                row_grids: list = None,
                col_grids: list = None,
                row_labels_as_color: bool = False,
                col_labels_as_color: bool = False,
                plot_row_ticks: bool = False,
                plot_col_ticks: bool = False,
                plot_cbar: bool = False,
                plot_grid: bool = False,
                ax: axes.Axes = None,
                **kwargs):
    """
    Plot heatmap.

    Parameters
    ----------
    :param num_data: an array with numerical data for heatmap
    :param txt_data: an array with text data for annotation, the size of texts_data should be the same as num_data
    :param cmap: colormap used in matplotlib
    :param xlabel: x label of axis
    :param ylabel: y label of axis
    :param row_labels: a list of row labels
    :param col_labels: a list of col labels
    :param row_grids: a list of row grids
    :param col_grids: a list of col grids
    :param row_labels_as_color: whether to use row labels as color
    :param col_labels_as_color: whether to use col labels as color
    :param plot_row_ticks: whether to plot row ticks
    :param plot_col_ticks: whether to plot col ticks
    :param plot_cbar: whether to plot color bar
    :param plot_grid: whether to plot grid
    :param ax: axes object to plot on (default: None)
    """
    # get axes object
    if ax is None:
        ax = plt.gca()

    # cmap = mpl.colormaps[cmap].resampled(3)
    ax_divider = make_axes_locatable(ax)

    # plot data
    im = ax.imshow(num_data, cmap=cmap, aspect='auto', vmin=0, vmax=2, **kwargs)
    # add texts
    threshold = im.norm(np.max(num_data)) / 2
    txt_colors = ['black', 'white']
    if txt_data is not None:
        txt_data = txt_data
        for i in range(num_data.shape[0]):
            for j in range(num_data.shape[1]):
                ax.text(j, i, txt_data[i, j],
                        ha='center', va='center',
                        color=txt_colors[int(im.norm(num_data[i, j]) > threshold)])
    # add ColorBar
    if plot_cbar:
        cax = ax.inset_axes([1.04, 0.2, 0.03, 0.6])
        bounds = [0, 1, 2, 3]
        norm = mpl.colors.BoundaryNorm(bounds, 3)
        cbar = ax.figure.colorbar(
            mpl.cm.ScalarMappable(cmap=mpl.colormaps[cmap].resampled(3), norm=norm),
            ax=ax,
            cax=cax,
            boundaries=bounds,
        )
        cbar.set_ticks([0.5, 1.5, 2.5], labels=['0', '1', '2'])
        # hide color bar ticks
        cbar.ax.tick_params(size=0)
    # add x, y label of axis
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # add column labels
    ax.set_xticks([])
    if col_labels is not None:
        if col_labels_as_color:
            ax_col = ax_divider.append_axes("top", size="5%", pad="1%", sharex=ax)
            # transform col labels to numerical data
            codes, uniques = pd.Series(col_labels).factorize()
            ax_col.imshow(codes.reshape(1, -1),
                          cmap=mpl.colormaps[cmap].resampled(len(uniques)),
                          aspect='auto')
            ax_col.set_axis_off()
        else:
            ax.set_xticks(np.arange(len(col_labels)), labels=col_labels)

        if not plot_col_ticks:
            ax.tick_params(axis='x', which='both', bottom=False)
    # add row labels
    ax.set_yticks([])
    if row_labels is not None:
        if row_labels_as_color:
            ax_row = ax_divider.append_axes("left", size="5%", pad="1%", sharey=ax)
            # transform row labels to numerical data
            codes, uniques = pd.Series(row_labels).factorize()
            ax_row.imshow(codes.reshape(-1, 1),
                          cmap=mpl.colormaps[cmap].resampled(len(uniques)),
                          aspect='auto')
            ax_row.set_axis_off()
        else:
            ax.set_yticks(np.arange(len(row_labels)), labels=row_labels)

        if not plot_row_ticks:
            ax.tick_params(axis='y', which='both', left=False)

    # hide spines
    ax.spines[:].set_visible(False)
    # add grid
    if plot_grid:
        if row_grids is not None:
            ax.set_yticks(np.asarray(row_grids) - .5, minor=True)
        else:
            ax.set_yticks(np.arange(num_data.shape[0] + 1) - .5, minor=True)

        if col_grids is not None:
            ax.set_xticks(np.asarray(col_grids) - .5, minor=True)
        else:
            ax.set_xticks(np.arange(num_data.shape[1] + 1) - .5, minor=True)
        ax.grid(which='minor',
                linestyle='-',
                color='w',
                linewidth=2)
        ax.tick_params(which='minor',
                       bottom=False,
                       left=False)

    return ax


def GeneWithHapHeatmap(df1: pd.DataFrame,
                       df2: pd.DataFrame,
                       cmap: str = 'PuBu',
                       fig: plt.Figure = None):
    """
    Plot gene structure and genotype heatmap.

    Parameters
    ----------
    :param df1: a dataFrame with gene structure data
    :param df2: a dataFrame with genotype data, columns are a multiIndex with two levels (chrom, pos)
    :param cmap: colormap used in matplotlib, e.g., Oranges, Reds, viridis, PuRd, PuBu, PuBuGn, ...
    :param fig: figure object to plot on (default: None)
    :return: a dictionary including all axes objects
    """
    # get figure object
    if fig is None:
        fig = plt.gcf()

    # create GridSpec with four axes
    gs = fig.add_gridspec(4, 1,
                          height_ratios=(1, 3, .5, .5),
                          hspace=0.1, wspace=0.02)
    axs = {
        'gene': fig.add_subplot(gs[0]),
        'heatmap': fig.add_subplot(gs[1]),
        'ref': fig.add_subplot(gs[2]),
        'alt': fig.add_subplot(gs[3])
    }
    # plot gene structure
    cax = GeneStrucPlot(df1, ax=axs['gene'])
    # plot genotype heatmap
    HeatmapPlot(
        df2.values,
        row_labels=df2.index.tolist(),
        cmap=cmap,
        plot_cbar=True,
        ax=axs['heatmap'],
        interpolation='none')
    # get ref and alt allele, e.g., A, T, G, C, AAAA, CCCC
    ref_allele = df2.columns.get_level_values(2)
    alt_allele = df2.columns.get_level_values(3)
    ref_allele_df = pd.Series(ref_allele).str.split('', expand=True).iloc[:, 1:-1].T
    HeatmapPlot(
        np.ones(ref_allele_df.shape),
        txt_data=ref_allele_df.values,
        cmap=cmap,
        ylabel='Ref allele',
        plot_grid=True,
        ax=axs['ref']
    )
    alt_allele_df = pd.Series(alt_allele).str.split('', expand=True).iloc[:, 1:-1].T
    HeatmapPlot(
        np.ones(alt_allele_df.shape),
        txt_data=alt_allele_df.values,
        cmap=cmap,
        ylabel="Alt allele",
        plot_grid=True,
        ax=axs['alt']
    )
    # connect two axes with connection patch
    # chromosome position must be int number
    pos = df2.columns.get_level_values(1).to_numpy(dtype=np.int64)
    points0 = [(i, .3) for i in pos]
    points1 = [(i, 0) for i in pos]
    points2 = [(i, 1) for i in range(len(pos))]

    transA = transforms.blended_transform_factory(cax.transData, cax.transAxes)
    transB = transforms.blended_transform_factory(axs['heatmap'].transData, axs['heatmap'].transAxes)
    if len(points1) == len(points2):
        for i in range(len(points1)):
            fig.add_artist(
                ConnectionPatch(
                    xyA=points0[i],
                    xyB=points1[i],
                    coordsA=transA,
                    coordsB=transA,
                    axesA=cax,
                    axesB=cax,
                    arrowstyle="-",
                    color="darkred",
                    linestyle="--",
                    linewidth=0.5,
                    zorder=0.1
                )
            )
            fig.add_artist(
                ConnectionPatch(
                    xyA=points1[i],
                    xyB=points2[i],
                    coordsA=transA,
                    coordsB=transB,
                    axesA=cax,
                    axesB=axs['heatmap'],
                    arrowstyle="-",
                    color="darkred",
                    linewidth=0.5,
                )
            )
    else:
        raise ValueError('Length of points1 and points2 are not equal.')

    return axs

