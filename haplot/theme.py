# -*- coding: utf-8 -*-
# @Time    : 2023/6/28 10:56
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : theme.py


from matplotlib import figure, axes
import matplotlib.pyplot as plt
from itertools import cycle


class Theme(object):
    def __str__(self):
        return "Theme"

    __repr__ = __str__

    def __call__(self, *args, **kwargs):
        return self.apply(*args, **kwargs)

    @staticmethod
    def apply(ax: axes.Axes = None):
        if ax is None:
            ax = plt.gca()

        # Set edge color for top and right spines
        ax.spines['top'].set_edgecolor('lightgrey')
        ax.spines['right'].set_edgecolor('lightgrey')

        # ax.spines[['top', 'right']].set_visible(False)

        ax.yaxis.grid(
            True,
            linestyle='--',
            which='major',
            color='lightgrey',
            alpha=.5
        )

        # Hide these grid behind plot objects
        ax.set_axisbelow(True)

        return ax


class BoxTheme(object):
    """
    Boxplot theme
    """
    def __str__(self):
        return "BoxTheme"

    __repr__ = __str__

    def __call__(self, *args, **kwargs):
        return self.apply(*args, **kwargs)

    @staticmethod
    def apply(ax: axes.Axes):
        if ax is None:
            ax = plt.gca()

        # color palette:
        if len(ax.patches) == 0:
            print("No patches in boxplot.")
            return ax

        box_face_colors = cycle(['skyblue', 'steelblue'])
        box_edge_colors = cycle(['black'])

        for p in ax.patches:
            p.set_facecolor(next(box_face_colors))
            p.set_edgecolor(next(box_edge_colors))

        # Add a horizontal grid to the boxplot, but make it very light in color
        ax.yaxis.grid(
            True,
            linestyle='--',
            which='major',
            color='lightgrey',
            alpha=.5
        )

        # Hide these grid behind plot objects
        ax.set_axisbelow(True)

        # Set edge color for top and right spines
        ax.spines['top'].set_edgecolor('lightgrey')
        ax.spines['right'].set_edgecolor('lightgrey')
        ax.spines['bottom'].set_position(('outward', 10))

        ax.spines['left'].set_bounds(ax.get_ybound())
        ax.spines['right'].set_bounds(ax.get_ybound())

        # Rotate the tick labels and set their alignment.
        # ax.tick_params(axis='x', colors='green', size=10, width=5, rotation=45)

        # Add custom rectangle annotation under axes
        for label in ax.xaxis.get_ticklabels():
            # label is a Text instance
            label.set_fontsize(16)

        tick_face_colors = cycle(['skyblue', 'steelblue'])
        tick_edge_colors = cycle(['skyblue', 'steelblue'])

        for tick in ax.xaxis.get_major_ticks():
            # line is a Line2D instance
            tick.tick1line.set_markerfacecolor(next(tick_face_colors))
            tick.tick1line.set_markeredgecolor(next(tick_edge_colors))

        ax.tick_params(axis='x',
                       pad=5,
                       size=5,  # tick height
                       width=ax.get_window_extent().width / len(ax.xaxis.get_major_ticks()) * .72,  # tick width
                       rotation=0)
        return ax


class ManhattanTheme(object):
    def __str__(self):
        return "ManhattanTheme"

    __repr__ = __str__

    def __call__(self, *args, **kwargs):
        return self.apply(*args, **kwargs)

    @staticmethod
    def apply(ax: axes.Axes = None,
              title: str = None,
              xlabel: str = None,
              ylabel: str = None,
              xspan: tuple = None,
              yspan: tuple = None,
              xticks: list = None,
              xticks_label: list = None,
              xticks_label_rotation: float = 0,
              ):
        if ax is None:
            ax = plt.gca()

        Theme.apply(ax)
        if title is not None:
            ax.set_title(title)
        if xticks is not None and xticks_label is not None:
            ax.set_xticks(xticks, xticks_label)

        if xlabel is not None:
            ax.set_xlabel(xlabel)
        else:
            ax.set_xlabel("Chromosome")

        if ylabel is not None:
            ax.set_ylabel(ylabel)
        else:
            ax.set_ylabel(r"$\mathrm{-log}_{10}(\mathit{p})$")

        if xspan is not None:
            ax.axvspan(*xspan, facecolor='lightgray')
        if yspan is not None:
            ax.axhspan(*yspan, facecolor='lightgray')

        ax.set_xmargin(0.01)
        ax.set_ylim(bottom=0, top=None)
        ax.tick_params(axis='x', rotation=xticks_label_rotation)
        return ax


class QQTheme(object):
    def __str__(self):
        return "QQTheme"

    __repr__ = __str__

    def __call__(self, *args, **kwargs):
        return self.apply(*args, **kwargs)

    @staticmethod
    def apply(ax: axes.Axes = None):
        if ax is None:
            ax = plt.gca()
        Theme.apply(ax)
        ax.set_xlabel(r"Expected $\mathrm{-log}_{10}(\mathit{p})$")
        ax.set_ylabel(r"Observed $\mathrm{-log}_{10}(\mathit{p})$")
        return ax

