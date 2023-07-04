# -*- coding: utf-8 -*-
# @Time    : 2023/6/28 10:56
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : theme.py


from matplotlib import axes
import matplotlib.pyplot as plt
from itertools import cycle


class Theme(object):
    def __str__(self):
        return "Theme"

    __repr__ = __str__

    def __call__(self, *args, **kwargs):
        return self.apply(*args, **kwargs)

    @staticmethod
    def apply(ax: axes.Axes):
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


class BarTheme(object):
    def __init__(self):
        return


class LineTheme(object):
    def __init__(self):
        return
