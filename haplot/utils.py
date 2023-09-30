# -*- coding: utf-8 -*-
# @Time    : 2023/8/1 11:54
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : utils.py
from abc import ABC

import numpy as np
from matplotlib.offsetbox import AnchoredOffsetbox, DrawingArea
from mpl_toolkits.axes_grid1.inset_locator import (
    BboxPatch, BboxConnector, BboxConnectorPatch)
from matplotlib.transforms import Bbox, TransformedBbox
from matplotlib.lines import Line2D
from matplotlib.text import Text


def adjust_text(texts, ax):
    """
    Automatically adjust text position to avoid overlap.

    Parameters
    ----------
    :param texts: a list of Text object
    :param ax: Axes object
    :return: None
    """
    return


def text_fit_patch(text, patch, ax):
    """
    Automatically adjust the font size of text to fit the patch.

    Parameters
    ----------
    :param text: Text object
    :param patch: Patch object
    :param ax: Axes object
    :return: None
    """
    font_size = text.get_fontsize()
    patch_width = patch.get_window_extent(ax.figure.canvas.get_renderer()).width
    text_width = text.get_window_extent(ax.figure.canvas.get_renderer()).width
    while text_width > patch_width and font_size > 1:
        font_size -= 1
        text.set_fontsize(font_size)


class AnchoredSizeLegend(AnchoredOffsetbox, ABC):
    def __init__(self,
                 size,
                 label,
                 label_size=10,
                 *args,
                 **kwargs):
        """
        Draw a size legend with given size in points.

        Parameters
        ----------
        :param size: a list of size in points
        :param label: a list of label
        :param label_size: font size of label
        :param args: other args of AnchoredOffsetbox
        :param kwargs: other kwargs of AnchoredOffsetbox
        """
        if isinstance(size, int):
            size = np.sqrt([size])
        else:
            size = np.sqrt(size)

        if isinstance(label, str):
            label = [label]

        # self.box = AuxTransformBox(transform)
        self.box = DrawingArea(np.max(size) * 2.2, np.max(size) * 1.2, 0, 0)

        for i in range(len(size)):
            x, y = np.max(size) / 2, size[i] / 2
            self.box.add_artist(
                Line2D([x], [y], marker='o', color='w',
                       markerfacecolor='none',
                       markeredgecolor='k',
                       markersize=size[i])
            )
            self.box.add_artist(
                Line2D([x, x + size[i]], [size[i], size[i]], color="k", linewidth=1)
            )
            self.box.add_artist(
                Text(x + size[i], size[i],
                     label[i], size=label_size,
                     color='k', ha="left", va="center")
            )
        super().__init__(*args, child=self.box, **kwargs)


def connect_bbox(ax1, ax2,
                 x_min1, x_max1, x_min2, x_max2,
                 plot_bbox1=True, plot_bbox2=True,
                 loc1a=3, loc2a=2, loc1b=4, loc2b=1,
                 prop_lines=None, prop_patches=None):
    """
    Connect two axes with bbox.

    Parameters
    ----------
    :param ax1: an Axes object
    :param ax2: an Axes object
    :param x_min1: x min of ax1 with data coordinates
    :param x_max1: x max of ax1 with data coordinates
    :param x_min2: x min of ax2 with data coordinates
    :param x_max2: x max of ax2 with data coordinates
    :param plot_bbox1: whether to plot bbox1
    :param plot_bbox2: whether to plot bbox2
    :param loc1a: loc 1a of BboxConnector
    :param loc2a: loc 2a of BboxConnector
    :param loc1b: loc 1b of BboxConnector
    :param loc2b: loc 2b of BboxConnector
    :param prop_lines: properties of lines
    :param prop_patches: properties of patches
    :return: None
    """
    if prop_lines is None:
        prop_lines = {
            "color": "gray",
        }

    if prop_patches is None:
        prop_patches = {
            **prop_lines,
            "alpha": prop_lines.get("alpha", 1) * 0.2,
            "clip_on": False,
        }

    bbox1 = TransformedBbox(Bbox.from_extents(x_min1, 0, x_max1, 1), ax1.get_xaxis_transform())
    bbox2 = TransformedBbox(Bbox.from_extents(x_min2, 0, x_max2, 1), ax2.get_xaxis_transform())
    c1 = BboxConnector(
        bbox1, bbox2, loc1=loc1a, loc2=loc2a, clip_on=False, **prop_lines)
    c2 = BboxConnector(
        bbox1, bbox2, loc1=loc1b, loc2=loc2b, clip_on=False, **prop_lines)

    bbox_patch1 = BboxPatch(bbox1, **prop_patches)
    bbox_patch2 = BboxPatch(bbox2, **prop_patches)

    p = BboxConnectorPatch(bbox1, bbox2,
                           loc1a=loc1a, loc2a=loc2a, loc1b=loc1b, loc2b=loc2b,
                           **prop_patches)

    if plot_bbox1:
        ax1.add_patch(bbox_patch1)
    if plot_bbox2:
        ax2.add_patch(bbox_patch2)

    ax2.add_patch(c1)
    ax2.add_patch(c2)
    ax2.add_patch(p)
