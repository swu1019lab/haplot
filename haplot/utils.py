# -*- coding: utf-8 -*-
# @Time    : 2023/8/1 11:54
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : utils.py
from abc import ABC

import numpy as np
from matplotlib.offsetbox import AnchoredOffsetbox, DrawingArea
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
