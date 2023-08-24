# -*- coding: utf-8 -*-
# @Time    : 2023/8/1 11:54
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : utils.py


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
