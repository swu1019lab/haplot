# -*- coding: utf-8 -*-
# @Time    : 2023/8/1 11:54
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : utils.py


def auto_text_size(text, patch, ax):
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
