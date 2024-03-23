# -*- coding: utf-8 -*-
# @Time    : 2024/1/19 14:50
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : base.py

# basic visual class for haplot
import numpy as np
import pandas as pd
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import axes, patches


class Gene(object):
    """
    Gene class for haplot
    """

    def __init__(self, name: str, chrom: str, start: float, end: float, strand: str,
                 exons_data: str | list = None, exons_color: str | list = 'C5', exons_shape: str | list = 'rect',
                 utr_data: str | list = None, utr_color: str | list = 'lightgray', utr_shape: str | list = 'rect',
                 introns_shape: str | list = 'line'):
        """
        Init a gene object

        :param name: gene name
        :param chrom: the chromosome that gene located
        :param start: gene start position with bp unit
        :param end: gene end position with bp unit
        :param strand: gene strand, + or -
        :param exons_data: exons list of gene, format: [[start1, end1], [start2, end2], ...] or
        exons string of gene, format: "start1-end1,start2-end2,..."
        :param exons_color: exons color, default is 'C5'
        :param exons_shape: exons shape, default is 'rect', can be 'rect', 'circle', 'ellipse', 'polygon', etc.
        :param utr_data: utr list of gene, format: [[start1, end1], [start2, end2], ...] or
        utr string of gene, format: "start1-end1,start2-end2,..."
        :param utr_color: utr color, default is 'lightgray'
        :param utr_shape: utr shape, default is 'rect', can be 'rect', 'circle', 'ellipse', 'polygon', etc.
        :param introns_shape: introns shape, default is 'line', can be 'line', 'polyline', 'arrow', etc.
        """
        self.name = name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand

        if isinstance(exons_data, str):
            exons_list = [list(map(float, exon.split('-')))
                          for exon in exons_data.replace(' ', '').split(',')]
        elif isinstance(exons_data, list):
            exons_list = exons_data
        else:
            raise ValueError("exons_data of gene ({}) is none or not supported".format(name))

        if isinstance(utr_data, str):
            utr_list = [list(map(float, utr.split('-')))
                        for utr in utr_data.replace(' ', '').split(',')]
        elif isinstance(utr_data, list):
            utr_list = utr_data
        else:
            utr_list = self.get_utr_from_exons(exons_list)

        self.exons = np.asarray(exons_list)
        self.exons.sort(axis=0)
        self.exons = self.exons.astype(np.float64)

        self.utr = np.asarray(utr_list)
        self.utr.sort(axis=0)
        self.utr = self.utr.astype(np.float64)

        # set default options
        self._options = {
            'gene': {
                'color': 'C3',
                'height': 0.03,
                'zorder': 5,
                'alpha': 1
            },
            'utr': {
                'color': [],
                'shape': [],
                'height': 0.03,
                'zorder': 3,
                'alpha': 1
            },
            'exons': {
                'color': [],
                'shape': [],
                'height': 0.03,
                'zorder': 4,
                'alpha': 1
            },
            'introns': {
                'color': 'teal',
                'shape': [],
                'height': 0.01,
                'zorder': 2,
                'alpha': 1
            },
            'strand': {
                'color': 'teal',
                'zorder': 2,
                'alpha': 1
            },
        }

        # define the center and height of gene
        self.center = 0.5
        self.height = self.get_height()

        # whether to set the start position to zero
        self._zero_start = False

        # define parameters for stacked plot
        self._stacked_y0 = None
        self._stacked_y1 = None

        # set options
        self.set_exons_color(exons_color)
        self.set_exons_shape(exons_shape)
        self.set_utr_color(utr_color)
        self.set_utr_shape(utr_shape)
        self.set_introns_shape(introns_shape)

    def __str__(self):
        return "Gene({})".format(self.name)

    def __repr__(self):
        return "Gene({})".format(self.name)

    def __len__(self):
        return self.end - self.start + 1

    def get_height(self):
        """
        Get gene height

        :return: gene height
        """
        return self._options['gene']['height']

    def get_length(self):
        """
        Get gene length

        :return: gene length
        """
        return self.end - self.start + 1

    def get_introns(self):
        """
        Get introns location of gene

        :return: intron array
        """
        starts = self.exons[:-1, 1]
        ends = self.exons[1:, 0]
        return np.column_stack((starts, ends))

    def get_exons(self):
        """
        Get exons location of gene

        :return: exon array
        """
        return self.exons

    def get_utr_from_exons(self, exons: list):
        """
        Get utr location from exons

        :param exons: exons list
        :return: utr array
        """
        if self.strand == '+':
            utr = [[self.start, exons[0][0]], [exons[-1][1], self.end]]
        else:
            utr = [[self.start, exons[0][0]], [exons[-1][1], self.end]]
        return np.asarray(utr)

    def reverse_strand(self):
        """
        Reverse the direction of strand (from 3'->5' to 5'->3') to get new position

        :return: None
        """
        if self.strand == '-':
            self.strand = '+'
            self.exons = self.end - self.exons[::-1][:, ::-1]
            self.exons += self.start
            self._options['exons']['color'] = self._options['exons']['color'][::-1]

    def set_offset(self, offset_x=None, offset_y=None):
        """
        Set offset of gene

        :param offset_x: offset value in x-axis
        :param offset_y: offset value in y-axis
        :return: None
        """
        self.start += offset_x
        self.end += offset_x
        self.exons += offset_x
        self.utr += offset_x
        self.center += offset_y

    def set_zero_start(self):
        """
        Set start position to zero

        :return: None
        """
        self.set_offset(-self.start, 0)
        self._zero_start = True

    def scale(self, min_value=0, max_value=1):
        """
        Scale gene to a new range

        :param min_value: min value of new range
        :param max_value: max value of new range
        :return: None
        """
        old_range = self.end - self.start
        new_range = max_value - min_value
        self.exons = (self.exons - self.start) / old_range * new_range + min_value
        self.utr = (self.utr - self.start) / old_range * new_range + min_value
        self.start = min_value
        self.end = max_value

    def set_exons_shape(self, shape):
        """
        Set exons shape

        :param shape: shape value, can be 'rect', 'ellipse', 'polygon', etc.
        :return: None
        """
        if isinstance(shape, str):
            self._options['exons']['shape'] = [shape] * len(self.exons)
        elif isinstance(shape, list) and len(shape) == len(self.exons):
            self._options['exons']['shape'] = shape
        else:
            raise ValueError("the length of shape is not equal to exons")

    def set_utr_shape(self, shape):
        """
        Set utr shape

        :param shape: shape value, can be 'rect', 'ellipse', 'polygon', etc.
        :return: None
        """
        if isinstance(shape, str):
            self._options['utr']['shape'] = [shape] * len(self.utr)
        elif isinstance(shape, list) and len(shape) == len(self.utr):
            self._options['utr']['shape'] = shape
        else:
            raise ValueError("the length of shape is not equal to utr")

    def set_utr_color(self, color):
        """
        Set utr color

        :param color: color value
        :return: None
        """
        if isinstance(color, str):
            self._options['utr']['color'] = [color] * len(self.utr)
        elif isinstance(color, list) and len(color) == len(self.utr):
            self._options['utr']['color'] = color
        else:
            raise ValueError("the length of color is not equal to utr")

    def set_exons_color(self, color):
        """
        Set exons color

        :param color: color value
        :return: None
        """
        if isinstance(color, str):
            self._options['exons']['color'] = [color] * len(self.exons)
        elif isinstance(color, list) and len(color) == len(self.exons):
            self._options['exons']['color'] = color
        else:
            raise ValueError("the length of color is not equal to exons")

    def set_exons_height(self, height):
        """
        Set exons height

        :param height: exon height within 0-1
        :return: None
        """
        self._options['exons']['height'] = height

    def set_introns_shape(self, shape):
        """
        Set introns shape

        :param shape: shape value, can be 'line', 'polyline', 'arrow'
        :return: None
        """
        if isinstance(shape, str):
            self._options['introns']['shape'] = [shape] * len(self.get_introns())
        elif isinstance(shape, list) and len(shape) == len(self.get_introns()):
            self._options['introns']['shape'] = shape
        else:
            raise ValueError("the length of shape is not equal to exons")

    def set_introns_color(self, color):
        """
        Set introns color

        :param color: color value
        :return: None
        """
        self._options['introns']['color'] = color

    def set_introns_height(self, height):
        """
        Set introns height

        :param height: intron height within 0-1
        :return: None
        """
        self._options['introns']['height'] = height

    def set_center(self, center: float):
        """
        Set center of gene

        :param center: center value in y-axis
        :return: None
        """
        self.center = center

    def set_options(self, options: dict):
        """
        Set options of gene

        :param options: options dict
        :return: None
        """
        self._options = options

    def get_options(self):
        """
        Get options of gene

        :return: options dict
        """
        return self._options

    def draw_primers(self, ax: axes.Axes, primers: list[list[float]], options: dict = None):
        """
        Draw primers on ax

        :param ax: matplotlib ax object
        :param primers: primers list, format: [[start, len, end, len], ...].
         start and end should be relative to gene
        :param options: options dict for plot
        :return: None
        """
        primers = np.asarray(primers)
        if self._zero_start:
            primers[:, [0, 2]] -= 1
        else:
            primers[:, [0, 2]] += self.start - 1

        _options = {
            'color': 'C1',
            'height': 0.03,
            'padding': 0.03,
            'zorder': 6,
            'alpha': 1,
            'stacked': True
        }

        if options is not None:
            _options.update(options)

        y0 = self.center + self.get_height() / 2 + _options['padding']
        y1 = y0 + _options['height']

        if _options['stacked']:
            y0 = self._stacked_y0 if self._stacked_y0 is not None else y0
            y1 = self._stacked_y1 if self._stacked_y1 is not None else y1

        for i, p in enumerate(primers):
            # draw forward primer
            ax.arrow(p[0], (y0 + y1) / 2, p[1], 0, color=_options['color'], zorder=_options['zorder'],
                     alpha=_options['alpha'], head_width=0.02, head_length=10, shape='right')
            # draw reverse primer
            ax.arrow(p[2], (y0 + y1) / 2, -p[3], 0, color=_options['color'], zorder=_options['zorder'],
                        alpha=_options['alpha'], head_width=0.02, head_length=10, shape='right')
            # draw line between forward and reverse primer
            ax.plot([p[0], p[2]], [(y0 + y1) / 2, (y0 + y1) / 2], color='gray', linestyle='--')
            # annotation number between forward and reverse primer
            ax.annotate(str(i + 1), xy=((p[0] + p[2]) / 2, y0), xycoords='data',
                        xytext=(0, 5), textcoords="offset points",
                        color='black', fontsize=10, ha='center', va='bottom')
            # update y0 and y1
            y0 += _options['height'] + _options['padding']
            y1 += _options['height'] + _options['padding']

        if _options['stacked']:
            self._stacked_y0 = y0
            self._stacked_y1 = y1

    def draw_gene(self, ax: axes.Axes, options: dict = None):
        """
        Draw gene on ax

        :param ax: matplotlib ax object
        :param options: options dict for plot
        :return: None
        """
        xy = (self.start, self.center - options['gene']['height'] / 2)
        width = self.end - self.start
        height = options['gene']['height']
        color = options['gene']['color']
        zorder = options['gene']['zorder']
        alpha = options['gene']['alpha']
        ax.add_patch(patches.Rectangle(xy, width, height, color=color, zorder=zorder, alpha=alpha))

    def draw_gene_label(self, ax: axes.Axes, options: dict = None):
        """
        Draw gene label on ax

        :param ax: matplotlib ax object
        :param options: options dict for plot
        :return: None
        """
        ax.annotate(self.name, xy=(self.start, self.center), xycoords='data',
                    xytext=(-10, 0), textcoords="offset points",
                    color='black', fontsize=10, ha='right', va='center')

    def draw_utr(self, ax: axes.Axes, options: dict = None):
        """
        Draw utr on ax

        :param ax: matplotlib ax object
        :param options: options dict for plot
        :return: None
        """
        for i, utr in enumerate(self.utr):
            xy = (utr[0], self.center - options['utr']['height'] / 2)
            width = utr[1] - utr[0]
            height = options['utr']['height']
            color = options['utr']['color'][i]
            zorder = options['utr']['zorder']
            alpha = options['utr']['alpha']
            ax.add_patch(patches.Rectangle(xy, width, height, color=color, zorder=zorder, alpha=alpha))

    def draw_exons(self, ax: axes.Axes, options: dict = None):
        """
        Draw exons on ax

        :param ax: matplotlib ax object
        :param options: options dict for plot
        :return: None
        """
        for i, exon in enumerate(self.exons):
            xy = (exon[0], self.center - options['exons']['height'] / 2)
            width = exon[1] - exon[0]
            height = options['exons']['height']
            color = options['exons']['color'][i]
            zorder = options['exons']['zorder']
            alpha = options['exons']['alpha']
            if options['exons']['shape'][i] == 'rect':
                ax.add_patch(patches.Rectangle(xy, width, height, color=color, zorder=zorder, alpha=alpha))
            elif options['exons']['shape'][i] == 'ellipse':
                # center of ellipse
                xy_c = (exon[0] + width / 2, self.center)
                ax.add_patch(patches.Ellipse(xy_c, width, height, color=color, zorder=zorder, alpha=alpha))
            elif options['exons']['shape'][i] == 'fancy_box':
                ax.add_patch(patches.FancyBboxPatch(xy, width, height, color=color, zorder=zorder, alpha=alpha,
                                                    boxstyle=patches.BoxStyle("Round", pad=10),
                                                    mutation_scale=1, mutation_aspect=height / width))
            else:
                raise ValueError("shape {} is not supported".format(options['exons']['shape'][i]))

    def draw_introns(self, ax: axes.Axes, options: dict = None):
        """
        Draw introns on ax

        :param ax: matplotlib ax object
        :param options: options dict for plot
        :return: None
        """
        for i, intron in enumerate(self.get_introns()):
            color = options['introns']['color']
            zorder = options['introns']['zorder']
            alpha = options['introns']['alpha']
            if options['introns']['shape'][i] == 'polyline':
                con = patches.ConnectionPatch(xyA=(intron[0], self.center), xyB=(intron[1], self.center),
                                              coordsA="data", coordsB="data", axesA=ax, axesB=ax,
                                              color=color, zorder=zorder, alpha=alpha,
                                              connectionstyle="angle,angleA=30,angleB=150,rad=0")
                ax.add_patch(con)
            elif options['introns']['shape'][i] == 'arrow':
                if self.strand == '+':
                    ax.add_patch(patches.FancyArrowPatch((intron[0], self.center),
                                                         (intron[1], self.center),
                                                         color=color, zorder=zorder, alpha=alpha,
                                                         arrowstyle='->', mutation_scale=10,
                                                         shrinkA=0, shrinkB=0))
                else:
                    ax.add_patch(patches.FancyArrowPatch((intron[1], self.center),
                                                         (intron[0], self.center),
                                                         color=color, zorder=zorder, alpha=alpha,
                                                         arrowstyle='<-', mutation_scale=10,
                                                         shrinkA=0, shrinkB=0))
            else:
                ax.plot(intron, [self.center, self.center], color=color, zorder=zorder, alpha=alpha)

    @staticmethod
    def draw_strand(ax: axes.Axes, options: dict = None):
        """
        Draw strand on ax

        :param ax: matplotlib ax object
        :param options: options dict for plot
        :return: None
        """
        ax.annotate("5'", xy=(0, 0), xycoords='axes fraction',
                    xytext=(-10, 0), textcoords="offset points",
                    color=options['strand']['color'], fontsize=10, ha='right', va='center')
        ax.annotate("3'", xy=(1, 0), xycoords='axes fraction',
                    xytext=(10, 0), textcoords="offset points",
                    color=options['strand']['color'], fontsize=10, ha='left', va='center')
        # only keep bottom spine
        ax.spines[['left', 'top', 'right']].set_visible(False)
        # set bottom spine color and width
        ax.spines['bottom'].set_color(options['strand']['color'])
        ax.spines['bottom'].set_linewidth(2)
        # set ticks parameters
        ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
        ax.tick_params(axis='x', which='both', color=options['strand']['color'],
                       labelcolor=options['strand']['color'])
        # set x-axis ticks locators
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=5))
        ax.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(n=5))
        # add grid lines in x-axis
        ax.grid(axis='x', which='major', color='gray', alpha=0.2, linewidth=1, linestyle='--')

    def draw_bezier_links(self, ax: axes.Axes, gene: 'Gene', options: dict = None):
        """
        Draw links between two genes

        :param ax: matplotlib ax object
        :param gene: another gene object
        :param options: options dict for link
        :return: None
        """
        _options = {
            'color': 'grey',
            'fc': 'lightblue',
            'ec': 'none',
            'lw': 1,
            'zorder': 2.5,
            'padding': 0.01
        }

        if options is not None:
            _options.update(options)

        if self.center > gene.center:
            yA = self.center - self._options['exons']['height'] / 2 - _options['padding']
            yB = gene.center + gene.get_options()['exons']['height'] / 2 + _options['padding']
        else:
            yA = self.center + self._options['exons']['height'] / 2 + _options['padding']
            yB = gene.center - gene.get_options()['exons']['height'] / 2 - _options['padding']
        xyA_from = (self.start, yA)
        xyA_to = (self.end, yA)
        xyB_from = (gene.start, yB)
        xyB_to = (gene.end, yB)

        # calculate the angle (in radians) of two genes
        angleA = np.arctan2(xyA_to[1] - xyA_from[1], xyA_to[0] - xyA_from[0])
        angleB = np.arctan2(xyB_to[1] - xyB_from[1], xyB_to[0] - xyB_from[0])

        # calculate the center point of two genes
        cenAB_from = [(xyA_from[0] + xyB_from[0]) / 2, (xyA_from[1] + xyB_from[1]) / 2]
        cenAB_to = [(xyA_to[0] + xyB_to[0]) / 2, (xyA_to[1] + xyB_to[1]) / 2]

        # calculate the control point of two genes
        conA_from = [xyA_from[0] + np.tan(angleA) * (cenAB_from[1] - xyA_from[1]), cenAB_from[1]]
        conB_from = [xyB_from[0] + np.tan(angleB) * (cenAB_from[1] - xyB_from[1]), cenAB_from[1]]

        conA_to = [xyA_to[0] + np.tan(angleA) * (cenAB_to[1] - xyA_to[1]), cenAB_to[1]]
        conB_to = [xyB_to[0] + np.tan(angleB) * (cenAB_to[1] - xyB_to[1]), cenAB_to[1]]

        # draw Bézier curve
        path_data = [
            (patches.Path.MOVETO, xyA_from),
            (patches.Path.CURVE4, conA_from),
            (patches.Path.CURVE4, conB_from),
            (patches.Path.CURVE4, xyB_from),
            (patches.Path.LINETO, xyB_to),
            (patches.Path.CURVE4, conB_to),
            (patches.Path.CURVE4, conA_to),
            (patches.Path.CURVE4, xyA_to),
            (patches.Path.CLOSEPOLY, xyA_from)
        ]
        codes, vertices = zip(*path_data)
        path = patches.Path(vertices, codes)
        patch = patches.PathPatch(path, facecolor=_options['fc'], edgecolor=_options['ec'],
                                  lw=_options['lw'], zorder=_options['zorder'])
        ax.add_patch(patch)

    def plot(self, ax: axes.Axes, options: dict = None, draw_gene=False, draw_exons=True, draw_introns=True,
             draw_strand=True, draw_label=True, draw_utr=True, link_gene=None, link_options=None):
        """
        Plot gene on ax

        :param ax: matplotlib ax object
        :param options: options dict for plot
        :param draw_gene: whether plot gene
        :param draw_exons: whether plot exons
        :param draw_introns: whether plot introns
        :param draw_strand: whether plot strand
        :param draw_label: whether plot label
        :param draw_utr: whether plot utr
        :param link_gene: link gene with Bézier curve, should be a Gene object
        :param link_options: options dict for link
        :return: None
        """
        if ax is None:
            raise ValueError("ax is None")

        if options is None:
            options = self.get_options()

        ax.autoscale(enable=True, axis='both', tight=True)
        if draw_gene:
            self.draw_gene(ax, options)
            print(ax.get_ylim())
        if draw_exons:
            self.draw_exons(ax, options)
            print(ax.get_ylim())
        if draw_introns:
            self.draw_introns(ax, options)
            print(ax.get_ylim())
        if draw_strand:
            self.draw_strand(ax, options)
            print(ax.get_ylim())
        if draw_label:
            self.draw_gene_label(ax, options)
            print(ax.get_ylim())
        if draw_utr:
            self.draw_utr(ax, options)
            print(ax.get_ylim())
        if link_gene is not None:
            self.draw_bezier_links(ax, link_gene, link_options)
            print(ax.get_ylim())


class Chromosome(object):
    def __init__(self, name: str, chrom: str, start: float, end: float):
        self.name = name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.offset = [0, 0]
        self._options = {
            'color': 'lightgrey',
            'height': 0.05,
            'zorder': 1,
            'alpha': 1
        }

    def __str__(self):
        return "Chromosome ({})".format(self.name)

    def __repr__(self):
        return "Chromosome ({})".format(self.name)

    def __len__(self):
        return self.end - self.start + 1

    def get_length(self):
        return self.end - self.start + 1

    def set_x_offset(self, offset):
        self.offset[0] = offset

    def set_y_offset(self, offset):
        self.offset[1] = offset

    def set_color(self, color):
        self._options['color'] = color

    def draw_chromosome(self, ax: axes.Axes, options: dict = None):
        if options is not None:
            self._options.update(options)
        xy = (self.start + self.offset[0], 0 + self.offset[1])
        width = self.end - self.start
        ax.add_patch(patches.Rectangle(xy,
                                       width,
                                       self._options['height'],
                                       color=self._options['color'],
                                       zorder=self._options['zorder'],
                                       alpha=self._options['alpha']))

    def plot(self, ax: axes.Axes, options: dict = None):
        if ax is None:
            raise ValueError("ax is None")
        ax.autoscale(enable=True, axis='both')

        self.draw_chromosome(ax, options)


class Heatmap(object):
    """
    Heatmap class for haplot
    """

    def __init__(self, data: np.ndarray, col_labels: list[str] = None, row_labels: list[str] = None,
                 cmap: str = 'PuOr'):
        """
        Init a heatmap object

        :param data: heatmap data
        :param col_labels: column labels
        :param row_labels: row labels
        :param cmap: color map
        """
        self.data = data
        self.text = data
        self.norm = mpl.colors.Normalize(vmin=np.min(self.data), vmax=np.max(self.data))
        self.cmap = mpl.colormaps.get_cmap(cmap)
        if row_labels is None:
            self.row_labels = np.arange(self.data.shape[0])
        else:
            self.row_labels = row_labels
        if col_labels is None:
            self.col_labels = np.arange(self.data.shape[1])
        else:
            self.col_labels = col_labels
        self.ax_divider = None

    def __str__(self):
        return "Heatmap"

    def __repr__(self):
        return "Heatmap"

    def set_cmap(self, cmap: str):
        self.cmap = mpl.colormaps.get_cmap(cmap)

    def set_row_labels(self, row_labels: list[str]):
        if len(row_labels) != len(self.row_labels):
            raise ValueError("the length of row_labels is not equal to data")
        self.row_labels = row_labels

    def set_col_labels(self, col_labels: list[str]):
        if len(col_labels) != len(self.col_labels):
            raise ValueError("the length of col_labels is not equal to data")
        self.col_labels = col_labels

    def draw_row_labels(self, ax: axes.Axes, cmap: str = 'RdBu', labels_as_color: bool = False,
                        show_labels: bool = True, show_ticks: bool = True):
        """
        Draw row labels on ax

        :param ax: matplotlib ax object
        :param cmap: color map
        :param labels_as_color: whether to use labels as color
        :param show_labels: whether show labels
        :param show_ticks: whether show ticks
        """
        if labels_as_color:
            _ax = self.ax_divider.append_axes("left", size="5%", pad="1%", sharey=ax)
            codes, uniques = pd.Series(self.row_labels).factorize()
            _ax.imshow(codes.reshape(-1, 1),
                       cmap=mpl.colormaps.get_cmap(cmap).resampled(len(uniques)),
                       aspect='auto')
            _ax.grid(which='minor',
                     linestyle='-',
                     color='w',
                     linewidth=2)
            _ax.spines[:].set_visible(False)
            _ax.tick_params(which='minor', bottom=False, left=False)
            _ax.tick_params(which='major', bottom=False, left=show_ticks, labelbottom=False, labelleft=show_labels)
            ax.tick_params(axis='y', which='both', left=False, labelleft=False)
            ax.set_yticks(np.arange(len(self.row_labels)), labels=self.row_labels)
        else:
            ax.set_yticks(np.arange(len(self.row_labels)), labels=self.row_labels)

    def draw_col_labels(self, ax: axes.Axes, cmap: str = 'RdBu', labels_as_color: bool = False,
                        show_labels: bool = True, show_ticks: bool = True):
        """
        Draw column labels on ax

        :param ax: matplotlib ax object
        :param cmap: color map
        :param labels_as_color: whether to use labels as color
        :param show_labels: whether show labels
        :param show_ticks: whether show ticks
        """
        if labels_as_color:
            _ax = self.ax_divider.append_axes("top", size="5%", pad="1%", sharex=ax)
            codes, uniques = pd.Series(self.col_labels).factorize()
            _ax.imshow(codes.reshape(1, -1),
                       cmap=mpl.colormaps.get_cmap(cmap).resampled(len(uniques)),
                       aspect='auto')
            _ax.grid(which='minor',
                     linestyle='-',
                     color='w',
                     linewidth=2)
            _ax.spines[:].set_visible(False)
            _ax.tick_params(which='minor', bottom=False, left=False)
            _ax.tick_params(which='major', bottom=False, left=False, labelbottom=False, labelleft=False,
                            top=show_ticks, labeltop=show_labels)
            _ax.set_xticks(np.arange(len(self.col_labels)), labels=self.col_labels, rotation=30, ha='left',
                           rotation_mode='anchor')
            ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
        else:
            ax.set_xticks(np.arange(len(self.col_labels)), labels=self.col_labels)

    def draw_row_dendrogram(self, ax: axes.Axes):
        pass

    def draw_col_dendrogram(self, ax: axes.Axes):
        pass

    def draw_text(self, ax: axes.Axes, fmt: str = '.2f'):
        """
        Draw text on ax

        :param ax: matplotlib ax object
        :param fmt: text format
        """
        # text color threshold
        threshold = self.norm(np.max(self.data)) / 2
        text_colors = ['black', 'white']
        for i in range(self.data.shape[0]):
            for j in range(self.data.shape[1]):
                ax.text(j, i, format(self.data[i, j], fmt), ha='center', va='center',
                        color=text_colors[int(self.norm(self.data[i, j]) > threshold)])

    def draw_cbar(self, ax: axes.Axes, bounds: list[float] = None):
        """
        Draw color bar on ax

        :param ax: matplotlib ax object
        :param bounds: color bar bounds
        """
        if bounds is not None:
            cax = ax.inset_axes(bounds)
        else:
            cax = ax.inset_axes([1.04, 0.2, 0.03, 0.6])
        ax.figure.colorbar(mpl.cm.ScalarMappable(norm=self.norm, cmap=self.cmap), ax=ax, cax=cax,
                           orientation='vertical')

    def draw_grid(self, ax: axes.Axes):
        """
        Draw grid on ax

        :param ax: matplotlib ax object
        """
        ax.set_yticks(np.arange(self.data.shape[0] + 1) - .5, minor=True)
        ax.set_xticks(np.arange(self.data.shape[1] + 1) - .5, minor=True)
        ax.grid(which='minor',
                linestyle='-',
                color='w',
                linewidth=2)
        ax.tick_params(which='minor',
                       bottom=False,
                       left=False)

    def draw_gene(self, ax: axes.Axes, genes: pd.DataFrame):
        """
        Draw gene on ax

        :param ax: matplotlib ax object
        :param genes: a dataFrame with gene structure data
        """
        pass

    def draw_data(self, ax: axes.Axes):
        """
        Draw heatmap data

        :param ax: matplotlib ax object
        """
        ax.imshow(self.data, aspect='auto', cmap=self.cmap, norm=self.norm)

    def plot(self, ax: axes.Axes, draw_cbar=True, draw_grid=True,
             draw_row_labels=True, draw_col_labels=True, draw_tree=False, draw_text=True):
        """
        Plot heatmap on ax

        :param ax: matplotlib ax object
        :param draw_cbar: whether draw color bar
        :param draw_grid: whether draw grid
        :param draw_row_labels: whether draw row labels
        :param draw_col_labels: whether draw column labels
        :param draw_tree: whether draw tree
        :param draw_text: whether draw text
        """
        if ax is None:
            raise ValueError("ax is None")

        self.ax_divider = make_axes_locatable(ax)
        if draw_cbar:
            self.draw_cbar(ax)
        if draw_grid:
            self.draw_grid(ax)
        if draw_row_labels and draw_tree:
            self.draw_row_labels(ax, labels_as_color=True, show_labels=False, show_ticks=False)
        else:
            self.draw_row_labels(ax, labels_as_color=True)
        if draw_col_labels:
            self.draw_col_labels(ax, labels_as_color=True)
        # must draw data and text in the end
        if draw_text:
            self.draw_text(ax)
        self.draw_data(ax)
        ax.spines[:].set_visible(False)
