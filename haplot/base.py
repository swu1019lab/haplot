# -*- coding: utf-8 -*-
# @Time    : 2024/1/19 14:50
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : base.py

# basic visual class for haplot
import numpy as np
from matplotlib import axes, patches


class Gene(object):
    """
    Gene class for haplot
    """

    def __init__(self, name: str, chrom: str, start: float, end: float, strand: str,
                 exons_list: list = None, exons_color: str | list = 'C5', exons_shape: str | list = 'rect'):
        """
        Init a gene object

        :param name: gene name
        :param chrom: the chromosome that gene located
        :param start: gene start position with bp unit
        :param end: gene end position with bp unit
        :param strand: gene strand, + or -
        :param exons_list: exons list of gene, each exon contains start, end position
        :param exons_color: exons color, default is 'C5'
        :param exons_shape: exons shape, default is 'rect', can be 'rect', 'circle', 'ellipse', 'polygon', etc.
        """
        self.name = name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.exons = np.asarray(exons_list)
        self.exons.sort(axis=0)
        self.exons = self.exons.astype(np.float64)

        # set default options
        self.options = {
            'center': 0.5,
            'gene': {
                'color': 'red',
                'height': 0.10,
                'zorder': 3,
                'alpha': 1
            },
            'exons': {
                'color': [],
                'shape': [],
                'height': 0.10,
                'zorder': 2,
                'alpha': 1
            },
            'introns': {
                'color': 'C0',
                'height': 0.01,
                'zorder': 1,
                'alpha': 1
            },
            'strand': {
                'color': 'C0',
                'zorder': 1,
                'alpha': 1
            },
        }
        self.set_exons_color(exons_color)
        self.set_exons_shape(exons_shape)

    def __str__(self):
        return "Gene({})".format(self.name)

    def __repr__(self):
        return "Gene({})".format(self.name)

    def __len__(self):
        return self.end - self.start + 1

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
        intron_list = []
        for i in range(len(self.exons) - 1):
            intron_list.append([self.exons[i][1], self.exons[i + 1][0]])
        return np.asarray(intron_list)

    def get_exons(self):
        """
        Get exons location of gene

        :return: exon array
        """
        return self.exons

    def set_offset(self, offset):
        """
        Set offset of gene

        :param offset: offset value in x axis
        :return: None
        """
        self.start += offset
        self.end += offset
        for i in range(len(self.exons)):
            self.exons[i][0] += offset
            self.exons[i][1] += offset

    def set_zero_start(self):
        """
        Set start position to zero

        :return: None
        """
        self.set_offset(-self.start)

    def scale(self, min_value=0, max_value=1):
        """
        Scale gene to a new range

        :param min_value: min value of new range
        :param max_value: max value of new range
        :return: None
        """
        old_range = self.end - self.start
        new_range = max_value - min_value
        new_array = np.zeros_like(self.exons, dtype=np.float64)
        for i in range(len(self.exons)):
            new_array[i][0] = (self.exons[i][0] - self.start) / old_range * new_range + min_value
            new_array[i][1] = (self.exons[i][1] - self.start) / old_range * new_range + min_value
        self.exons = new_array
        self.start = min_value
        self.end = max_value

    def set_exons_shape(self, shape):
        """
        Set exons shape

        :param shape: shape value, can be 'rect', 'ellipse', 'polygon', etc.
        :return: None
        """
        if isinstance(shape, str):
            self.options['exons']['shape'] = [shape] * len(self.exons)
        elif isinstance(shape, list) and len(shape) == len(self.exons):
            self.options['exons']['shape'] = shape
        else:
            raise ValueError("the length of shape is not equal to exons")

    def set_exons_color(self, color):
        """
        Set exons color

        :param color: color value
        :return: None
        """
        if isinstance(color, str):
            self.options['exons']['color'] = [color] * len(self.exons)
        elif isinstance(color, list) and len(color) == len(self.exons):
            self.options['exons']['color'] = color
        else:
            raise ValueError("the length of color is not equal to exons")

    def set_exons_height(self, height):
        """
        Set exons height

        :param height: exon height within 0-1
        :return: None
        """
        self.options['exons']['height'] = height

    def set_introns_color(self, color):
        """
        Set introns color

        :param color: color value
        :return: None
        """
        self.options['introns']['color'] = color

    def set_center(self, center: float):
        """
        Set center of gene

        :param center: center value range from 0 to 1
        :return: None
        """
        self.options['center'] = center

    def set_options(self, options: dict):
        """
        Set options of gene

        :param options: options dict
        :return: None
        """
        self.options = options

    def get_options(self):
        """
        Get options of gene

        :return: options dict
        """
        return self.options

    def draw_gene(self, ax: axes.Axes, options: dict = None):
        """
        Draw gene on ax

        :param ax: matplotlib ax object
        :param options: options dict for plot
        :return: None
        """
        xy = (self.start, options['center'] - options['gene']['height'] / 2)
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
        ax.annotate(self.name, xy=(self.start, options['center']), xycoords='data',
                    xytext=(-10, 0), textcoords="offset points",
                    color='black', fontsize=10, ha='right', va='center')

    def draw_exons(self, ax: axes.Axes, options: dict = None):
        """
        Draw exons on ax

        :param ax: matplotlib ax object
        :param options: options dict for plot
        :return: None
        """
        for i, exon in enumerate(self.exons):
            xy = (exon[0], options['center'] - options['exons']['height'] / 2)
            width = exon[1] - exon[0]
            height = options['exons']['height']
            color = options['exons']['color'][i]
            zorder = options['exons']['zorder']
            alpha = options['exons']['alpha']
            if options['exons']['shape'][i] == 'rect':
                ax.add_patch(patches.Rectangle(xy, width, height, color=color, zorder=zorder, alpha=alpha))
            elif options['exons']['shape'][i] == 'ellipse':
                # center of ellipse
                xy_c = (exon[0] + width / 2, options['center'])
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
        for intron in self.get_introns():
            xy = (intron[0], options['center'] - options['introns']['height'] / 2)
            width = intron[1] - intron[0]
            height = options['introns']['height']
            color = options['introns']['color']
            zorder = options['introns']['zorder']
            alpha = options['introns']['alpha']
            ax.add_patch(patches.Rectangle(xy, width, height, color=color, zorder=zorder, alpha=alpha))

    def draw_strand(self, ax: axes.Axes, options: dict = None):
        """
        Draw strand on ax

        :param ax: matplotlib ax object
        :param options: options dict for plot
        :return: None
        """
        if self.strand == '+':
            arrow_style = '->'
            posA = (self.start, options['center'])
            posB = (self.end, options['center'])
        else:
            arrow_style = '<-'
            posA = (self.start, options['center'])
            posB = (self.end, options['center'])
        ax.add_patch(patches.FancyArrowPatch(posA, posB, arrowstyle=arrow_style,
                                             color=options['strand']['color'],
                                             zorder=options['strand']['zorder'],
                                             alpha=options['strand']['alpha'],
                                             mutation_scale=10,
                                             shrinkA=0, shrinkB=0))

    def plot(self, ax: axes.Axes, options: dict = None, draw_gene=False, draw_exons=True, draw_introns=True,
             draw_strand=True):
        """
        Plot gene on ax

        :param ax: matplotlib ax object
        :param options: options dict for plot
        :param draw_gene: whether plot gene
        :param draw_exons: whether plot exons
        :param draw_introns: whether plot introns
        :param draw_strand: whether plot strand
        :return: None
        """
        if ax is None:
            raise ValueError("ax is None")

        if options is None:
            options = self.options

        ax.autoscale(enable=True, axis='both')
        if draw_gene:
            self.draw_gene(ax, options)
        if draw_exons:
            self.draw_exons(ax, options)
        if draw_introns:
            self.draw_introns(ax, options)
        if draw_strand:
            self.draw_strand(ax, options)
        self.draw_gene_label(ax, options)

        # only keep bottom spine
        ax.spines[['left', 'top', 'right']].set_visible(False)

        # hide ticks
        ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
        ax.tick_params(axis='x', which='both', width=2, labelsize='large')

        # add grid lines in x-axis
        ax.grid(axis='x', which='major', color='gray', alpha=0.2, linewidth=1, linestyle='--')
