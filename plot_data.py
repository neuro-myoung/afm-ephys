import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import collections as mc
import numpy as np
from math import log10, floor


class PlotData(object):
    """
    This class takes as input a file path, filename, and region of interest to
    make a subsetted afm-ephys data object that can be plotted using the
    defined methods.
    """

    def __init__(self, file_path, rois):
        self.fullpath = file_path
        self.start = rois[0]
        self.end = rois[1]
        self.dat = pd.read_hdf(self.fullpath + '_augmented.h5')
        self.params = pd.read_csv(self.fullpath + '_params.csv')

        self.grps = self.dat.groupby('sweep')
        self.dat_sub = self.dat[((self.dat['ti'] >= self.start)
                                 & (self.dat['ti'] <= self.end))]
        self.grps_sub = self.dat_sub.groupby('sweep')

        self.lab_dict = {'i_blsub': ' pA', 'force': ' nN',
                         'work': ' fJ', 'position': ' nm',
                         'ti': 'ms', 'tin0': 'ms', 'tz': 'ms'}

    def plot_sweep(self, sweep, scalebars=False, scalelabs=False):
        """
        This function will take as an argument a single sweep of pre-processed
        afm-ephys data and create a representative trace.

        The traces shown by default are position, force, work, and current. The
        traces are stacked vertically and aligned on the time axis. By default
        there are no axis or labels save for the titles indicating the data
        shown in each plot.
        By setting scalebars = True and/or scalelabs =True scalebars are and
        labels are sized automatically and added to the plots.
        """
        self.plot_dat = self.grps_sub.get_group(sweep)
        self.fig, self.axs = plt.subplots(nrows=4, dpi=300, figsize=(2, 4),
                                          gridspec_kw={'height_ratios': [1.5, 1.5, 1.5, 5]})
        colors = ['g', 'b', 'm', 'k']
        titles = ['Position', 'Force', 'Work', 'Current']
        vars = [('tz', 'position'), ('tin0', 'force'),
                ('tin0', 'work'), ('ti', 'i_blsub')]

        for ax, color, var, title in zip(self.axs, colors, vars, titles):
            ax.plot(self.plot_dat[var[0]], self.plot_dat[var[1]],
                    color=color, linewidth=0.5)
            ax.set_title(title, size=10)
            self.plot_range = ax.get_ylim[1] - ax.get_ylim[0]
            self.plot_domain = ax.get_xlim[1] - ax.get_xlim[0]
            ax.set_ylim(np.min(self.dat_sub[var[1]]) - 0.05 * self.plot_range,
                        np.max(self.dat_sub[var[1]]) + 0.05 * self.plot_range)
            ax.axis('off')

            if scalebars is True:
                self.add_scalebars(ax, var[1])
            else:
                pass

            if scalelabs is True:
                self.add_scalelabs(ax, var[1])

        plt.tight_layout()
        plt.savefig(self.fullpath + '_ex-trace.pdf', dpi=300)
        plt.show()

    def add_scalebars(self, ax, var):
        """
        This function adds scalebars to the trace axis when used.

        By default it only adds a scalebar to the y-axis unless it is a current
        vs. time plot. Scales are sized automatically.
        """

        ylen = self.round_1_sf(0.2*self.plot_range)

        if var == 'i_blsub':
            [x, xend] = [np.min(self.dat_sub['ti']),
                         np.min(self.dat_sub['ti']) + 100]

            [y, yend] = [0.9 * np.max(self.dat_sub[var])-ylen,
                         0.9 * np.max(self.dat_sub[var])]
            lines = [[(x, y), (x, yend)],
                     [(x, y), (xend, y)]]
        else:
            x = np.min(self.dat_sub['tin0'])
            [y, yend] = [0.9 * np.max(self.dat_sub[var])-ylen*2,
                         0.9 * np.max(self.dat_sub[var])]

            lines = [[(x, y), (x, yend)]]
        lc = mc.LineCollection(lines, linewidths=0.5, colors='black')
        ax.add_collection(lc)

    def add_scalelabs(self, ax, var):
        """
        This function adds labels to the scalebars when called.
        """
        ylen = self.round_1_sf(0.2*self.plot_range)

        if var == 'i_blsub':
            x1 = np.min(self.dat_sub['ti']) + 0.02 * self.plot_domain
            y1 = 0.9 * np.max(self.dat_sub[var])-(ylen + 0.08*self.plot_range)
            x2 = np.min(self.dat_sub['ti']) + 0.02 * self.plot_domain
            y2 = 0.9 * np.max(self.dat_sub[var])-(ylen-0.02*self.plot_range)
            ax.text(x1, y1, '100 ms', fontsize=6)
            ax.text(x2, y2, str(int(ylen)) + self.lab_dict[var], fontsize=6)
        else:
            ylen *= 2
            x1 = np.min(self.dat_sub['tin0']) + 0.02 * self.plot_domain
            y1 = 0.9 * np.max(self.dat_sub[var])-(0.75*ylen)
            ax.text(x1, y1, str(int(ylen)) + self.lab_dict[var], fontsize=6)

    def round_1_sf(self, num):
        """
        This function will round a number to only 1 significant figure.
        """
        return(round(num, -int(floor(log10(abs(num))))))
