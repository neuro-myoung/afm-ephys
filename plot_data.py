import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import collections as mc
import numpy as np


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

    def plot_sweep(self, sweep, scalebars=False):
        """
        This function will take as an argument a single sweep of pre-processed
        afm-ephys data and create a representative trace.

        The traces shown by default are position, force, work, and current. The
        traces are stacked vertically and aligned on the time axis. By default
        there are no axis or labels save for the titles indicating the data
        shown in each plot.
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
            self.plot_range = (np.max(self.plot_dat[var[1]])
                               - np.min(self.plot_dat[var[1]]))
            ax.set_ylim(np.min(self.plot_dat[var[1]]) - 0.05 * self.plot_range,
                        np.max(self.plot_dat[var[1]]) + 0.05 * self.plot_range)
            ax.axis('off')

            if scalebars is True:
                self.add_scalebars(ax, var[1])
            else:
                pass

        plt.tight_layout()
        plt.show()
        plt.savefig(self.fullpath + '_ex-trace.pdf', dpi=300,
                    transparent=True)

    def add_scalebars(self, ax, var):
        print(var)
        if var == 'i_blsub':
            [x, xend] = [np.min(self.plot_dat['ti']),
                         np.min(self.plot_dat['ti']) + 100]
            [y, yend] = [0.9 * np.max(self.plot_dat[var])-0.2*self.plot_range,
                         0.9 * np.max(self.plot_dat[var])]
            lines = [[(x, y), (x, yend)],
                     [(x, y), (xend, y)]]
        else:
            x = np.min(self.plot_dat['tin0'])
            [y, yend] = [0.9 * np.max(self.plot_dat[var])-0.4*self.plot_range,
                         0.9 * np.max(self.plot_dat[var])]
            print(y)
            print(yend)
            lines = [[(x, y), (x, yend)]]

        lc = mc.LineCollection(lines, linewidths=0.5, colors='black')
        ax.add_collection(lc)


x = PlotData('example/test', (350, 1000))
x.plot_sweep(10, scalebars=True)
