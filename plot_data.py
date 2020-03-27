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

    def __init__(self, file_path):
        self.fullpath = file_path
        self.dat = pd.read_hdf(self.fullpath + '_augmented.h5')
        self.params = pd.read_csv(self.fullpath + '_params.csv')

        # Dictionaries for axis customization based on variable plotted.
        self.lab_dict = {'i_blsub': ' pA', 'force': ' nN',
                         'work': ' fJ', 'position': ' nm',
                         'ti': 'ms', 'tin0': 'ms', 'tz': 'ms'}
        self.title_dict = {'i_blsub': 'Current', 'force': 'Force',
                           'work': 'Work', 'position': 'Position',
                           'ti': 'Time', 'tin0': 'Time', 'tz': 'Time'}

        self.col_dict = {'i_blsub': 'k', 'force': 'b', 'work': 'm',
                         'position': 'g'}

        self.horiz_dict = {'i_blsub': 'ti', 'force': 'tin0',
                           'work': 'tin0', 'position': 'tz'}
        self.height_dict = {'i_blsub': 5, 'force': 1.5,
                            'work': 1.5, 'position': 1.5}

    def plot_sweep(self, sweep, vars, roi=None, scalebars=False,
                   scalelabs=False):
        """
        Arguments:
            sweep: identity of sweep number to be plotted
            vars: trace variables to be shown (as iterable)
            roi: region of interest to be plotted in ms (default: None)
            scalebars: logical as to whether or not to add automated add
                       scalebars automatically
            scalelabs: logical as to whether scalebars should be automatically
                       labeled
        Output:
            A high res plot of the specified traces aligned vertically in time
            saved as a .pdf file sharing the prefix of the file from which the
            data was derived.

        """
        if roi is not None:
            try:
                iter(roi)
            except TypeError:
                print('If an roi is give it must be an iterable!')

            self.dat_sub = self.dat[(self.dat['ti'] >= roi[0])
                                    & (self.dat['ti'] <= roi[1])]
            self.plot_dat = self.dat_sub.groupby('sweep').get_group(sweep)

        else:
            self.plot_dat = self.grps.get_group(sweep)

        colors = [self.col_dict[x] for x in vars]
        titles = [self.title_dict[x] for x in vars]
        xvals = [self.horiz_dict[x] for x in vars]
        heights = [self.height_dict[x] for x in vars]

        self.fig, self.axs = plt.subplots(len(vars), dpi=300,
                                          figsize=(1.5, 1.5*len(vars)),
                                          gridspec_kw={'height_ratios': heights})

        for ax, x, y, color, title in zip(self.axs, xvals, vars, colors, titles):
            ax.plot(self.plot_dat[x], self.plot_dat[y],
                    color=color, linewidth=0.5)
            ax.set_title(title, size=8)
            self.plot_range = ax.get_ylim()[1] - ax.get_ylim()[0]
            self.plot_domain = ax.get_xlim()[1] - ax.get_xlim()[0]
            ax.set_ylim(np.min(self.dat_sub[y]) - 0.05 * self.plot_range,
                        np.max(self.dat_sub[y]) + 0.05 * self.plot_range)
            ax.axis('off')

            if scalebars is True:
                self.add_scalebars(ax, y)
            else:
                pass

            if scalelabs is True:
                self.add_scalelabs(ax, y)

        plt.tight_layout()
        plt.savefig(self.fullpath + '_ex-trace.pdf', dpi=300)
        plt.show()

    def add_scalebars(self, ax, var):
        """
        This function adds scalebars to the trace axis when used.

        Arguments:
            ax: plot axis to have scalebar added
            var: variable corresponding to the specific axis

        Returns:
            automatically scaled and positions scalebars on the plot axis
        """

        ylen = self.round_1_sf(0.2*self.plot_range)

        if var == 'i_blsub':
            [x, xend] = [np.min(self.dat_sub['ti']),
                         np.min(self.dat_sub['ti']) + 100]

            if max(self.plot_dat['absi_blsub']) == max(self.plot_dat[var]):
                [y, yend] = [0.9 * np.max(self.dat_sub[var])-ylen,
                             0.9 * np.max(self.dat_sub[var])]
            else:
                [y, yend] = [0.9 * np.min(self.dat_sub[var]),
                             0.9 * np.min(self.dat_sub[var]) + ylen]
                print(y)
                print(yend)

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

        Arguments:
            ax: plot axis to have scalebar added
            var: variable corresponding to the specific axis

        Returns:
            Automatically positioned labels of appropriately representing the
            scalebars present.
        """
        ylen = self.round_1_sf(0.2*self.plot_range)

        if var == 'i_blsub':
            x1 = np.min(self.dat_sub['ti']) + 0.02 * self.plot_domain
            x2 = np.min(self.dat_sub['ti']) + 0.02 * self.plot_domain

            y1 = (0.9 * np.min(self.dat_sub['i_blsub'])
                  - (0.05*self.plot_range))
            y2 = (0.9 * np.min(self.dat_sub['i_blsub'])
                  + (0.02*self.plot_range))

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
