import matplotlib.gridspec as gridspec
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import collections as mc
import numpy as np
from math import log10, floor


class PlotData(object):
    """
    This class takes as input a file path, filename, and region of interest to make a subsetted afm-ephys data object
    that can be plotted using the defined methods.
    """

    def __init__(self, file_path):
        self.fullpath = file_path
        self.dat = pd.read_hdf(self.fullpath + '_augmented.h5')
        self.dat['position'] /= 1000
        self.params = pd.read_csv(self.fullpath + '_params.csv')
        self.grps = self.dat.groupby('sweep')

        # Dictionaries for axis customization based on variable plotted.
        self.lab_dict = {'i_blsub': ' pA', 'force': ' nN', 'work': ' fJ', 'position': ' um', 'ti': 'ms', 'tin0': 'ms',
                         'tz': 'ms'}
        self.title_dict = {'i_blsub': 'Current', 'force': 'Force', 'work': 'Work', 'position': 'Position',
                           'ti': 'Time', 'tin0': 'Time', 'tz': 'Time'}
        self.col_dict = {'i_blsub': 'k', 'force': 'b', 'work': 'm', 'position': 'g'}
        self.horiz_dict = {'i_blsub': 'ti', 'force': 'tin0', 'work': 'tin0', 'position': 'tz'}
        self.height_dict = {'i_blsub': 3, 'force': 1.5, 'work': 1.5, 'position': 1.5}

    def plot_sweep(self, sweep, vars, roi=None, scalebars=False, scalelabs=False):
        """
        This function will return a vertically stacked plot of traces in a single sweep.

        Arguments:
            sweep: identity of sweep number to be plotted
            vars: trace variables to be shown (as iterable)
            roi: region of interest to be plotted in ms (default: None)
            scalebars: logical as to whether or not to add automated add scalebars automatically
            scalelabs: logical as to whether scalebars should be automatically labeled

        Returns:
            A high res plot of the specified traces aligned vertically in time saved as a .pdf file sharing the prefix
            of the file from which the data was derived.
        """
        if roi is not None:
            try:
                iter(roi)
            except TypeError:
                print('If an roi is give it must be an iterable!')

            self.dat_sub = self.dat[(self.dat['ti'] >= roi[0]) & (self.dat['ti'] <= roi[1])]
            self.plot_dat = self.dat_sub.groupby('sweep').get_group(sweep)

        else:
            self.plot_dat = self.grps.get_group(sweep)

        colors = [self.col_dict[x] for x in vars]
        xvals = [self.horiz_dict[x] for x in vars]
        heights = [self.height_dict[x] for x in vars]

        self.fig, self.axs = plt.subplots(len(vars), dpi=300, figsize=(1.5, 0.75*len(vars)),
                                          gridspec_kw={'height_ratios': heights})

        for ax, x, y, color in zip(self.axs, xvals, vars, colors):
            ax.plot(self.plot_dat[x], self.plot_dat[y], color=color, linewidth=0.5)
            self.plot_range = ax.get_ylim()[1] - ax.get_ylim()[0]
            self.plot_domain = ax.get_xlim()[1] - ax.get_xlim()[0]
            ax.set_ylim(np.min(self.dat[y]) - 0.05 * self.plot_range,
                        np.max(self.dat[y]) + 0.05 * self.plot_range)
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
            Automatically scaled and positioned scalebars on the plot axis
        """
        ylen = self.round_1_sf(0.2*self.plot_range)

        if var == 'i_blsub':
            [x, xend] = [np.min(self.dat_sub['ti']), np.min(self.dat_sub['ti']) + 100]

            if max(self.plot_dat['absi_blsub']) == max(self.plot_dat[var]):
                [y, yend] = [0.9 * np.max(self.dat_sub[var])-ylen, 0.9 * np.max(self.dat_sub[var])]
            else:
                [y, yend] = [0.9 * np.min(self.dat_sub[var]), 0.9 *
                             np.min(self.dat_sub[var]) + ylen]
            lines = [[(x, y), (x, yend)], [(x, y), (xend, y)]]
        else:
            x = np.min(self.dat_sub['tin0'])
            [y, yend] = [0.9 * np.max(self.dat_sub[var])-ylen*2, 0.9 * np.max(self.dat_sub[var])]

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
            Automatically positioned labels of appropriately representing the scalebars present.
        """
        ylen = self.round_1_sf(0.2*self.plot_range)

        if var == 'i_blsub':
            x1 = np.min(self.dat_sub['ti']) + 0.02 * self.plot_domain
            x2 = np.min(self.dat_sub['ti']) + 0.02 * self.plot_domain

            y1 = (0.9 * np.min(self.dat_sub['i_blsub']) - (0.1*self.plot_range))
            y2 = (0.9 * np.min(self.dat_sub['i_blsub']) + (0.02*self.plot_range))

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

    def remove_sweep(self, sweeps):
        """
        This function will remove sweeps from a dataframe by sweep number.

        Arguments:
            sweep: The number of the sweep or list of numbers of sweeps to be removed.

        Returns:
            A dataframe with the passed sweeps removed.
        """
        if hasattr(sweeps, '__iter__') is True:
            self.dat = self.dat[self.dat['sweep'] not in sweeps]
            self.grps = self.dat.groupby('sweep')
        else:
            self.dat = self.dat[self.dat['sweep'] != sweeps]
            self.grps = self.dat.groupby('sweep')

    def plot_all_sweeps(self, vars, roi=None, scalebars=False):
        """
        This function will plot all the sweeps in a given experiment.

        Arguments:
            vars: trace variables to be shown (as iterable)
            roi: region of interest to be plotted in ms (default: None)
            scalebars: logical as to whether or not to add automated add scalebars automatically to the
                       first sweep.

        Returns:
            A high res plot of the specified traces for all sweeps aligned vertically in time within each sweep
            saved as a .pdf file sharing the prefix of the file from which the data was derived.
        """
        if roi is not None:
            try:
                iter(roi)
            except TypeError:
                print('If an roi is give it must be an iterable!')
            self.dat_sub = self.dat[(self.dat['ti'] >= roi[0]) & (self.dat['ti'] <= roi[1])]
        else:
            self.dat_sub = self.dat

        self.grps_sub = self.dat_sub.groupby('sweep')
        colors = [self.col_dict[x] for x in vars]
        xvals = [self.horiz_dict[x] for x in vars]
        heights = [self.height_dict[x] for x in vars]

        ncols = len(np.unique(self.dat['sweep']))
        nrows = len(vars)

        fig = plt.figure(dpi=300, figsize=(0.5 * ncols, 0.75*len(vars)))
        outer = gridspec.GridSpec(1, ncols)

        for i in range(ncols):
            inner = gridspec.GridSpecFromSubplotSpec(
                nrows, 1, subplot_spec=outer[i], height_ratios=heights)
            self.plot_dat = self.grps_sub.get_group(i+1)

            for j in range(nrows):
                ax = plt.Subplot(fig, inner[j])
                ax.plot(self.plot_dat[xvals[j]], self.plot_dat[vars[j]],
                        color=colors[j], linewidth=0.5)
                self.plot_range = max(self.dat[vars[j]]) - min(self.dat[vars[j]])
                ax.set_ylim(np.min(self.dat[vars[j]]) - 0.05 * self.plot_range,
                            np.max(self.dat[vars[j]]) + 0.05 * self.plot_range)
                ax.axis('off')

                if j == 0:
                    ax.set_title(i, fontsize=8)
                else:
                    pass
                if scalebars is True and i == 0:
                    self.add_scalebars(ax, vars[j])
                else:
                    pass

                fig.add_subplot(ax)

        plt.tight_layout()
        plt.savefig(self.fullpath + '_allsweeps.pdf', dpi=300)
        plt.show()
