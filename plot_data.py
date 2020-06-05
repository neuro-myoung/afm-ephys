import matplotlib.gridspec as gridspec
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import collections as mc
import numpy as np
from math import log10, floor
from scipy.optimize import curve_fit
from scipy.optimize import minimize, fmin, leastsq
import os


class PlotData(object):
    """
    This class takes as input a file path, filename, and region of interest to make a subsetted afm-ephys data object
    that can be plotted using the defined methods.
    """

    def __init__(self, folder, filename):
        self.filename = filename
        self.folder = folder
        self.dat = pd.read_hdf(self.folder + self.filename + '_augmented.h5')
        self.dat['position'] /= 1000
        self.params = pd.read_csv(self.folder+'/params/' + self.filename + '_params.csv')
        self.grps = self.dat.groupby('sweep')
        self.summary = pd.read_csv(self.folder + self.filename + '_summary.csv')

        # Dictionaries for axis customization based on variable plotted.
        self.lab_dict = {'i_blsub': ' pA', 'force': ' nN', 'work': ' fJ', 'position': ' um', 'ti': 'ms', 'tin0': 'ms',
                         'tz': 'ms'}
        self.title_dict = {'i_blsub': 'Current', 'force': 'Force', 'work': 'Work', 'position': 'Position',
                           'ti': 'Time', 'tin0': 'Time', 'tz': 'Time'}
        self.col_dict = {'i_blsub': 'k', 'force': 'b', 'work': 'm', 'position': 'g'}
        self.horiz_dict = {'i_blsub': 'ti', 'force': 'tin0', 'work': 'tin0', 'position': 'tz'}
        self.height_dict = {'i_blsub': 2, 'force': 1, 'work': 1, 'position': 0.5}

    def plot_sweep(self, sweep, vars, roi=None, scalebars=False, scalelabs=False, checksum=False):
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

        self.fig, self.axs = plt.subplots(len(vars), dpi=300, figsize=(3, 1.5*len(vars)),
                                          gridspec_kw={'height_ratios': heights})
        #fill_dat = self.plot_dat.iloc[:np.argmax(self.plot_dat['force'])]
        for ax, x, y, color in zip(self.axs, xvals, vars, colors):
            ax.plot(self.plot_dat[x], self.plot_dat[y], color=color, linewidth=0.5)
            ax.axvline(x=self.summary['tpeaki'][sweep-1], ymin=-1.5*self.summary['peaki'][sweep-1], color='black',
                       linewidth=0.5, clip_on=False)
            ax.axvline(x=self.summary['tpeakf'][sweep-1], color='blue',
                       linewidth=0.5, clip_on=False)
            self.plot_range = max(self.dat_sub[y]) - min(self.dat_sub[y])
            self.plot_domain = ax.get_xlim()[1] - ax.get_xlim()[0]
            #ax.fill_between(fill_dat[x], y1=fill_dat[y], y2=0, color='r', linewidth=0, alpha=0.5)
            ax.set_ylim(np.min(self.dat_sub[y]) - 0.05 * self.plot_range,
                        np.max(self.dat_sub[y]) + 0.05 * self.plot_range)
            ax.axis('off')

            if scalebars is True:
                self.add_scalebars(ax, y)
            else:
                pass

            if scalelabs is True:
                self.add_scalelabs(ax, y)

            if checksum is True:
                self.check_summary_data(ax, sweep-1, y)
            else:
                pass

        plt.tight_layout()
        plt.savefig(self.folder + self.filename + '_ex-trace.pdf', dpi=300)
        plt.show()

    def plot_sweep_movie(self, sweep, vars, gif_points, roi=None, scalebars=False, scalelabs=False):
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

        frame = 0
        for i in gif_points:
            self.fig, self.axs = plt.subplots(len(vars), dpi=300, figsize=(1.5, 0.75*len(vars)),
                                              gridspec_kw={'height_ratios': heights})
            gif_subdat = self.plot_dat[self.plot_dat['ti'] <= i]
            for ax, x, y, color in zip(self.axs, xvals, vars, colors):
                ax.plot(gif_subdat[x], gif_subdat[y], color=color, linewidth=0.5)
                self.plot_range = max(self.dat_sub[y]) - min(self.dat_sub[y])
                self.plot_domain = roi[1] - roi[0]
                ax.set_ylim(np.min(self.dat_sub[y]) - 0.05 * self.plot_range,
                            np.max(self.dat_sub[y]) + 0.05 * self.plot_range)
                ax.set_xlim(roi[0], roi[1])
                ax.axis('off')

                if scalebars is True:
                    self.add_scalebars(ax, y)
                else:
                    pass

                if scalelabs is True:
                    self.add_scalelabs(ax, y)

            plt.tight_layout()
            plt.savefig(self.folder + self.filename + '_frame' + str(frame) + '.png', dpi=300)
            frame += 1

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
        lc = mc.LineCollection(lines, linewidths=1, colors='black')
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

            if max(self.plot_dat['absi_blsub']) == max(self.plot_dat[var]):
                y2 = (0.9 * np.max(self.dat_sub['i_blsub']) - ylen + (0.2*ylen))
                y1 = (0.9 * np.max(self.dat_sub['i_blsub']) - ylen - (0.5*ylen))
            else:
                y1 = (0.9 * np.min(self.dat_sub['i_blsub']) - (0.1*self.plot_range))
                y2 = (0.9 * np.min(self.dat_sub['i_blsub']) + (0.02*self.plot_range))

            ax.text(x1, y1, '100 ms', fontsize=12)
            ax.text(x2, y2, str(int(ylen)) + self.lab_dict[var], fontsize=12)
        else:
            ylen *= 2
            x1 = np.min(self.dat_sub['tin0']) + 0.02 * self.plot_domain
            y1 = 0.9 * np.max(self.dat_sub[var])-(0.75*ylen)
            ax.text(x1, y1, str(int(ylen)) + self.lab_dict[var], fontsize=12)

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
        self.dat.to_hdf(self.folder + self.filename + '_augmented.h5', key='df', mode='w')

    def plot_all_sweeps(self, vars, roi=None, scalebars=False, checksum=False):
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
                self.plot_range = max(self.dat_sub[vars[j]]) - min(self.dat_sub[vars[j]])
                ax.set_ylim(np.min(self.dat_sub[vars[j]]) - 0.05 * self.plot_range,
                            np.max(self.dat_sub[vars[j]]) + 0.05 * self.plot_range)
                ax.axis('off')

                if j == 0:
                    ax.set_title(i, fontsize=8)
                else:
                    pass
                if scalebars is True and i == 0:
                    self.add_scalebars(ax, vars[j])
                else:
                    pass

                if checksum is True:
                    self.check_summary_data(ax, i, vars[j])
                else:
                    pass

                fig.add_subplot(ax)

        plt.tight_layout()
        plt.savefig(self.folder + self.filename + '_allsweeps.pdf', dpi=300)
        plt.show()

    def check_summary_data(self, ax, sweep, val):
        """
        This function will overlay some of the calculated summary data as a check that the analysis makes sense.

        Arguments:
            ax: axis on which to plot the checks
            sweep: the sweep from which to get the summary data from
            val: the trace identity on which to plot the checks

        Returns:
            Vertical and horizontal lines showing where calculated parameters tpeakf, tpeaki, thresh, and threshind
            fall in the data.
        """
        summary = pd.read_csv(self.folder + self.filename + '_summary.csv')

        ax.axvline(x=summary.loc[sweep, 'tpeaki'], color='k', linewidth=0.25)
        ax.axvline(x=summary.loc[sweep, 'tpeakf'], color='b', linewidth=0.25)
        ax.axvline(x=self.dat_sub['ti'].reset_index(drop=True)
                   [summary.loc[sweep, 'threshind']], color='r', linewidth=0.25)
        ax.axhline(y=0, color='r', linewidth=0.5, linestyle='dashed')

        if (val == 'i_blsub') and (summary.loc[sweep, 'threshind'] != 0):
            ax.axhline(y=summary.loc[sweep, 'thresh'], color='r', linewidth=0.25)
            ax.axhline(y=-1*summary.loc[sweep, 'thresh'], color='r', linewidth=0.25)
        else:
            pass

    def doubleExpDecay(self, t,a1,a2,k1,k2,plateau):
        y = a1*np.exp(-k1*t) + a2*np.exp(-k2*t) + plateau
        return(y)

    def amplitudeError(self, primary_param_num ,params, dparams):
        summed_errs = np.abs(dparams[0]) + np.abs(dparams[1]) + np.abs(dparams[4])
        summed_vals = params[0] + params[1] + params[4]
        err = dparams[primary_param_num]/params[primary_param_num] + summed_errs/summed_vals
        return(err)

    def RelaxationFit(self, sweep, vars, roi, windowSize):
        if roi is not None:
            try:
                iter(roi)
            except TypeError:
                print('If an roi is give it must be an iterable!')

            self.dat_sub = self.dat[(self.dat['ti'] >= roi[0]) & (self.dat['ti'] <= roi[1])]
            self.plot_dat = self.dat_sub.groupby('sweep').get_group(sweep)

        else:
            self.plot_dat = self.grps.get_group(sweep)

        forceFitStart = self.summary['tpeakf'][sweep-1]
        forceFitTrace = self.plot_dat['force'][(self.plot_dat['tin0'] >= forceFitStart) &
                                               (self.plot_dat['tin0'] <= forceFitStart + windowSize)]
        forceFitx = self.plot_dat['tin0'][(self.plot_dat['tin0'] >= forceFitStart) &
                                          (self.plot_dat['tin0'] <= forceFitStart + windowSize)]-forceFitStart
        forceInit = [0.1*self.summary['peakf'][sweep-1],
                     0.8*self.summary['peakf'][sweep-1],
                     0.05,
                     0.01,
                     min(forceFitTrace)]

        iFitStart = self.summary['tpeaki'][sweep-1]
        iFitTrace = self.plot_dat['absi_blsub'][(self.plot_dat['ti'] >= iFitStart) &
                                               (self.plot_dat['ti'] <= iFitStart + windowSize)]
        iFitx = self.plot_dat['ti'][(self.plot_dat['ti'] >= iFitStart) &
                                    (self.plot_dat['ti'] <= iFitStart + windowSize)]-iFitStart
        iInit = [0.1*self.summary['peaki'][sweep-1],
                     0.8*self.summary['peaki'][sweep-1],
                     0.05,
                     0.01,
                     min(iFitTrace)]

        try:
            poptf, pcovf = curve_fit(self.doubleExpDecay, forceFitx,
                                   forceFitTrace, forceInit, maxfev=10000)
        except RuntimeError:
            print("Error - curve_fit failed")

        try:
            popti, pcovi = curve_fit(self.doubleExpDecay, iFitx,
                                   iFitTrace, iInit, maxfev=10000)
        except RuntimeError:
            print("Error - curve_fit failed")

        perrf = np.sqrt(np.diag(pcovf))
        perri = np.sqrt(np.diag(pcovi))

        colors = [self.col_dict[x] for x in vars]
        xvals = [self.horiz_dict[x] for x in vars]
        heights = [self.height_dict[x] for x in vars]

        self.fig, self.axs = plt.subplots(len(vars), dpi=300, figsize=(3, 1.5*len(vars)),
                                          gridspec_kw={'height_ratios': heights})
        for ax, x, y, color in zip(self.axs, xvals, vars, colors):
            ax.plot(self.plot_dat[x], self.plot_dat[y], color=color, linewidth=0.5)
            if y == 'force':
                ax.plot(forceFitx+forceFitStart, forceFitTrace, color='green', linewidth=0.5)
                ax.plot(forceFitx+forceFitStart, self.doubleExpDecay(forceFitx, *poptf), color='red', linewidth=0.5)
            elif y == 'i_blsub':
                ax.plot(iFitx+iFitStart, iFitTrace, color='green', linewidth=0.5)
                ax.plot(iFitx+iFitStart, self.doubleExpDecay(iFitx, *popti), color='red', linewidth=0.5)
            else:
                pass
            ax.axis('off')

        plt.tight_layout()
        plt.savefig(self.folder + self.filename + '_fit' + '.png', dpi=300)
        plt.show()

        data = dict({'force_tauf':1/poptf[3],
                     'force_taus':1/poptf[2],
                     'force_tauf_Rerr':round(perrf[3]/poptf[3]*100,1),
                     'force_taus_Rerr':round(perrf[2]/poptf[2]*100,1),
                     'force_Afast':poptf[1]/(poptf[1]+poptf[0]+poptf[4]),
                     'force_Aslow':poptf[0]/(poptf[1]+poptf[0]+poptf[4]),
                     'force_Afast_Rerr':round((self.amplitudeError(1,poptf,perrf))*100,2),
                     'force_Aslow_Rerr':round((self.amplitudeError(0,poptf,perrf))*100,2),
                     'i_tauf':1/popti[3],
                     'i_taus':1/popti[2],
                     'i_tauf_Rerr':round(perri[3]/popti[3]*100,1),
                     'i_taus_Rerr':round(perri[2]/popti[2]*100,1),
                     'i_Afast':popti[1]/(popti[1]+popti[0]+popti[4]),
                     'i_Aslow':poptf[0]/(popti[1]+popti[0]+popti[4]),
                     'i_Afast_Rerr':round((self.amplitudeError(1,popti,perri))*100,2),
                     'i_Aslow_Rerr':round((self.amplitudeError(0,popti,perri))*100,2)
                    })

        self.summary = pd.concat([self.summary, pd.DataFrame(data, index=[0])], axis=1)
        print(self.summary)
        self.summary.to_csv(self.folder + self.filename + '_summary.csv', sep=',', index=False)

def plot_all_folder(folder):
    rmstring = '_augmented.h5'
    file_list = [f.replace(rmstring, '')
                 for f in os.listdir(folder) if f.endswith('.h5')]
    print(file_list)
    for i in file_list:
        x = PlotData(folder, i)
        x.plot_all_sweeps(['position', 'force', 'work', 'i_blsub'],
                          roi=[150, 650], scalebars=True, checksum=True)
