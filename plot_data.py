import pandas as pd
import matplotlib.pyplot as plt


class PlotData(object):
    """
    This creates an afm-ephys data object that can be plotted using the defined
    methods.
    """

    def __init__(self, path, filename, roi, scalebars=False,
                 scalelabs=False, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filename = filename
        self.fullpath = path + filename
        self.roi = roi

        self.dat = pd.read_hdf(self.fullpath + '_augmented.h5')
        self.dat_sub = self.dat[((self.dat['ti'] >= self.roi[0])
                                 & (self.dat['ti'] <= self.roi[1]))]

        self.grps = self.dat_sub.groupby('sweep')

    def plot_sweep(self, sweep):
        """
        This function will take as an argument a single sweep of pre-processed
        afm-ephys data and create a representative trace.

        The traces shown by default are position, force, work, and current. The
        traces are stacked vertically and aligned on the time axis. By default
        there are no axis or labels save for the titles indicating the data
        shown in each plot.
        """
        self.plot_dat = self.grps.get_group(sweep)
        self.fig, self.axs = plt.subplots(nrows=4, dpi=300, figsize=(2, 4),
                                          gridspec_kw={'height_ratios': [0.5, 1.5, 1.5, 5]})
        colors = ['g', 'b', 'm', 'k']
        titles = ['Position', 'Force', 'Work', 'Current']
        vars = [('tz', 'position'), ('tin0', 'force'),
                ('tin0', 'work'), ('ti', 'i_blsub')]

        for ax, color, var, title in zip(self.axs, colors, vars, titles):
            ax.plot(self.plot_dat[var[0]], self.plot_dat[var[1]],
                    color=color, linewidth=0.5)
            ax.set_title(title, size=10)
            ax.axis('off')
        plt.tight_layout()
        plt.show()
        plt.savefig(self.fullpath + '_ex-trace.pdf', dpi=300,
                    transparent=True)

    def add_scalebars(self):
        pass

    def add_scalelabs(self):
        pass
