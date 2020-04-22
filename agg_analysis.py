import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector
from scipy.optimize import curve_fit


class aggFile(object):

    def __init__(self, file_path):
        self.agg_path = file_path

    def split_sweep(self, sweep):
        """
        This function finds the index of peak force and splits the sweep into
        an approach phase and a retract phase based on this index.
        """
        peak = sweep['force'].idxmax()
        approach = sweep[1:peak]
        retract = sweep[peak:np.shape(sweep)[0]]
        return(approach, retract)

    def make_sweepfile(self):
        """
        This function takes all the sweeps in the aggregate file and concatenates their approach traces into a single
        file with unique identifiers for later analysis and plotting.
        """
        agg = pd.read_csv(self.agg_path)
        appended_df = []
        for i, row in agg.iterrows():
            folder = row['construct']
            sweep_file = 'modulators/' + folder + '/' + \
                str(row['date']) + '_hek293t_' + folder + '_c' + str(row['cell']) + '_augmented.h5'
            print(sweep_file)
            t = pd.read_hdf(sweep_file)

            v = t[(t['sweep'] == row['sweep']) & (t['ti'] >= 450)
                  & (t['ti'] <= 1250)].reset_index(drop=True)

            [approach, retract] = self.split_sweep(v)
            approach = approach.assign(uniqueID=np.repeat(row['uniqueID'], np.shape(approach)[0]),
                                       construct=np.repeat(row['construct'], np.shape(approach)[0]),
                                       position_adj=(approach['position']
                                                     - approach['position'][np.max(np.where(approach['force'] <= 25))]),
                                       peaki=np.repeat(row['peaki'], np.shape(approach)[0]))

            appended_df.append(
                approach[['uniqueID', 'construct', 'position', 'absi_blsub', 'force', 'work']])

        appended_df = pd.concat(appended_df)
        output_filename = input('What would you like to call the aggregate sweep file? ')
        appended_df.to_csv(output_filename + '.csv')

    def find_slopes(self):
        """
        This function will take a summary file and allow you to fit the slopes to all the associated approach traces.
        The slopes will then be added to the summary data file.
        """
        agg = pd.read_csv(self.agg_path)

        def fitselect(xmin, xmax):
            """
            This function will fit a line to a manually selected region of the work-current plot. The fit will be
            added to the plot for visualization and the slope will be returned in pA/fJ.
            """
            indmin, indmax = np.searchsorted(self.approach['work'], (xmin, xmax))
            indmax = min(len(self.approach['work']) - 1, indmax)

            subx = self.approach['work'][indmin:indmax]
            suby = self.approach['absi_blsub'][indmin:indmax]
            try:
                poptLin = curve_fit(linFit, subx, suby)[0]
            except RuntimeError:
                print("Error - curve fit failed")

            from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
            at = AnchoredText("Sensitivity: " + str(round(poptLin[0], 3))+" (pA/fJ)",
                              prop=dict(size=8), frameon=True,
                              loc=2,)

            at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
            ax.add_artist(at)

            fitY = poptLin[0]*subx + poptLin[1]
            ax.plot(subx, fitY, '--', color='red')
            ax.axvspan(self.approach['work'][indmin],
                       self.approach['work'][indmax],
                       color='grey', alpha=0.25)
            slopes.append(round(poptLin[0], 3))
            fig.canvas.draw()

        def linFit(x, a, b):
            """
            This function defines the linear fit.
            """
            return a*x + b

        slopes = []
        for i, row in agg.iterrows():
            plt.close('all')
            folder = row['construct']
            sweep_file = 'modulators/' + folder + '/' + \
                str(row['date']) + '_hek293t_' + folder + '_c' + str(row['cell']) + '_augmented.h5'
            print(sweep_file)
            t = pd.read_hdf(sweep_file)

            v = t[(t['sweep'] == row['sweep']) & (t['ti'] >= 450)
                  & (t['ti'] <= 950)].reset_index(drop=True)

            [self.approach, self.retract] = self.split_sweep(v)

            fig, ax = plt.subplots(figsize=(15, 7))
            ax.plot(self.approach['work'], self.approach['absi_blsub'])
            ax.set_xlabel('Work(fJ)')
            ax.set_ylabel('Current (A)')

            span = SpanSelector(ax, fitselect, 'horizontal', useblit=True,
                                rectprops=dict(alpha=0.5, facecolor='red'))

            plt.draw()
            while not plt.waitforbuttonpress():
                pass
            new_column = pd.DataFrame({'slope': slopes})
        agg = agg.merge(new_column, left_index=True, right_index=True)
        agg.to_csv(self.agg_path)
        return(agg)
