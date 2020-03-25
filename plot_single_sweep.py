import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

os.chdir('C:/Users/HAL/afm-ephys')

filename = 'test'
roi_start = 450
roi_end = 950

dat = pd.read_hdf(filename + '_augmented.h5')
dat_sub = dat[(dat['ti'] >= roi_start) & (dat['ti'] <= roi_end)]
grps = dat_sub.groupby('sweep')

plot_dat = grps.get_group(10)


def plot_single_sweep(sweep):
    fig, axs = plt.subplots(nrows=4, dpi=300, figsize=(2, 4),
                            gridspec_kw={'height_ratios': [0.5, 1.5, 1.5, 5]})
    colors = ['g', 'b', 'm', 'k']
    titles = ['Position', 'Force', 'Work', 'Current']
    vars = [('tz', 'position'), ('tin0', 'force'),
            ('tin0', 'work'), ('ti', 'i_blsub')]

    for ax, color, var, title in zip(axs, colors, vars, titles):
        ax.plot(plot_dat[var[0]], plot_dat[var[1]],
                color=color, linewidth=0.5)
        ax.set_title(title, size=10)
        ax.axis('off')
    plt.tight_layout()
    fig.show()
    plt.savefig(filename + '_ex-trace.pdf', dpi=300, transparent=True)


plot_single_sweep(filename)
#    positLims = findLimits(exTrace, 'position', 1.1)
#    forceLims = findLimits(exTrace, 'force', 1.1)
#    workLims = findLimits(exTrace, 'work', 1.1)
#    currentLims = findLimits(exTrace, 'i_blsub', 1.1)

#    plot_position(exTrace, ax1, 1, positLims)
#    plot_force(exTrace, ax2, 1,  forceLims)
#    plot_work(exTrace, ax3, 1, workLims)
#    plot_current(exTrace, ax4, 1, currentLims)

#    fig.subplots_adjust(wspace=0, hspace=-0.4)

# Labels for Plots
#    if labels == True:

#        axl1 = makeAxisLabel('Position', 0, 2, 0, 8)
#        axl2 = makeAxisLabel('Force', 7, 2, 0, 8)
#    axl3 = makeAxisLabel('Work', 17, 2, 0, 8)
#        axl4 = makeAxisLabel('Current', 27, 2, 0, 8)

#    else:
#        next
#    from matplotlib import collections as mc

# Scalebars for plots
#    if scalebars == True:
#        line1 = [[(455, 3000), (455, 6000)]]
#        ax1.text(457, 6000, str(3) + ' um', fontsize=6,
#         horizontalalignment='left', verticalalignment='top')

#        lc = mc.LineCollection(line1, colors='black', linewidths=0.5)
#        ax1.add_collection(lc)

#        line2 = [[(455, forceLims[1]-200), (455, forceLims[1])]]
#        ax2.text(457, forceLims[1], str(200) + ' nN', fontsize=6,
#         horizontalalignment='left', verticalalignment='top')

#    lc2 = mc.LineCollection(line2, colors='black', linewidths=0.5)
#        ax2.add_collection(lc2)

#    line3 = [[(455, workLims[1]-2e-13), (455, workLims[1])]]
#        ax3.text(457, workLims[1], str(200) + ' fJ', fontsize=6,
#                 horizontalalignment='left', verticalalignment='top')
#
#        lc3 = mc.LineCollection(line3, colors='black', linewidths=0.5)
#        ax3.add_collection(lc3)

#        line4 = [[(455, 0.9*currentLims[0]+200), (455, 0.9*currentLims[0])],
#                 [(455, 0.9*currentLims[0]), (505, 0.9*currentLims[0])]]
#        ax4.text(457, 0.9*currentLims[0]+20, str(200) + ' pA',
#                 fontsize=6, horizontalalignment='left')
#        ax4.text(457, 0.9*currentLims[0]-20, str(50) + ' ms', fontsize=6,
#             horizontalalignment='left', verticalalignment='top')

#        lc4 = mc.LineCollection(line4, colors='black', linewidths=0.5)
#        ax4.add_collection(lc4)

#    else:
#    next

#    plt.savefig(filename + '_exTrace.png', bbox_inches='tight')
