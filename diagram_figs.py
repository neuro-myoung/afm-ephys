import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

axlab_dict = {'i_blsub': 'Current', 'force': 'Force',
              'work': 'Work', 'position': 'Position',
              'ti': 'Time', 'tin0': 'Time', 'tz': 'Time'}

axcol_dict = {'i_blsub': 'k', 'force': 'b', 'work': 'm', 'position': 'g'}

axhoriz_dict = {'i_blsub': 'ti', 'force': 'tin0',
                'work': 'tin0', 'position': 'tz'}

axparam_dict = {'labs': axlab_dict, 'cols': axcol_dict, 'xvals': axhoriz_dict}


def loadFile(file_path):
    """
    This function will load preprocessed data in the for of an h5 file.
    """
    dat = pd.read_hdf(file_path + '_augmented.h5')
    grps = dat.groupby('sweep')
    return(grps)


def split_sweep(sweep):
    """
    This function finds the index of peak force and splits the sweep into
    an approach phase and a retract phase based on this index.
    """
    peak = sweep['force'].idxmax()
    approach = sweep[1:peak]
    retract = sweep[peak:np.shape(sweep)[0]]
    return(approach, retract)


def reptrace_plot(filename, groups, sweep, vars, rois=None,
                  axparams=axparam_dict, peaklines=False):
    """
    This function will create a simplified representative figure of a given
    sweep identified by a grouped data frame and sweep number. The variables to
    be plotted are passed in the vars argument and will be vertically stacked.
    The data can be subsetted to zoom in on plot areas by passing a list of
    containing the start and end interval into the variable rois.If vertical
    peak annotation lines for force and current are to be displayed
    peaklines should be set to True.
    """
    plot_dat = groups.get_group(sweep)
    if rois is not None:
        plot_dat = (plot_dat[(plot_dat['ti'] >= rois[0])
                             & (plot_dat['ti'] <= rois[1])])
    else:
        pass
    fig, axs = plt.subplots(len(vars), dpi=300, figsize=(1, 0.75*len(vars)))
    if isinstance(axs, np.ndarray):
        colors = [axparams['cols'][x] for x in vars]
        xvalues = [axparams['xvals'][x] for x in vars]
        for ax, color, y, x in zip(axs, colors, vars, xvalues):
            ax.plot(plot_dat[x], plot_dat[y], color=color,
                    linewidth=0.5)

            if vars.index(y) == 0:
                axis_arrows(fig, ax, x, y, labels=axparams['labs'],
                            hide='x')
                if peaklines is True:
                    peak_lines(plot_dat, ax, ymin=-1)
                else:
                    pass
            else:
                axis_arrows(fig, ax, x, y, labels=axparams['labs'])
                if peaklines is True:
                    peak_lines(plot_dat, ax, ymin=0)
                else:
                    pass
    else:
        y = vars[0]
        color = axparams['cols'][y]
        x = axparams['xvals'][y]
        axs.plot(plot_dat[x], plot_dat[y], color=color,
                 linewidth=0.5)

        axis_arrows(fig, axs, x, y, labels=axparams['labs'])
    plt.tight_layout()
    plt.savefig(filename + '.pdf', format="pdf", dpi=300)
    plt.show()


def xy_plot(filename, groups, sweep, vars, integrate=False,
            axparams=axparam_dict):
    """
    This function will take as an argument a sweep defined by a grouped
    dataframe and sweep number and the variables to be plotted against one
    another to return a simple plot of the two variables. If integrate is set
    to true a fill will be introduced under the line.
    """
    plot_dat = groups.get_group(sweep).reset_index(drop=True)
    fig, ax = plt.subplots(nrows=1, dpi=300, figsize=(1.5, 1.5))
    [approach, retract] = split_sweep(plot_dat)
    ax.plot(approach[vars[0]], approach[vars[1]], 'k-', linewidth=0.5)
    if integrate is True:
        ax.fill_between(approach[vars[0]], 0, approach[vars[1]],
                        facecolor='m', alpha=0.6)

        ax.text(0.68, 0.15, 'Work', fontsize=8,
                transform=ax.transAxes)
    else:
        pass
    axis_arrows(fig, ax, vars[0], vars[1], axparams['labs'])
    plt.tight_layout()
    plt.savefig(filename + '.pdf', format="pdf", dpi=300)
    plt.show()


def hide_ax(fig, ax):
    """
    This function takes an axis and figure as an argument and removes the
    default axis.
    """
    for side in ['bottom', 'right', 'top', 'left']:
        ax.spines[side].set_visible(False)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')


def axis_arrows(fig, ax, x, y, labels, hide=None):
    """
    This function takes as arguments a figure, axis, and list of axis labels
    to return default axis replaced by arrows. It also has the option to hide
    either of the axis for stacked or side-by-side plots by setting hide='x' or
    hide = 'y'.
    """
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    hide_ax(fig, ax)
    dps = fig.dpi_scale_trans.inverted()
    bbox = ax.get_window_extent().transformed(dps)
    width, height = bbox.width, bbox.height

    hw = 1./20.*(ymax-ymin)
    hl = 1./20.*(xmax-xmin)
    ohg = 0

    yhw = hw/(ymax-ymin)*(xmax-xmin) * height/width
    yhl = hl/(xmax-xmin)*(ymax-ymin) * width/height

    if hide == 'y':
        ax.arrow(xmin, ymin, xmax-xmin, 0., fc='k', ec='k', lw=1,
                 head_width=hw, head_length=hl, overhang=ohg,
                 length_includes_head=True, clip_on=False)
        ax.set_xlabel(label_dict=labels[x], fontsize=8)
    elif hide == 'x':
        ax.arrow(xmin, ymin, 0., ymax-ymin, fc='k', ec='k', lw=1,
                 head_width=yhw, head_length=yhl, overhang=ohg,
                 length_includes_head=True, clip_on=False)
        ax.set_ylabel(labels[y], fontsize=8)
    else:
        ax.arrow(xmin, ymin, xmax-xmin, 0., fc='k', ec='k', lw=1,
                 head_width=hw, head_length=hl, overhang=ohg,
                 length_includes_head=True, clip_on=False)
        ax.set_xlabel(labels[x], fontsize=8)
        ax.arrow(xmin, ymin, 0., ymax-ymin, fc='k', ec='k', lw=1,
                 head_width=yhw, head_length=yhl, overhang=ohg,
                 length_includes_head=True, clip_on=False)
        ax.set_ylabel(labels[y], fontsize=8)


def peak_lines(df, ax, ymin):
    """
    This function hard codes adding vertical lines that pass through the plot
    axis at the peak force and peak current to illustrate the delay parameter.
    """
    loc_lst = [df['tin0'][df['force'].idxmax()],
               df['ti'][df['absi_blsub'].idxmax()]]
    for i in loc_lst:
        ax.axvline(x=i, ymin=ymin, ymax=1.2, c="red", linewidth=0.5, zorder=0,
                   clip_on=False)
