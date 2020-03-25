import pandas as pd
import matplotlib.pyplot as plt

path = 'example/'
filename = 'test'
roi_start = 350
roi_end = 1000

dat = pd.read_hdf(path + filename + '_augmented.h5')
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
    plt.savefig(path + filename + '_ex-trace.pdf', dpi=300, transparent=True)


plot_single_sweep(filename)
