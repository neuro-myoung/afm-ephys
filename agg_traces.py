import pandas as pd
import numpy as np


def split_sweep(sweep):
    """
    This function finds the index of peak force and splits the sweep into
    an approach phase and a retract phase based on this index.
    """
    peak = sweep['force'].idxmax()
    approach = sweep[1:peak]
    retract = sweep[peak:np.shape(sweep)[0]]
    return(approach, retract)


filename = 'agg/agg_bychannel.csv'

x = pd.read_csv(filename)

appended_df = []
for i, row in x.iterrows():
    folder = row['construct']
    sweep_file = 'bychannel/' + folder + '/' + \
        str(row['date']) + '_hek293t_' + folder + '_c' + str(row['cell']) + '_augmented.h5'
    t = pd.read_hdf(sweep_file)

    v = t[(t['sweep'] == row['sweep']) & (t['ti'] >= 450) & (t['ti'] <= 950)].reset_index(drop=True)

    [approach, retract] = split_sweep(v)
    approach = approach.assign(uniqueID=np.repeat(row['uniqueID'], np.shape(approach)[0]),
                               construct=np.repeat(row['construct'], np.shape(approach)[0]),
                               position_adj=(approach['position']
                                             - approach['position'][np.max(np.where(approach['force'] <= 25))]),
                               peaki=np.repeat(row['peaki'], np.shape(approach)[0]))

    appended_df.append(
        approach[['uniqueID', 'construct', 'position', 'position_adj', 'force', 'work']])
appended_df = pd.concat(appended_df)
np.shape(appended_df)
appended_df.to_csv('bychannel_reptraces.csv')
