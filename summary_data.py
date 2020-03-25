"""
The following set of functions will use an augmented file and parameter file
in order to create a summary data file.

The function takes as an argument the file name and loads the parameter files
and augmented data file from the experiment folder. Aggregate data is
calculated for individual sweeps by splitting the data and output into a single
summary csv file. Experimental parameters from the _param.csv file are carried
over into the summary file for later analysis.

The summary file contains the following calculated values:
    -peak force (nN)
    -peak current (pA)
    -peak work (J)
    -time of peak current (ms)
    -time of peak force (ms)
    -time of peak work (ms)
    -delay between peak current and peak force (ms)
    -current offset during voltage step (pA)
    -holding voltage (mV)
    -step voltage (mV)
    -leak current (pA)
    -standard deviation of the baselie current (pA)
    -seal quality (Gohm)
    -threshold (pA)
    -threshold index
    -force threshold (nN)
    -work threshold (J)

The following experimental parameters are carried over:
    - unique experiment ID (date and cell number)
    - date
    - cell number
    - protocol
    - construct
    - series resistance (Mohm)
    - series resistance compensation (%)
    - capacitance (pF)
    - velocity (um/s)
    - cantilever stiffness (N/m)
    - uncertainty in cantilever stiffness (N/m)
    - osmolality (mOsm)
"""
import numpy as np
import pandas as pd

path = 'example/'
filename = 'test'
roi_start = 450
roi_end = 950
blsub_start = 50
blsub_end = blsub_start + 100


def window_val(df, col, window_start, window_end, func):
    """
    This function will perform an aggregate function (std or mean currently)
    on a column over a predefined time window.
    """
    window = df[col][(df['ti'] >= window_start)
                     & (df['ti'] <= window_end)]
    if func == 'mean':
        val = np.mean(window)
    elif func == 'sd':
        val = np.std(window)
    else:
        print('This function is not currently supported.')
    return(val)


def find_threshvals(pd_groups, sum_df):
    """
    This function will iterate over groups and use predefined current
    thresholds to find the index of the threshold and corresponding stimulus
    (force and work) thresholds.
    """
    val_lst = []
    fthresh_lst = []
    wthresh_lst = []
    for ind, row in sum_df.iterrows():
        df = pd_groups.get_group(ind+1)
        if row['peaki'] < (row['thresh'] + np.abs(row['offset'])):
            val = 0
            fthresh_lst.append(None)
            wthresh_lst.append(None)
        else:
            val = np.amax(np.where((df['absi_blsub'] <
                                    row['thresh'] + row['offset'])
                                   & (df['ti'] < row['tpeaki'])))
            fthresh_lst.append(df['force'].reset_index(drop=True)[val])
            wthresh_lst.append(df['work'].reset_index(drop=True)[val])
        val_lst.append(val)

    return([pd.Series(val_lst),
            pd.Series(fthresh_lst),
            pd.Series(wthresh_lst)])


def summarize(filename):
    """
    The following function will use an augmented file and parameter file in
    order to make a summary data file.

    The function takes as an argument the file name and loads the parameter
    files and augmented data file from the experiment folder. Aggregate data is
    calculated for individual sweeps by splitting the data and aggregated into
    a single csv file. Experimental parameters from the _param.csv file are
    carried over into the summary file for later analysis.

    The summary file contains the following calculated values:
        -peak force (nN)
        -peak current (pA)
        -peak work (J)
        -time of peak current (ms)
        -time of peak force (ms)
        -time of peak work (ms)
        -delay between peak current and peak force (ms)
        -current offset during voltage step (pA)
        -holding voltage (mV)
        -step voltage (mV)
        -leak current (pA)
        -standard deviation of the baselie current (pA)
        -seal quality (Gohm)
        -threshold (pA)
        -threshold index
        -force threshold (nN)
        -work threshold (J)

    The following experimental parameters are carried over:
        - unique experiment ID (date and cell number)
        - date
        - cell number
        - protocol
        - construct
        - series resistance (Mohm)
        - series resistance compensation (%)
        - capacitance (pF)
        - velocity (um/s)
        - cantilever stiffness (N/m)
        - uncertainty in cantilever stiffness (N/m)
        - osmolality (mOsm)
    """

    # Read in required files
    param = pd.read_csv(path + filename + '_params.csv').set_index('param')
    dat = pd.read_hdf(path + filename + '_augmented.h5')

    dat_sub = (dat.loc[(dat['ti'] >= roi_start) & (dat['ti'] <= roi_end)]
               .reset_index())
    dat_sub['absi_blsub'] = dat_sub['i_blsub'].transform(abs)

    # Define groups for full dataframe and subsetted dataframe.
    grps = dat.groupby('sweep')
    grps_sub = dat_sub.groupby('sweep')

    # Add columns for calculated experimental values into summary dataframe.
    agg_df = (grps.agg({'force': max,
                        'i_blsub': lambda x: abs(max(x)),
                        'work': max}).reset_index(drop=True))
    agg_df.rename(columns={"force": "peakf",
                           "i_blsub": "peaki",
                           "work": "peakw"}, inplace=True)

    agg_df = agg_df.assign(
        tpeakf=dat_sub['tin0'][(grps_sub['force']
                                .idxmax())].reset_index(drop=True),
        tpeaki=dat_sub['ti'][(grps_sub['absi_blsub']
                              .idxmax())].reset_index(drop=True),
        tpeakw=dat_sub['tin0'][(grps_sub['work']
                                .idxmax())].reset_index(drop=True),
        leak=(grps.apply(window_val, 'i', blsub_start, blsub_end, 'mean')
              .reset_index(drop=True)),
        offset=(grps.apply(window_val, 'i_blsub', roi_start-50, roi_start,
                           'mean').reset_index(drop=True)),
        stdev=(grps.apply(window_val, 'i', roi_start-100, roi_start, 'sd')
               .reset_index(drop=True)),
        vhold=(grps.apply(window_val, 'v', blsub_start, blsub_end, 'mean')
               .reset_index(drop=True)),
        vstep=(grps.apply(window_val, 'v', roi_start, roi_end, 'mean')
               .reset_index(drop=True))
    )

    agg_df = agg_df.assign(
        delay=agg_df['tpeaki'] - agg_df['tpeakf'],
        seal=agg_df['vhold']/agg_df['leak'],
        thresh=agg_df['offset'] + 5*agg_df['stdev']
    )

    col_lst = find_threshvals(grps_sub, agg_df)

    # Add columns for experimental parameters to summary file.
    agg_df = agg_df.assign(
        threshind=col_lst[0],
        fthresh=col_lst[1],
        wthresh=col_lst[2],
        date=np.repeat(param.loc['date'][0], 10),
        cell=np.repeat(param.loc['cell #'][0], 10),
        Rs=np.repeat(param.loc['Rs'][0], 10),
        Cm=np.repeat(param.loc['Cm'][0], 10),
        Rscomp=np.repeat(param.loc['Rscomp'][0], 10),
        kcant=np.repeat(param.loc['kcant'][0], 10),
        dkcant=np.repeat(param.loc['dkcant'][0], 10),
        protocol=np.repeat(param.loc['protocol'][0], 10),
        velocity=np.repeat(param.loc['velocity'][0], 10),
        construct=np.repeat(param.loc['construct'][0], 10),
        osm=np.repeat(param.loc['mosm'][0], 10),
        uniqueID=np.repeat("_".join([param.loc['date'][0],
                                     param.loc['cell #'][0]]), 10)

    )

    # Reorder columns of dataframe.
    agg_df = agg_df[['uniqueID', 'date', 'construct', 'cell', 'protocol',
                     'velocity', 'kcant', 'dkcant', 'osm', 'Rs', 'Rscomp',
                     'Cm', 'seal', 'vhold', 'vstep', 'peaki', 'tpeaki',
                     'peakf', 'tpeakf', 'peakw', 'tpeakw', 'leak', 'offset',
                     'stdev', 'seal', 'delay', 'thresh', 'threshind',
                     'fthresh', 'wthresh']]

    return(agg_df)


summary_file = summarize(filename)
summary_file.to_csv(path + filename + '_summary.csv', sep=',', index=False)
