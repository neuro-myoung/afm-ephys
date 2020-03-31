import numpy as np
import pandas as pd
import os


def window_val(df, col, window_start, window_end, func):
    """
    This function will perform an aggregate function (std or mean currently) on a column over a predefined time window.

    Arguments:
        df: dataframe object
        col: column name as string
        window: iterable of window boundaries
        func: aggregate function to be performed over window

    Returns:
        Float
    """
    window = df[col][(df['ti'] >= window_start) & (df['ti'] <= window_end)]
    if func == 'mean':
        val = np.mean(window)
    elif func == 'sd':
        val = np.std(window)
    else:
        print('This function is not currently supported.')
    return(val)


def find_threshvals(pd_groups, sum_df):
    """
    This function will iterate over groups and use predefined current thresholds to find the index of the threshold
    and corresponding stimulus (force and work) thresholds.

    Arguments:
        pd_groups: A pandas grouped dataframe object.
        sum_df: Partial summary dataframe.

    Returns:
        A list of pandas series corresponding to the val_lst, force thresholds, and work thresholds.
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
            val = np.amax(np.where((df['absi_blsub'] < row['thresh'] +
                                    row['offset']) & (df['ti'] < row['tpeaki'])))
            fthresh_lst.append(df['force'].reset_index(drop=True)[val])
            wthresh_lst.append(df['work'].reset_index(drop=True)[val])
        val_lst.append(val)

    return([pd.Series(val_lst),
            pd.Series(fthresh_lst),
            pd.Series(wthresh_lst)])


def summarize(fullpath, roi=None, blsub=[50, 150]):
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
    rmstring = '_scan-80'
    fullpath = fullpath.replace(rmstring, '')
    param = pd.read_csv(fullpath + '_params.csv', header=0, index_col=0)
    dat = pd.read_hdf(fullpath + '_augmented.h5')

    nsweeps = int(param.loc['nsweeps', 'val'])

    if roi is not None:
        dat_sub = (dat.loc[(dat['ti'] >= roi[0]) & (dat['ti'] <= roi[1])].reset_index())
    else:
        dat_sub = dat

    dat_sub['absi_blsub'] = dat_sub['i_blsub'].transform(abs)

    # Define groups for full dataframe and subsetted dataframe.
    grps = dat.groupby('sweep')
    grps_sub = dat_sub.groupby('sweep')

    # Add columns for calculated experimental values into summary dataframe.
    agg_df = grps_sub.agg({'force': max, 'absi_blsub': max, 'work': max}).reset_index(drop=True)
    agg_df.rename(columns={"force": "peakf", "absi_blsub": "peaki", "work": "peakw"}, inplace=True)

    agg_df = agg_df.assign(
        tpeakf=dat_sub['tin0'][(grps_sub['force'].idxmax())].reset_index(drop=True),
        tpeaki=dat_sub['ti'][(grps_sub['absi_blsub'].idxmax())].reset_index(drop=True),
        tpeakw=dat_sub['tin0'][(grps_sub['work'].idxmax())].reset_index(drop=True),
        leak=grps.apply(window_val, 'i', blsub[0], blsub[1], 'mean').reset_index(drop=True),
        offset=grps.apply(window_val, 'i_blsub', roi[0]-50, roi[0], 'mean').reset_index(drop=True),
        stdev=grps.apply(window_val, 'i', roi[0]-100, roi[0], 'sd').reset_index(drop=True),
        vhold=grps.apply(window_val, 'v', blsub[0], blsub[1], 'mean').reset_index(drop=True),
        vstep=grps.apply(window_val, 'v', roi[0], roi[1], 'mean').reset_index(drop=True),
        sweep=np.unique(dat_sub['sweep'])
    )

    agg_df = agg_df.assign(
        delay=agg_df['tpeaki'] - agg_df['tpeakf'],
        seal=agg_df['vhold']/agg_df['leak'],
        thresh=agg_df['offset'] + 10*agg_df['stdev']
    )

    col_lst = find_threshvals(grps_sub, agg_df)

    # Add columns for experimental parameters to summary file.
    agg_df = agg_df.assign(
        threshind=col_lst[0],
        fthresh=col_lst[1],
        wthresh=col_lst[2],
        date=np.repeat(param.loc['date'][0], nsweeps),
        cell=np.repeat(param.loc['cell#'][0], nsweeps),
        Rs=np.repeat(param.loc['Rs'][0], nsweeps),
        Cm=np.repeat(param.loc['Cm'][0], nsweeps),
        Rscomp=np.repeat(param.loc['Rscomp'][0], nsweeps),
        kcant=np.repeat(param.loc['kcant'][0], nsweeps),
        dkcant=np.repeat(param.loc['dkcant'][0], nsweeps),
        protocol=np.repeat(param.loc['protocol'][0], nsweeps),
        velocity=np.repeat(param.loc['velocity'][0], nsweeps),
        construct=np.repeat(param.loc['construct'][0], nsweeps),
        osm=np.repeat(param.loc['mosm'][0], nsweeps),
        uniqueID=np.repeat(param.loc['uniqueID'][0], nsweeps)

    )

    # Reorder columns of dataframe.
    agg_df = agg_df[['uniqueID', 'date', 'construct', 'cell', 'protocol', 'sweep', 'velocity', 'kcant', 'dkcant',
                     'osm', 'Rs', 'Rscomp', 'Cm', 'seal', 'vhold', 'vstep', 'peaki', 'tpeaki', 'peakf', 'tpeakf',
                     'peakw', 'tpeakw', 'leak', 'offset', 'stdev', 'seal', 'delay', 'thresh', 'threshind',  'fthresh',
                     'wthresh']]

    agg_df.to_csv(fullpath + '_summary.csv', sep=',', index=False)
    return(agg_df)


def summarize_all(folder, protocol, roi=None, blsub=[50, 150]):
    """
    Apply summarize function to all .h5 extension dataframes in a folder.

        Arguments:
            folder: Folder containing files to be analyzed.
            protocol: The name of the protocol used in the dataset.
            roi: Subset of the data containing the stimulus window (default none will not subset the sweep).
            blsub: region of the trace to be used for baseline subtraction in ms (default [50,150]).

        Returns:
            A series of files with the summary data relevant to each dataframe with same base name and suffix
            'summary'.
:        """
    rmstring = '_' + protocol + '.asc'
    file_list = [f.replace(rmstring, '') for f in os.listdir(folder) if f.endswith('.' + 'asc')]

    for i in file_list:
        print(i)
        summarize(folder + i, roi=roi, blsub=blsub)
