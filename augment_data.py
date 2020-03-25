#! /usr/bin/env python
"""
Convert HEKA exported .asc file into a .csv file used for later analysis.

The following code take as input a pseudo-raw HEKA exported .asc file and
reformat it into an augmented .csv file. All values of time will be converted
into ms, current into pA, and voltage into mV. The column headers are currently
set-up for my own channel configuration but can be changed as needed.
Additional calculated parameters will be included in the new file.

The sensitivity file should share a base filename with the data appended with
the suffix _sensitivity.csv.

The parameter file should share a base filename with the data appended with the
suffix _params.csv.

The additional calculated parameters include the following:
    - piezoscanner position
    - cantilever deflection
    - baseline subtracted current
    - baseline subtracted photodetector voltage
    - force
    - relative error in the force
    - work
    - relative error in the work

The final output is a long-form .h5 file that can be subjected to later
analysis.

"""
import numpy as np
import pandas as pd
import scipy.integrate as it

filename = 'test'
nsweeps = 10
header_list = [
    'index', 'ti', 'i', 'tv', 'v',
    'tin0', 'in0', 'tz', 'z', 'tlat', 'lat'
]

blsub_start = 50
blsub_end = blsub_start + 100


def load_file(filename, nsweeps, headers):
    """
    This function will load a pseudo-raw HEKA .asc file, reformat it for later
    analysis.
    """
    dat = pd.read_csv(filename + '.asc', sep=",", header=None,
                      names=headers)

    time_cols = [col for col in dat if col.startswith('t')]
    dat[time_cols] *= 1e3
    dat['i'] *= 1e12
    dat['v'] *= 1e3

    dat['sweep'] = np.repeat(list(range(1, nsweeps + 1)), len(dat)/nsweeps)

    return(dat)


def V2nm(V):
    """
    This function will convert the piezoscanner voltage signal from a Digital
    Instruments Bioscope into distance.

    The values for the calibration are based on the most recent piezoscanner
    calibration parameters and 20x high-voltage amplifier in series with the
    waveform generator.
    """
    gain = 20
    calibration = 15.21  # nm per volt
    dist = gain * calibration * V
    return(dist)


def bl_subtraction(df, col, window_start, window_end):
    """
    This function will baseline subtract either the photodetector signal or the
    current by subtracting the mean of the designated time window.
    """
    if col == 'i':
        t = 'ti'
    elif col == 'in0':
        t = 'tin0'
    else:
        print('No time series selected!')

    base = np.mean(df[col][(df[t] >= window_start) & (df[t] <= window_end)])
    bl_sub = df[col] - base
    return(bl_sub)


def calc_work(x, y):
    """
    This function will calculate the cumulative integral of force as a function
    of piezoscanner distance traveled (work).
    """
    work = it.cumtrapz(y,  x=x,
                       initial=0) * 1e-18
    work = pd.DataFrame(work)

    return(work)


def augment_file(filename, nsweeps, window_start, window_end):
    """
    This function will use pseudo-raw HEKA data, a sensitivity calibration file
    , and some experimental meta-data to create an augmented .h5 file with
    additional calculated parameters.

    The additional calculated parameters include the following:
        - piezoscanner position
        - cantilever deflection
        - baseline subtracted current
        - baseline subtracted photodetector voltage
        - force
        - relative error in the force
        - work
        - relative error in the work

    """
    sensitivity_dat = pd.read_csv(filename + '_sensitivity.csv', sep=",",
                                  header=None)

    mean_sensitivity = np.mean(sensitivity_dat).values[0]
    std_sensitivity = np.std(sensitivity_dat).values[0]

    param_dat = pd.read_csv(filename + '_params.csv')
    kcant = float(param_dat['val'][param_dat['param'] == 'kcant'].values[0])
    dkcant = float(param_dat['val'][param_dat['param'] == 'dkcant'].values[0])

    augmented_dat = load_file(filename, headers=header_list,
                              nsweeps=nsweeps)

    grps = augmented_dat.groupby('sweep')

    # The following lines will perform the majority of the processing to
    # calculate force and work as parameters and add them to the dataframe.
    i_blsub = (grps.apply(bl_subtraction, 'i', 50, 150)
               .reset_index(drop=True))
    in0_blsub = (grps.apply(bl_subtraction, 'in0', 50, 150)
                 .reset_index(drop=True))
    deflection = in0_blsub * mean_sensitivity
    position = V2nm(augmented_dat['z']) - min(V2nm(augmented_dat['z']))
    - deflection
    force = deflection * kcant
    work = calc_work(position, force)
    rel_error = np.sqrt((std_sensitivity / mean_sensitivity) ** 2
                        + (dkcant / kcant) ** 2)

    augmented_dat = augmented_dat.assign(
        i_blsub=i_blsub,
        in0_blsub=in0_blsub,
        deflection=deflection,
        position=position,
        force=force,
        work=work,
        dforce=force * rel_error,
        dwork=work * rel_error
    )
    return(augmented_dat)


augmented_dat = augment_file('test', nsweeps, blsub_start, blsub_end)
augmented_dat.to_hdf(filename+'_augmented.h5', key='df', mode='w')
