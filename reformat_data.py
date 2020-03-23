#! /usr/bin/env python
"""
Convert HEKA exported .asc file into a .csv file used for later analysis.

The following code will load raw HEKA exported data and reformat it into a long-
form csv file. All values of time will be converted into ms, current into pA,
and voltage into mV. The column headers are currently set-up for my own channel
configuration but can be changed as needed.
"""
import numpy as np
import pandas as pd

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
    dat = pd.read_csv(filename + '.asc', sep = ",", header=None,
        names = headers)

    time_cols = [col for col in dat if col.startswith('t')]
    dat[time_cols] = dat[time_cols].multiply(1000, axis="index")
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
    calibration = 15.21                                #nm per volt
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

    base = np.mean(df[col][(df[t] >= window_start)&(df[t] <= window_end)])
    bl_sub = df[col] - base
    return(bl_sub)


def augment_file(filename, nsweeps,
    window_start, window_end):

    sensitivity_dat = pd.read_csv(filename + '_sensitivity.csv', sep = ",",
    header=None)

    mean_sensitivity = np.mean(sensitivity_dat).values[0]
    std_sensitivity = np.std(sensitivity_dat).values[0]

    augmented_dat = load_file(filename, headers = header_list,
        nsweeps = nsweeps)

    grps = augmented_dat.groupby('sweep')

    augmented_dat = augmented_dat.assign(
        i_blsub = grps
            .apply(bl_subtraction, 'i', 50, 150)
            .reset_index(drop=True),
        in0_blsub = grps
            .apply(bl_subtraction, 'in0', 50, 150)
            .reset_index(drop=True))

    return(augmented_dat)

augmented_dat = augment_file('test', nsweeps, blsub_start, blsub_end)
print(augmented_dat.head())
