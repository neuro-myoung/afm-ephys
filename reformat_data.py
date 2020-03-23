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


#reformatted_dat = load_file(filename, headers = header_list, nsweeps = nsweeps)
#reformatted_dat.to_csv(filename + '_reformatted.csv',sep=',', index = False)
