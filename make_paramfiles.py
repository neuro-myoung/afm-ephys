import os
import pandas as pd

keys = ['date', 'cell#',  'construct', 'protocol', 'nsweeps', 'velocity',
        'kcant', 'dkcant', 'Rp', 'Rm', 'Cm', 'Rs', 'Rscomp', 'mosm']


def make_paramfiles(folder, protocol, keys=keys):
    """
    This function creates a file containing specific experimental parameters
    that need to be input manually.

    Arguments:
        -folder: The path to a folder containing the files that need parameter
                 files made.
        -protocol: The specific protocol used in these experiments as a string
                   so that it can be removed from the filename.
        -keys: A list of keys for the parameter file. Default contains the
               following:
                date, cell#, construct, protocol, nsweeps, velocity, kcant,
                dkcant, Rp, Rm, Cm, Rs, Rscomp, mosm

    """
    rmstring = '_' + protocol + '.asc'
    file_list = [f.replace(rmstring, '') for f in os.listdir(folder) if f.endswith('.' + 'asc')]

    for i in file_list:
        param_str = ' '.join(keys)
        print(i)
        input_str = input(
            "Enter a space separated list of the following parameters: \n {}\n".format(param_str))
        input_vals = input_str.split(' ')
        if len(input_vals) != len(keys):
            input_str = input("Your input was missing values!")
        else:
            param_dict = {k: v for (k, v) in zip(keys, input_vals)}
        param_dict['uniqueID'] = '-'.join(param_dict['date'],
                                          param_dict['cell'],
                                          param_dict['protocol'])
        x = pd.DataFrame.from_dict(param_dict, orient='index', columns=['val'])
        x.to_csv(folder + i + '_params.csv')
