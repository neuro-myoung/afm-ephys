import os
import pandas as pd

protocol = 'scan-80'
rmstring = '_' + protocol + '.asc'
keys = ['date', 'cell#',  'construct', 'protocol', 'nsweeps', 'velocity',
        'kcant', 'dkcant', 'Rp', 'Rm', 'Cm', 'Rs', 'Rscomp', 'mosm']


def make_paramfiles(folder, keys=keys, rmstring=rmstring):
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

        x = pd.DataFrame.from_dict(param_dict, orient='index', columns=['val'])
        x.to_csv(folder + i + '_params.csv')
