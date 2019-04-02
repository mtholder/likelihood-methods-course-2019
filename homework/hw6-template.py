#!/usr/bin/env python
"""
Recommended prereq installation:

    virtualenv -p $(which python3) ~/env/pandas
    source  ~/env/pandas/bin/activate
    pip install pandas
    pip install scipy
    pip install xlrd

To use:

    source  ~/env/pandas/bin/activate
    python hw6-template.py Mutt_Gamete_data.xlsx


"""

from scipy import optimize
import pandas as pd
from math import log, exp
import logging
import random
import sys
import os


_LOG = logging.getLogger(__file__)
if True or 'DEBUG' in os.environ:
    lvl = logging.DEBUG
else:
    lvl = logging.INFO
_LOG.setLevel(lvl)
_LOG.addHandler(logging.StreamHandler())


WORST_LN_L = float('-inf')
def ln_likelihood(parameter, data):
    '''Here we need to calculate the log likelihood for
    any valid point in parameters space.
    '''
    r = parameter
    ln_l = 0.0
    try:
        for event_type, num_bases in data:
            msg = 'Calculate ln Pr(event={} after {} bases | r={}) HERE and add it to ln_l'
            _LOG.debug(msg.format(event_type, num_bases, r))
    except ValueError: # we get this if we take a log 0
        return WORST_LN_L
    return ln_l

def estimate_global_MLE(event_list):
    min_r = 0.0
    some_guess_for_r = 0.001
    max_r = 0.5
    bracket = (min_r, some_guess_for_r, max_r)
    # Now we call the optimizer.
    # It will minimize the function passed in so we pass
    #    in a simple "adaptor" function that negates the lnL
    #    so that we maximize the lnL by minimizing -lnL
    param_opt = optimize.brent(scipy_ln_likelihood, brack=bracket)
    return param_opt, -scipy_ln_likelihood(param_opt)

DATA = None
def scipy_ln_likelihood(parameters):
    '''This is our simple adaptor to make a minimizer maximize.'''
    negLnL = -ln_likelihood(parameters, DATA)
    _LOG.debug('ln_likelihood({p:12.6f}) = {l}'.format(p=parameters, l=-negLnL))
    return negLnL

def get_events_from_file(data_filepath):
    df = pd.read_csv(data_filepath, sep='\t')
    series = pd.Series([1], dtype='int64')
    data = []
    for n, row in df.iterrows():
        data.append((row[0], int(row[1])))
    return data

# DON'T CALL THIS. This is the code I used to clean the data
# def mth_data_cleaning_get_chromosomes_from_file(data_filepath):

#     df = pd.read_excel(data_filepath)
#     series = pd.Series([2, 3], dtype='int64')
#     data = []
#     current_chromosome = None
#     for n, row in df.iterrows():
#         chromosome_num = row[0]
#         if not pd.isna(chromosome_num):
#             if current_chromosome:
#                 data.append(current_chromosome)
#             current_chromosome = []
#         current_chromosome.append(int(row[3]))
#     if current_chromosome:
#         data.append(current_chromosome)
    
#     with open('mutt_gamete_event.tsv', 'w') as outp:
#         header = ['event', 'num_base_pairs_to_prev']
#         outp.write('{}\n'.format('\t'.join(header)))
#         for c in data:
#             num_events = len(c)
#             prev_pos = 0
#             for index, position in enumerate(c):
#                 distance = str(position - prev_pos)
#                 if index + 1 == num_events:
#                     data_row = ['E', distance]
#                 else:
#                     data_row = ['R', distance]
#                 prev_pos = position
#                 outp.write('{}\n'.format('\t'.join(data_row)))
#     return data

def main(data_filepath):
    global DATA
    events = get_events_from_file(data_filepath)
    DATA = events
    mle, lnL = estimate_global_MLE(events)
    print('mle = {:12.6f}'.format(mle))
    print('lnL = {:12.6f}'.format(lnL))

      


if __name__ == '__main__':
    # User interface
    try:
        filename = sys.argv[1]
    except:
        sys.exit('Expecting 1 argument: the filepath to the data .xlsx file.')
    if not os.path.isfile(filename):
        sys.exit('Data file "{}" does not exist\n'.format(filename))
    main(filename)
