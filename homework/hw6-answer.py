#!/usr/bin/env python
"""
Recommended prereq installation:

    virtualenv -p $(which python3) ~/env/pandas
    source  ~/env/pandas/bin/activate
    pip install pandas
    pip install scipy

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
import decimal
decimal.getcontext().prec = 50

_LOG = logging.getLogger(__file__)
if True or 'DEBUG' in os.environ:
    lvl = logging.DEBUG
else:
    lvl = logging.INFO
_LOG.setLevel(lvl)
_LOG.addHandler(logging.StreamHandler())

_LOG_ONE_HALF = log(0.5)

WORST_LN_L = float('-inf')


class ProbRecomb:
    @staticmethod
    def to_logs(r_param):
        if r_param <= 0.0 or r_param > 1.0:
            return WORST_LN_L, WORST_LN_L
        return log(r_param), log(1 - r_param)

    @staticmethod
    def recomb_prob(r_param):
        return r_param

class UsingDecimal:
    @staticmethod
    def to_logs(r_param):
        if r_param <= 0.0 or r_param > 1.0:
            return WORST_LN_L, WORST_LN_L
        dr = decimal.Decimal(r_param)
        domr = 1 - dr
        return log(r_param), float(domr.ln())

    @staticmethod
    def recomb_prob(r_param):
        return r_param

transformation = UsingDecimal
NUM_CALLS = 0
def ln_likelihood(parameter, data):
    '''Here we need to calculate the log likelihood for
    any valid point in parameters space.
    '''
    global NUM_CALLS
    NUM_CALLS += 1
    r = parameter
    log_r, log_omr = transformation.to_logs(parameter)
    ln_l = 0.0
    try:
        prev_was_recomb = False
        for event_type, num_bases in data:
            if event_type == 'R':
                ln_l += log_r + (num_bases - 1)*log_omr
            else:
                assert event_type == 'E'
                ln_l += _LOG_ONE_HALF # account for each begining at the end of the chromosome
                ln_l += (num_bases - 1)*log_omr
            #msg = 'Sum ln Pr(event={} after {} bases | r={}) HERE and add it to ln_l'
            #_LOG.debug(msg.format(event_type, num_bases, r))
    except ValueError: # we get this if we take a log 0
        return WORST_LN_L
    return ln_l


def ln_likelihood_interference(parameter, data):
    '''Make the step function lnL continuous
    '''

    r, w = parameter
    round_down_w = int(w)
    round_up_w = 1 + round_down_w
    rdw_lnl = ln_likelihood_interference_int(r, round_down_w, data)
    ruw_lnl = ln_likelihood_interference_int(r, round_up_w, data)
    frac_up = w - round_down_w
    frac_down = 1.0 - frac_up
    fake_lnl = frac_up * ruw_lnl + frac_down * rdw_lnl
    return fake_lnl

def ln_likelihood_interference_int(r_param, w, data):
    global NUM_CALLS
    NUM_CALLS += 1
    if w < 0:
        return WORST_LN_L
        return (-w)*ln_likelihood(r_param, data)
    log_r, log_omr = transformation.to_logs(r_param)
    ln_l = 0.0
    try:
        prev_was_recomb = False
        for event_type, num_bases in data:
            num_rec_opp = num_bases - 1
            curr_event_recomb = event_type == 'R'
            if (prev_was_recomb or curr_event_recomb) and num_rec_opp < w:
                return WORST_LN_L
            if prev_was_recomb:
                if curr_event_recomb:
                    if num_rec_opp > 2*w:
                        num_no_rec = num_rec_opp - 2*w
                    else:
                        num_no_rec = 0
                else:
                    num_no_rec = max(0, num_rec_opp - w)
            else:
                if curr_event_recomb:
                    num_no_rec = max(0, num_rec_opp - w)
                else:
                    num_no_rec = num_rec_opp
            if curr_event_recomb:
                ln_l += log_r + num_no_rec*log_omr
            else:
                assert event_type == 'E'
                ln_l += _LOG_ONE_HALF # account for each begining at the end of the chromosome
                ln_l += num_no_rec*log_omr
            prev_was_recomb = curr_event_recomb
                
            #msg = 'sum ln Pr(event={} after {} bases | r={}, w={} ) HERE and add it to ln_l'
            #_LOG.debug(msg.format(event_type, num_bases, r, w))
    except ValueError: # we get this if we take a log 0
        return WORST_LN_L
    return ln_l

def estimate_global_MLE(event_list, bracket):
    # Now we call the optimizer.
    # It will minimize the function passed in so we pass
    #    in a simple "adaptor" function that negates the lnL
    #    so that we maximize the lnL by minimizing -lnL
    func_adapt = lambda p : scipy_ln_likelihood(p, ln_likelihood)
    param_opt = optimize.brent(func_adapt, brack=bracket)
    return param_opt, -func_adapt(param_opt)

def estimate_global_MLE_2(event_list, params):
    func_adapt = lambda p : scipy_ln_likelihood(p, ln_likelihood_interference)
    param_opt = optimize.fmin(func_adapt, params, xtol=1e-8, disp=False)
    return param_opt, -func_adapt(param_opt)

DATA = None
def scipy_ln_likelihood(parameters, fn):
    '''This is our simple adaptor to make a minimizer maximize.'''
    negLnL = -fn(parameters, DATA)
    # _LOG.debug('ln_likelihood({p}) = {l}'.format(p=repr(parameters), l=-negLnL))
    return negLnL

def get_events_from_file(data_filepath):
    df = pd.read_csv(data_filepath, sep='\t')
    series = pd.Series([1], dtype='int64')
    data = []
    for n, row in df.iterrows():
        data.append((row[0], int(row[1])))
    return data

REALLY_BAD_LNL = -999999999999999999999999
def main(data_filepath):
    global DATA
    events = get_events_from_file(data_filepath)
    DATA = events
    min_r = 0.0
    some_guess_for_r = 0.001
    max_r = 0.5
    bracket = (min_r, some_guess_for_r, max_r)
    
    mle_1, lnL_1 = estimate_global_MLE(events, bracket)
    one_param_num_calls = NUM_CALLS
    b2, b2lnl = None, WORST_LN_L
    w = 2
    while True:
        ini = [some_guess_for_r, w]
        mle_2, lnL_2 = estimate_global_MLE_2(events, ini)
        _LOG.debug('started w={w} ln_likelihood({p}) = {l}'.format(w=w, p=repr(mle_2), l=lnL_2))
        if lnL_2 < REALLY_BAD_LNL:
            break
        if lnL_2 > b2lnl:
            b2, b2lnl = mle_2, lnL_2
        w *= 2

    mle_2, lnL_2 = b2, b2lnl
    two_param_num_calls = NUM_CALLS
    r1 = transformation.recomb_prob(mle_1)
    print('mle: r={}'.format(r1))
    print('lnL(r) = {:12}'.format(lnL_1))
    print('# calls = {:12}'.format(one_param_num_calls))

    r2 = transformation.recomb_prob(mle_2[0])
    print('mles: r={}, w={}'.format(r2, mle_2[1]))
    print('lnL(r, w) = {:12}'.format(lnL_2))
    print('# calls = {:12}'.format(two_param_num_calls))

      


if __name__ == '__main__':
    # User interface
    try:
        filename = sys.argv[1]
    except:
        sys.exit('Expecting 1 argument: the filepath to the data .xlsx file.')
    if not os.path.isfile(filename):
        sys.exit('Data file "{}" does not exist\n'.format(filename))
    main(filename)
