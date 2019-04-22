#!/usr/bin/env python 
'''
'''
from __future__ import absolute_import, division, print_function, unicode_literals
import random
import math
import sys



def calc_ln_likelihood(data, theta):
    mu, sigma = theta
    ln_sigma = math.log(sigma)
    ln_like = -len(data)*ln_sigma
    s = 0.0
    for h in data:
        s -= (h - mu)**2
    s /= (2*sigma*sigma) 
    ln_like += s
    return ln_like

WINDOW_SIZE = 0
def main(data, theta):
    ln_likelihood = calc_ln_likelihood(data, theta)

    mu_window = WINDOW_SIZE
    sigma_window = WINDOW_SIZE
    # This is MCMC using the Metropolis algorithm:
    out = sys.stdout
    out.write("mu\tsigma\n")
    n_prop_mu = 0
    n_prop_sigma = 0
    n_accept_mu = 0
    n_accept_sigma = 0

    sum_mean  = 0.0
    sum_sd = 0.0
    num_samples = 0

    # Aim to collect 10000 samples
    if num_it < 10000:
        sample_freq = 1
    else:
        sample_freq = num_it // 10000

    # This is the Metropolis-Hastings algorithm
    for i in range(num_it):
        prev_theta = list(theta)
        if (1 + i) % sample_freq == 0:
            num_samples += 1
            sum_mean += theta[0]
            sum_sd += theta[1]
            out.write("{}\t{}\n".format(theta[0], theta[1]))
        prev_ln_likelhood = ln_likelihood
        if random.random() < 0.5:
            # change mu
            u = random.random() - 0.5
            diff = mu_window * u
            theta[0] += diff
            proposed_mu = True
            n_prop_mu += 1
        else:
            # change sigma
            u = random.random() - 0.5
            diff = sigma_window * u
            theta[1] += diff
            if theta[1] < 0.0:
                theta[1] = -theta[1]
            proposed_mu = False
            n_prop_sigma += 1
        # Prior ratio is 1.0 if we use (improper) uniform priors, so we could ignore it...
        prior_ratio = 1.0
        ln_prior_ratio = 0.0

        ln_likelihood = calc_ln_likelihood(data, theta)
        ln_likelihood_ratio = ln_likelihood - prev_ln_likelhood

        ln_posterior_ratio = ln_likelihood_ratio + ln_prior_ratio

        # Hastings ratio is 1.0, so we can ignore it...
        hastings_ratio = 1.0
        ln_hastings_ratio = 0.0

        ln_acceptance_ratio = ln_posterior_ratio + ln_hastings_ratio
        if math.log(random.random()) < ln_acceptance_ratio:
            if proposed_mu:
                n_accept_mu += 1
            else:
                n_accept_sigma += 1
        else:
            theta = prev_theta
            ln_likelihood = prev_ln_likelhood

    sys.stderr.write('mu accept %    = {:6.5f}\n'.format(n_accept_mu/n_prop_mu))
    sys.stderr.write('sigma accept % = {:6.5f}\n'.format(n_accept_sigma/n_prop_sigma))
    sys.stderr.write('Posterior mean of mu    = {}\n'.format(sum_mean/num_samples))
    sys.stderr.write('Posterior mean of sigma = {}\n'.format(sum_sd/num_samples))

if __name__ == '__main__':
    try:
        num_it = int(sys.argv[1])
        assert(num_it > 0)
        WINDOW_SIZE = float(sys.argv[2])
        theta = [float(sys.argv[3]), float(sys.argv[4])]
        datafn = sys.argv[5]
    except:
        sys.exit('''Script to estimate mean and standard deviation of the normal distribution using improper uniform priors.\n
Expecting:\npython continuous-mcmc.py <# iter> <window size> <start mu> <start sigma> <data filename>
Data should just be a column of real numbers (one per line).
''')
    try:
        data = [float(i.strip()) for i in open(datafn) if i.strip()]
        assert len(data) > 0
    except IOError:
        sys.exit('Data file "{}" does not exist.'.format(datafn))
    except ValueError:
        sys.exit('Data file "{}" was expected to have one number per line.'.format(datafn))
    except AssertionError:
        sys.exit('Data file "{}" was empty.'.format(datafn))
    main(data, theta)
