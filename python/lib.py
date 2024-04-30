
import numpy as np

# Settings
nsamp = 10000
regtype = 'beta'

def get_mean_post_prop(success, total, nsamp = 100):
    prop_succ = success / total
    post_alpha = success + 0.5
    post_beta = total-success + 0.5
    post_mean = post_alpha / (post_alpha + post_beta)
    post_var = post_alpha*post_beta / (np.square(post_alpha+post_beta)*(post_alpha+post_beta+1))
    post_sd = np.sqrt(post_var)
    samp = np.random.beta(post_alpha, post_beta, size=nsamp)
    return {'est' : prop_succ, 'samp' : samp}
