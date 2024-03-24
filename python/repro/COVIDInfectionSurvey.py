#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  COVIDInfectionSurvey.py Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 02.03.2024

import numpy as np

week = np.array([0, 5, 12])
prevalence = np.array([0.5, 0.2, 0.1])

#import jax.numpy as jnp
#import jax
#import matplotlib.pyplot as plt
#trytorepro = False
### Read in data.
#if trytorepro:
#
#    ## Define cost function
#    pred = lambda week, constant, power: constant*jnp.exp(power*week)
#    mse = lambda x: jnp.sum(jnp.square(pred(week, x[0], x[1])-prevalence))
#
#    ng = 100
#
#    cgrid = np.linspace(0,1,num=ng)
#    pgrid = np.linspace(-1,0,num=ng)
#
#    err = np.zeros([ng,ng])
#    for i,c in enumerate(cgrid):
#        for j,p in enumerate(pgrid):
#            err[i,j] = mse([c,p])
#
#    fig = plt.figure()
#    plt.imshow(np.log10(err))
#    plt.savefig("err.pdf")
#    plt.close()
#
#    opt_ind = np.unravel_index(np.argmin(err), err.shape)
#    params = jnp.array([cgrid[opt_ind[0]], pgrid[opt_ind[1]]])
#
#    grad = jax.grad(mse)
#    hess = jax.jacobian(grad)
#    damp = 0.5
#
#    for i in range(100):
#        g = grad(params)
#        H = hess(params)
#        sd = jnp.linalg.solve(H, g)
#        params -= damp*sd
#        print(mse(params))
#
#    constant = float(params[0])
#    period = float(params[1])
#
#    # I am not able to recover their exact parameters; in terms of in-sample MSE the estimates obtained here are superior
#    print(mse(their_params) / mse(params))

# Ohh!! they used a linarizeation :)
lin = np.linalg.lstsq(np.stack([np.ones(3), week], axis = 1), np.log(prevalence))[0]
Constant = np.round(np.exp(lin[0]), 4)
PowerTerm = np.round(lin[1], 3)

#their_params = [0.4548,	-0.132]
#Constant = their_params[0]
#PowerTerm = their_params[1]
