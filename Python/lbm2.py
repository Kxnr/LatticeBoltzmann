#!/usr/bin/python
##########
# Copyright (C) 2013 FlowKit Ltd, Lausanne, Switzerland
# E-mail contact: contact@flowkit.com
#
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License, either
# version 3 of the License, or (at your option) any later version.

#
# 2D flow around a cylinder
#
# Modified by Connor Keane to support arbitrary flow
#
##########

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

###### Flow definition #########################################################

SAVE_INTERVAL = 10
MAX_TIMES = 200

WIDTH = 200
HEIGHT = 200

l0 = 1.0 # meters
u0 = .1  # m/s
t0 = l0 / u0 # s
v = .0005 # kinematic viscosity in m**2 / s


dx = 1.0 / max(WIDTH, HEIGHT)  # lattice steps
dt = dx**2                     # time step
steps = 1.0 / dt               # simulation steps in 1 characteristic time
Re = u0 * l0 / v               # Reynolds number.
q = 9                          # Lattice dimensions and populations.
vlb = dt / (Re * dx**2)
uLB = dt / dx                  # Velocity in lattice units.
omega = 1.0 / (vlb / 3 + .5)  # Relaxation parameter.

###### Lattice Constants #######################################################
c = np.array([(x, y) for x in [0, -1, 1] for y in [0, -1, 1]]) # lattice velocities
t = 1/36 * np.ones(q)                                   # Lattice weights.
t[np.asarray([np.linalg.norm(ci) < 1.1 for ci in c])] = 1 / 9
t[0] = 4 / 9
noslip = [c.tolist().index((-c[i]).tolist()) for i in range(q)]
i1 = np.arange(q)[np.asarray([ci[0] < 0 for ci in c])]  # Unknown on right wall.
i2 = np.arange(q)[np.asarray([ci[0] == 0 for ci in c])]  # Vertical middle.
i3 = np.arange(q)[np.asarray([ci[0] > 0 for ci in c])]  # Unknown on left wall.

###### Function Definitions ####################################################

def curl(image):
    # naiive method is super slow!
    # maybe try to optimize with convolution
    vorticity = np.zeros((WIDTH, HEIGHT))
    for x in range(WIDTH):
        for y in range(HEIGHT):
            vorticity[x, y] = (image[1, (x+1) % WIDTH, y] -
                               image[1, (x-1) % WIDTH, y]) - \
                              (image[0, x, (y+1) % HEIGHT] -
                               image[0, x, (y-1) % HEIGHT])

    return vorticity


def gauss(x, u, s):
    return np.exp(-((x-u)**2)/s)


def sumpop(fin): return np.sum(fin, axis=0)


def equilibrium(rho, u):              # Equilibrium distribution function.
    cu = 3.0 * np.dot(c, u.transpose(1, 0, 2))
    usqr = 3/2*(u[0]**2+u[1]**2)
    feq = np.zeros((q, WIDTH, HEIGHT))
    for i in range(q):
        feq[i, :, :] = rho*t[i]*(1.+cu[i]+0.5*cu[i]**2-usqr)
    return feq


vel = np.zeros((2, WIDTH, HEIGHT))
vel = np.random.uniform(-1, 1, size=(2, WIDTH, HEIGHT))*uLB

rho = np.ones((WIDTH, HEIGHT))
#rho += fromfunction(lambda x, y: gauss(x, WIDTH/2, 10)
                  # * gauss(y, HEIGHT/2, 10), (WIDTH, HEIGHT))
#rho /= amax(rho) * 2

feq = equilibrium(rho, vel)
fin = feq.copy()
fout = feq.copy()
u = vel.copy()

###### Main time loop ##########################################################
for time in range(int(steps * MAX_TIMES)):
    '''
    if (time % SAVE_INTERVAL == 0):  # Visualization
        plt.clf()
        plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(),
                   cmap=cm.Reds, vmin=0, vmax=.175)
        # plt.imshow(rho.transpose(), cmap=cm.Reds, vmin=.2, vmax=.5)
        plt.colorbar()
        plt.savefig("turbulence."+str(time/SAVE_INTERVAL).zfill(4)+".png")
    '''
    if (time % SAVE_INTERVAL == 0):  # Visualization
        plt.clf()
        plt.imshow(curl(u).transpose(), cmap=cm.seismic)
        plt.colorbar()
        plt.savefig("curl."+str(time/SAVE_INTERVAL).zfill(4)+".png")

    for i in range(q):  # Streaming step.
        fin[i, :, :] = np.roll(np.roll(fout[i, :, :], c[i, 0], axis=0), c[i, 1], axis=1)

    rho = sumpop(fin)
    u = np.dot(c.transpose(), fin.transpose((1, 0, 2)))/rho
    feq = equilibrium(rho, u)

    fout = fin - omega * (fin - feq)  # Collision step.
