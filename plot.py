#! /bin/python

import numpy as np
import matplotlib.pyplot as plt
from myplot import *
from glob import glob

plt.rcParams.update({'font.size': 16})

def analyt(xx, yy, t):
    
    i = np.where(t >= (xx+2*yy))

    res = np.zeros_like(xx)
    res[i] = 0.5*np.sqrt(-1 + np.sqrt(1 + 16*(t - xx[i] - 2*yy[i])))

    return res



if __name__ == "__main__":

    filenames = glob("data/*.dat")

    n = len(filenames)

    for i in range(n):
        name = f"data/{i:03d}.dat"

        t = np.loadtxt(name, max_rows = 1)
        x = np.loadtxt(name, skiprows = 1, max_rows = 1)
        y = np.loadtxt(name, skiprows = 2, max_rows = 1)

        xx, yy = np.meshgrid(x, y)

        u = np.loadtxt(name, skiprows = 3).T

        fig, ax = plt.subplots(1, 2, figsize = (12, 6))

        fig.suptitle(f"Temperature distribution t = {t:.2e}")

        ax[0].imshow(u, interpolation = "bilinear", \
                extent = [x.min(), x.max(), y.min(), y.max()],\
                vmin = 0, vmax = 25,
                aspect = "auto", cmap = "Oranges_r")

        ax[0].set_xlabel ("x-axis")
        ax[0].set_ylabel ("y-axis")

        ax[1].plot(x, u[10], "-o", label = f"x = {x[10]:.2f}")
        ax[1].plot(x, u[100], "-o", label = f"x = {x[100]:.2f}")
        ax[1].plot(x, u[-10], "-o", label = f"x = {x[-10]:.2f}")

        ax[1].legend()

        ax[1].set_ylim(-1, 26)

        ax[1].set_xlabel ("x-axis")
        ax[1].set_ylabel ("Tempetature")



    save_image("123.pdf")
