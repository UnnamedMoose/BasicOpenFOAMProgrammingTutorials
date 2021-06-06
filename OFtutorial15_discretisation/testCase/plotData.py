#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  6 11:04:35 2021

@author: artur
"""

import pandas
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)

# %% Read the data.
refData = {f.split("_")[1].replace(".csv", ""): pandas.read_csv(os.path.join("refData", f))
           for f in os.listdir("./refData")}

# %% Plot solutions obtained using each scheme.
fig, ax = plt.subplots(1)
ax.set_xlabel("x [m]")
ax.set_ylabel(r"$\phi$ [-]")
for s in refData:
    ax.plot(refData[s]["x"], refData[s]["T"], "o--", lw=1, ms=4, mew=2, label=s)
xlim = ax.get_xlim()
ax.hlines([0, 1], xlim[0], xlim[1], color="k", ls="--", alpha=0.5, zorder=-10)
ax.set_xlim(xlim)
ax.legend(loc="lower left", framealpha=1)
plt.tight_layout()
