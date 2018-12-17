#!/usr/bin/env python
"""
Usage:
./plot_gangstr_bamstats.py <gangstr_prefix>
"""

import matplotlib as mpl
mpl.use('Agg')
import warnings
warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
from scipy.stats import norm

try:
    gangstr_prefix = sys.argv[1]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

color = "darkblue"
normcolor = "maroon"

#### Load data ####
insdata = pd.read_csv(gangstr_prefix+".insdata.tab", delim_whitespace=True, names=["Sample","val","pdf","cdf"])
sampdata = pd.read_csv(gangstr_prefix+".samplestats.tab", delim_whitespace=True)

#### Plot GC cov ####
def GetVal(val):
    if val == -1: return None
    else: return val

gccols = [item for item in sampdata.columns if "Coverage-GCBin" in item]
binsize = 1.0/len(gccols)
xvals = [i*binsize+binsize/2 for i in range(len(gccols))]
for i in range(sampdata.shape[0]):
    sample = sampdata["Sample"].values[i]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    meancov = sampdata["MeanCoverage"].values[i]
    ax.axhline(y=meancov, color=color, linestyle="dashed")
    if len(gccols) > 0:
        gcdata = [GetVal(sampdata[col].values[i]) for col in gccols]
        ax.plot(xvals, gcdata, color=color, marker="o")
    ax.set_xlim(left=-0.05, right=1.05)
    ax.set_xlabel("GC content", size=15)
    ax.set_ylabel("Coverage", size=15)
    fig.savefig(gangstr_prefix+"."+sample+".gccov.png")
    
#### Plot ins size dist ####
for sample in set(insdata["Sample"]):
    x = insdata[(insdata["Sample"]==sample) & (insdata["pdf"]>0)]
    mean = sampdata[sampdata["Sample"]==sample]["InsMean"].values[0]
    sdev = sampdata[sampdata["Sample"]==sample]["InsSdev"].values[0]
    normval = norm(mean, sdev)
    normdist = [normval.pdf(val) for val in range(min(x["val"]), max(x["val"]+1))]
    normdist = [item/sum(normdist) for item in normdist]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x["val"], x["pdf"], color=color)
    ax.plot(range(min(x["val"]), max(x["val"]+1)), normdist, color=normcolor)
    ax.set_xlabel("Ins. size", size=15)
    ax.set_ylabel("PDF", size=15)
    fig.savefig(gangstr_prefix+"."+sample+".inspdf.png")
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x["val"], x["cdf"], color=color)
    ax.plot(range(min(x["val"]), max(x["val"]+1)), np.cumsum(normdist), color=normcolor)
    ax.set_xlabel("Ins. size", size=15)
    ax.set_ylabel("CDF", size=15)
    fig.savefig(gangstr_prefix+"."+sample+".inscdf.png")
