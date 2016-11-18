import matplotlib.pyplot as plt
import pandas as pd
import random
import os
import matplotlib as mpl
import numpy as np
import scipy.stats
mpl.rcParams['pdf.fonttype'] = 42

# Function to plot simulations
def PlotSimulations(dfs, outfile=None, estcolors=["blue","red","green"]):
    paramcolor = "black"
    prefixes = ["_x","_y",""]
    if len(dfs) == 1: prefixes = [""]
    data = reduce(lambda x, y: pd.merge(x, y, on=["chrom","start","end","mu","beta","strsd"]), dfs)
    data = data.sort(columns=["mu","beta","strsd"], ascending=False)
    order = range(data.shape[0])
    fig = plt.figure()
    fig.set_size_inches((10,4))
    numplots = 2
    # mu
    ax = fig.add_subplot(numplots, 1, 1)
    ax.plot(order, np.log10(data.mu), color=paramcolor, linewidth=1)
    for i in range(len(dfs)):
        prefix = prefixes[i]
        mvals = data["log_estmu%s"%prefix]
        ax.scatter(order, mvals, color=estcolors[i], s=5)
    ax.axhline(-8, linestyle="dashed", color="gray")
    ax.axhline(np.log10(0.05), linestyle="dashed", color="gray")
    ax.set_ylabel("Log mu")
    ax.set_xticks([])
    ax.set_xlim(left=0, right=len(order))
    ax.set_ylim(bottom=-8.5, top=-1.5)
    ax.set_yticks(range(-8,0))
    # Eff beta
    ax = fig.add_subplot(numplots, 1, 2)
    ax.plot(order, data.apply(lambda x: x.beta/x.strsd**2, 1), color=paramcolor, linewidth=1)
    for i in range(len(dfs)):
        prefix = prefixes[i]
        bvals = data.apply(lambda x: x["estbeta%s"%prefix]/((2-x["estp%s"%prefix])/x["estp%s"%prefix]**2), 1)
        ax.scatter(order, bvals, color=estcolors[i], s=5)
    ax.set_ylabel("Effective beta")
    ax.set_xticks([])
    ax.set_xlim(left=0, right=len(order))
    ax.set_ylim(bottom=-0.1, top=1)
    ax.axhline(0, linestyle="dashed", color="gray")
    ax.axhline(0.9, linestyle="dashed", color="gray")
    ax.set_yticks(np.arange(0, 1.1, 0.2))
    fig.tight_layout()
    if outfile is not None: fig.savefig(outfile)

def LoadMLData(datafile):
    data = pd.read_csv(datafile, sep="\t", names=["chrom","start","end","est_logmu_ml","est_beta_ml","est_pgeom_ml","stderr_ml","numsamples_ml","filter"])
    data["est_beta_eff_ml"] = data.apply(lambda x: x.est_beta_ml/((2-x.est_pgeom_ml)/x.est_pgeom_ml**2), 1)
    return data[["chrom","start","end","est_logmu_ml","est_beta_eff_ml","est_beta_ml", "est_pgeom_ml","stderr_ml","numsamples_ml"]]
            
def LoadML(estfile, truthfile, strsdfile, minmu=10e-9):
    ests = pd.read_csv(estfile, sep="\t", names=["chrom","start","end","log_estmu","estbeta","estp","stderr","numsamples"])
    truth = pd.read_csv(truthfile, sep="\t", names=["chrom","start","end","mu","beta"])
    strsd = pd.read_csv(strsdfile, sep="\t", names=["chrom","start","end","strsd"])
    merged = pd.merge(ests[["chrom","start","end","log_estmu","estbeta","estp","stderr"]],
                      truth, on=["chrom","start","end"])
    merged = pd.merge(merged, strsd, on=["chrom","start","end"])
    merged = merged[merged["mu"]>=minmu]
    return merged    

# Functions to plot mu and beta comparisons

def CompareMu(ydata, col1, col2, figname=None, plot=True):
    x = ydata[col1]
    y = ydata[col2]
    keep = ~np.isnan(x) & ~np.isnan(y)
    y = y[keep]
    x = x[keep]
    print "mu: %s vs. %s"%(col1, col2), scipy.stats.pearsonr(x, y), sum(keep)
    if not plot: return
    
    mutrange = [-5, -1.5]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(ydata[col1],
               ydata[col2])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left();
    ax.set_ylabel(col2, size=15)
    ax.set_xlabel(col1, size=15)
    
    ax.plot(mutrange, mutrange, color="gray", linestyle="dashed")
    ax.set_xlim(left=mutrange[0], right=mutrange[1])
    ax.set_xticklabels(ax.get_xticks(), size=12)
    ax.set_yticklabels(ax.get_yticks(), size=12);

    if figname is not None: fig.savefig(figname)
    return fig, ax

def CompareBeta(ydata, col1, col2, figname=None, minval=0):
    betarange = [minval, 1]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(ydata[col1],
               ydata[col2])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left();
    ax.set_xlabel(col1, size=15)
    ax.set_ylabel(col2, size=15)
    
    ax.plot(betarange, betarange, color="gray", linestyle="dashed")
    ax.set_ylim(bottom=betarange[0], top=betarange[1])
    ax.set_xlim(left=betarange[0], right=betarange[1])
    ax.set_xticklabels(ax.get_xticks(), size=12)
    ax.set_yticklabels(ax.get_yticks(), size=12);

    x = ydata[col1]
    y = ydata[col2]
    keep = ~np.isnan(x) & ~np.isnan(y)
    y = y[keep]
    x = x[keep]
    print "beta: %s vs. %s"%(col1, col2), scipy.stats.pearsonr(x, y), sum(keep)
    if figname is not None: fig.savefig(figname)
    return fig, ax
