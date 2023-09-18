from collections import defaultdict
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib as mpl

mpl.rcParams['figure.dpi'] = 300


def format_data(name, mu, pd):
    data = {}
    for ncells in [256, 512, 1024, 2048]:
        with open(f'results/{name}_{ncells}_mu_{mu}_pd_{pd}.run') as f:
            rfs, times = [], []
            for line in f.readlines():
                rf, t = line.split('\t')
                rfs.append(float(rf))
                times.append(float(t))
            data[ncells] = ((np.mean(rfs), np.std(rfs)), (np.mean(times), np.std(times)))
    x = [256, 512, 1024, 2048]
    rf_y = [data[i][0] for i in x]
    times_y = [data[i][1] for i in x]

    return x, rf_y, times_y


def plot_vals(x, ax, arr, color, label):
    means = np.array([mean for mean, var in arr])
    vars = np.array([var for mean, var in arr])
    ax.plot(x, means, '-o', label=label, c=color)
    ax.fill_between(x, means - vars, means + vars, alpha=0.2, facecolor=color)


def plot_results(mu, pd):
    # Plot RF dists
    fig, ax = plt.subplots()
    datasets = ['linrace', 'dclear', 'lintimat', 'cas']
    labels = ['LinRace', 'DCLEAR-kmer', 'LinTIMaT', 'Cassiopeia-greedy']

    for i, label in enumerate(datasets):
        x, rfs, times = format_data(label, mu, pd)
        plot_vals(x, ax, rfs, clrs[i], labels[i])
    plt.legend()
    plt.xticks(x)
    plt.title(f'RF Distance', fontsize=20)
    plt.xlabel('N Cells', fontsize=15)
    plt.ylabel('RF Distance', fontsize=15)
    plt.ylim((0, 1))
    plt.show()

    # Plot Times
    fig, ax = plt.subplots()

    for i, label in enumerate(datasets):
        x, rfs, times = format_data(label, mu, pd)
        plot_vals(x, ax, times, clrs[i], labels[i])
    plt.legend()
    plt.xticks(x)
    # plt.title(f'Time', fontsize=20)
    plt.xlabel('N Cells', fontsize=20)
    plt.ylabel('Time (seconds)', fontsize=20)
    plt.show()


clrs = [(15 / 255, 182 / 255, 189 / 255), (123 / 255, 160 / 255, 22 / 255), (187 / 255, 124 / 255, 236 / 255),
        (250 / 255, 134 / 255, 114 / 255)]
mu, pd = 0.1, 1

plot_results(mu, pd)
