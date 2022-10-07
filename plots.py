data_lr_dynamic = [(256, 0.295, 11.84), (512, 0.42, 43), (1024, 0.44, 102), (2048, 0.48, 463)]
data_lr_fixed = [(256, 0.253, 92.93), (512, 0.377, 243), (1024, 0.38, 1026), (2048, 0.422, 5290)]
data_lm = [(256, 0.5652, 264), (512, 0.6404, 562), (1024, 0.6483, 1298), (2048, 0.6039, 3008.22)]
data_dclear = [(256, 0.229, 33.5), (512, 0.367, 52), (1024, 0.3359, 72), (2048, 0.43, 122)]
from matplotlib import pyplot as plt

plt.title("Cells v. RF Distance")
plt.xlabel("Num Iterations")
plt.ylabel("RF Distance")

def plot_data(arr, y_idx, labels, y_ax):
    for i, data in enumerate(arr):
        y = [float(item[y_idx]) for item in data]
        x = [float(ncells) for ncells, rf, t in data]

        plt.plot(x, y, '-o', label=labels[i])
    plt.xlabel("N Cells")
    plt.ylabel(y_ax)
    plt.title(y_ax, fontsize=20)
    plt.xticks(x)
    plt.legend()
    plt.show()

labels = ['LinRace-dynamic', 'LinRace-fixed', 'LinTIMaT', 'DCLEAR']

plot_data([data_lr_dynamic, data_lr_fixed, data_lm, data_dclear], 1, labels, 'RF Distance')
plot_data([data_lr_dynamic, data_lr_fixed, data_lm, data_dclear], 2, labels, 'Time')