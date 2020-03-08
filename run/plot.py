import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import collections as mc

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python plot.py MAT NUM")
        exit(0)
    mat = sys.argv[1]
    num = int(sys.argv[2])

    carrier_name = f'../data/carriers/{mat}/{num}.csv'
    measure_name = f'../data/measures/{mat}/{num}.csv'

    for fname in [carrier_name, measure_name]:
        if not os.path.exists(fname):
            print(f"File {fname} isn't available.")
            exit(1)

    carrier = pd.read_csv(carrier_name, header=None)
    lc = mc.LineCollection(carrier.values.reshape((-1, 2, 2)))
    measure = pd.read_csv(measure_name, names=['s', 'm'])
    
    fig, (ax_car, ax_mes) = plt.subplots(1, 2)
    plt.subplots_adjust(wspace=0.05)
    ax_car.add_collection(lc)
    ax_car.autoscale()
    ax_mes.plot(measure.m, measure.s)

    ax_mes.yaxis.set_visible(False)
    ax_car.set_ylabel('S')
    ax_car.set_xlabel('Energy')
    ax_mes.set_xlabel('measure')

    plt.show()
