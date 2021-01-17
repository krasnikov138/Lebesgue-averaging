import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import collections as mc
import matplotlib.lines as ml

from lebesgue.bindings import Research

carriers_path = '../data/carriers'

def find_carrier(mat, lower, upper):                                                                                                                                                                                                                  
    dirname = os.path.join(carriers_path, mat)                                             
    for fname in os.listdir(dirname):                                    
        fullname = os.path.join(dirname, fname)                                       
        df = pd.read_csv(fullname, header=None, names=['e1', 'cs1', 'e2', 'cs2'])   
        if ((df.e1 >= lower) & (df.e1 <= upper)).any():
            print(fullname)


def plot_carriers(config, scale=None):                                                                                                                                                                                                                
    fig, ax = plt.subplots()
    legend_handles = []
    for name, info in config.items():
        fname = os.path.join(carriers_path, name, f"{info['i']}.csv")                          
        carrier = pd.read_csv(fname, header=None)
        lc = mc.LineCollection(carrier.values.reshape((-1, 2, 2)), colors=info['c'])
        ax.add_collection(lc)
        legend_handles.append(ml.Line2D([], [], color=info['c'], label=name))
    plt.legend(handles=legend_handles, fontsize=15)                 
    ax.autoscale()                                            
    if scale:                                                              
        ax.set_yscale(scale)
    ax.set_xlabel('energy (eV)', fontsize=15)
    ax.set_ylabel('cross section (barn)', fontsize=15)
    ax.set_title('Carrier resonance partition', fontsize=15)      
    plt.show()


config = {
    'U238': {
        'cross_section_file': 'tape43', 
        'concentration': 4.5}, 
    'U235': {
        'cross_section_file': 'tape41', 
        'concentration': 22}, 
    'Fe56': {
        'cross_section_file': 'tape47', 
        'concentration': 1.5}, 
    'C12': {
        'cross_section_file': 'tape45', 
        'concentration': 14.8}
}
data_path = '../../NJOY21/bin'
for item in config.values():
    item['cross_section_file'] = os.path.join(data_path, item['cross_section_file'])
r = Research(config, 0.1, 1)

def plot_measure(mat, i, cs):
    fig, ax = plt.subplots()
    fname = os.path.join(carriers_path, mat, f'{i}.csv')
    carrier = pd.read_csv(fname, header=None)
    boundary_min = carrier.loc[:, 0].min()
    boundary_max = carrier.loc[:, 2].max()
    lc = mc.LineCollection(carrier.values.reshape((-1, 2, 2)), colors='b')
    ax.add_collection(lc)
    for p in cs:
        params = r.measures_const_section_params(mat, i, p)[0]
        line = [[boundary_min, boundary_max], [p, p]]
        ax.plot(line[0], line[1], 'r--')
        ax.scatter(params, [p] * len(params), c='r', s=15)
        ax.text(boundary_max, p, f'CS = {p:.2f}', fontsize=15)
    ax.set_xlabel('energy (eV)', fontsize=15)
    ax.set_ylabel('cross section (barn)', fontsize=15)
    ax.autoscale()
    ax.set_xlim([
        boundary_min - (boundary_max - boundary_min) * 0.02, 
        boundary_max + (boundary_max - boundary_min) * 0.08
    ])
    ax.set_title('Optimal measure grid', fontsize=15)
    plt.show()


plot_carriers({'U235': {'i': 60, 'c': 'b'}, 'U238': {'i': 19, 'c': 'r'}}) 
grid = r.measures_optimal_grid_rational('U235', 65, 5, 0.01)
plot_measure('U235', 65, grid[1:])
