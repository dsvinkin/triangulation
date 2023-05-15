
"""

"""
import matplotlib
#matplotlib.use('Agg')

import sys
import re
import numpy as np
import matplotlib.pyplot as pl 
from matplotlib.ticker import  MultipleLocator

# шрифт
pl.rc('font',family='serif')
pl.rc('font',serif='Arial')

# расположение панелей рисунка
left, width = 0.15, 0.7
height = 0.5
bottom = 0.2
rect = [left, bottom, width, height]

# размер шрифта подписей 
label_font_size = 14

def plot_cc(
    arr_t, 
    arr_rchi2,
    arr_sig3, # two elemets T_min, Tmax
    sig3,
    t_min,
    fig_file_name,
    caption=None
    ):

    fig = pl.figure(figsize=(8.27, 11.69), edgecolor='w', facecolor='w')
    ax = fig.add_axes(rect)
       
    ax.set_ylabel("$\chi^2_r$" , fontsize=label_font_size)
    ax.set_xlabel(r'$\delta T$ (s)', fontsize=label_font_size)

    ax.plot(arr_t, arr_rchi2, marker='s', color='k', linewidth=0.5)
    ax.axhline(sig3, color='k', linewidth=0.5)
    ax.vlines(arr_sig3, [0,0], [10,10], linestyles='dashed', color='k', linewidth=0.5)
    ax.axvline(x=t_min, linestyle='dashed', color='b', linewidth=0.5)
    
    ax.set_ylim(0, 10)
    dT_err = (arr_sig3[1]-arr_sig3[0])/2.0
    ax.set_xlim(arr_sig3[0]-dT_err, arr_sig3[1]+dT_err)
    #ax.set_xlim(0, 0.003)

    ax.set_yticks(np.arange(0, 11, 1))
   
    minorLocator_x = MultipleLocator(0.1)
    minorLocator_y = MultipleLocator(0.2)
    ax.xaxis.set_minor_locator(minorLocator_x)   # x minor ticks
    ax.yaxis.set_minor_locator(minorLocator_y)  # y minor ticks

    if(caption):
        ax.set_title(caption, fontsize=18)   

    pl.savefig(fig_file_name, dpi=100)
    pl.show()