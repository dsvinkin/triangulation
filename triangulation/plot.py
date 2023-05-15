import sys
import os 
import numpy as np

import matplotlib as mpl
#mpl.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator


def mpl_setup():
    '''Matplotlob settings'''
    line_width = 0.8

    mpl.rcParams["figure.facecolor"] = 'w'
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True 
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True 
    mpl.rcParams['xtick.labelsize'] = 12
    mpl.rcParams['ytick.labelsize'] = 12
    mpl.rcParams['axes.facecolor'] = 'w'  
    mpl.rcParams['axes.edgecolor'] = 'k' 
    mpl.rcParams['axes.linewidth'] = line_width
    
    mpl.rcParams['axes.labelsize'] = 14
    mpl.rcParams['axes.titlesize'] = 18
    
    mpl.rcParams['xtick.major.size'] = 8
    mpl.rcParams['xtick.minor.size'] = 4
    
    mpl.rcParams['xtick.major.width'] = line_width
    mpl.rcParams['xtick.minor.width'] = line_width

    mpl.rcParams['xtick.major.pad'] = 5
    mpl.rcParams['ytick.major.pad'] = 5
    
    mpl.rcParams['figure.figsize'] = (7, 11) #A4
    #mpl.rcParams['figure.figsize'] = (3, 4) #A4
    
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
    
    mpl.rcParams['lines.linewidth'] = line_width

def get_delta_y(y_min, y_max):

    n_ticks_max = 4
    dy = int(y_max - y_min)
    mult_dy = 0.99 #1.5

    lst_num = [2, 5, 10]

    delta_y = lst_num[0]
    n_ticks = int(y_max) // delta_y

    i = 0
    while(n_ticks > n_ticks_max):
        for n in lst_num:
           delta_y = int(n * 10**i)
           y_min_curr = (int(y_min) // delta_y) * delta_y
           y_max_curr = (int(y_max) // delta_y) * delta_y + mult_dy * delta_y
           dy_curr = int(y_max_curr - y_min_curr)
           n_ticks = dy_curr // delta_y + 1
           if n_ticks <= n_ticks_max:
               break
        i = i + 1
    
    y_max_int = (int(y_max) // delta_y) * delta_y + 1.01 * delta_y
    y_max_ = (int(y_max) // delta_y) * delta_y + mult_dy * delta_y
    y_min_ = int(y_min) // delta_y * delta_y

    return delta_y, y_min_, y_max_, y_max_int

def get_delta_x(x_min, x_max):

    dx = int(1000*(x_max - x_min))
    #print(f"dx={dx}")
    #sys.exit()
    lst_num = [10, 20, 50]

    delta_ = lst_num[0]
    n_ticks = dx // delta_

    i = 0
    while(n_ticks > 10):
        for n in lst_num:
           delta_ = int(n * 10**i)
           n_ticks = dx // delta_ 
           if n_ticks <= 5:
               break
        i = i + 1

    return delta_/1000

class Plot(object):

    data_types = (
        'kw_2', 'kw_16', 'kw_64', 'kw_256',
        'gbm_0.1','gbm_1', 'gbm_2', 'gbm_16', 'gbm_64', 'gbm_256',
        'gbm_0.1_kw_2'
    )

    def __init__(self, 
        data_type, 
        caption, 
        fig_name, 
        channel_cuts=None, 
        png_dpi=100, 
        make_eps=False):
        """
        channel_cuts - list of 4 channel cuts
        """

        mpl_setup()

        self._png_dpi = png_dpi
        self._make_eps = make_eps # eps or png

        self.caption = caption
        self.fig_name = fig_name

        self.fig = plt.figure()
        self.lst_axis = []

        if not data_type in self.data_types:
            print(f"{data_type} must be one of {self.data_types}")
            return

        self.res_us = float(data_type.split('_')[-1])*1000

        if channel_cuts is not None:
           self.channel_text = [
               f'{channel_cuts[i]:.1f}-{channel_cuts[j]:.1f} keV' for i,j in [(0,3), (0,1), (1,2), (2,3)] 
           ]
        else:
           self.channel_text = "Sum G1 G2 G3".split()

        self._set_layout(self.res_us)

        #self.fig.savefig(f'test_{data_type}.png')
        #sys.exit()

    def _set_layout(self, res_us, is_double_layout=False):

        # Set panel parameters
        left, width = 0.12, 0.8  # each pannel
        heigt_channel_sum = 0.25 # heigt of the pannel with sum of channels
        heigt_gap = 0.10         # heigt of the gap between sum and channels
        heigt_channels = 0.15    # heigt of pannel with each channels
        heigt_tot = 0.88         # top point of the pannel with sum of channels
        
        # rect [left, bottom, width, height] 
 
        lst_rect = []
        lst_rect.append([left, heigt_tot - heigt_channel_sum, width, heigt_channel_sum])
        for i in range(1, 4):
            lst_rect.append([left, heigt_tot - heigt_channel_sum - heigt_gap - i*heigt_channels, width, heigt_channels])
     
        
        for i in range(len(lst_rect)):
            self.lst_axis.append(self.fig.add_axes(lst_rect[i]))
            
                
        self.lst_axis[1].set_xticklabels([])
        self.lst_axis[2].set_xticklabels([])
        self.lst_axis[0].set_xlabel(r'T-T$_{0}$ (s)')
        self.lst_axis[-1].set_xlabel(r'T-T$_{0}$ (s)')

        self.lst_axis[0].set_ylabel('counts / {:.1f} ms'.format(res_us/1000))
        self.lst_axis[2].set_ylabel('counts / {:.1f} ms'.format(res_us/1000))

        self.fig.suptitle(self.caption, y=0.98, fontsize=18)

        for i in range(len(lst_rect)):
            self.lst_axis[i].text(0.95, 0.8, self.channel_text[i], 
                fontsize=12, transform=self.lst_axis[i].transAxes, horizontalalignment='right')

        if is_double_layout:

            self.lst_axis.append(self.lst_axis[0].twinx())
            self.lst_axis.append(self.lst_axis[1].twinx())
            self.lst_axis.append(self.lst_axis[2].twinx())
            self.lst_axis.append(self.lst_axis[3].twinx())


    def _set_xscale(self, xmin, xmax, dx):

        #print(f'dx={dx}')
        x_beg = int(xmin * 1000)//int(dx * 1000) * dx
        for ax in self.lst_axis:
            ax.set_xticks(np.arange(x_beg, xmax+2*dx, dx))
            ax.set_xlim(xmin, xmax)
            ax.xaxis.set_minor_locator(MultipleLocator(dx/2))
            

    def plot_th(self, 
        arr_ti,
        arr_counts, 
        arr_bg, 
        t_min_max=[-1, 250], 
        y_limits=None, 
        t_intervals=None, 
        doShow=False
        ):
        
        """
        arr_ti - bin start times (size N)
        arr_counts - bin counts  (size N x 4)
        arr_bg - background levels (size 4)
        
        t_min_max - list of [x_min, x_max]
        y_limits - list of 4 lists containing [y_min, y_max] for each channel
        t_intervals - list of times to plot vertical lines
        """
    
        bottom = np.zeros(8)
        top = np.zeros(8)
        delta_y = np.zeros(8)

        if not y_limits is None:
            for i in range(len(self.lst_axis)):
               bottom[i], top[i] = y_limits[i]
               delta_y[i], _, _, _ = get_delta_y(bottom[i], top[i]) # should be cheked!


        loc_dic = {'Sum':4, 'G1':1, 'G2':2, 'G3':3}

        
        arr_bool = np.logical_and(arr_ti >= t_min_max[0], arr_ti <=  t_min_max[1])
        x = arr_ti[arr_bool]

        for i in range(len(self.lst_axis)):
     
            y = arr_counts[arr_bool,i]

            self.lst_axis[i].step(x, y, where='post', c='k')

            self.lst_axis[i].hlines(y=arr_bg[i], xmin=-250, xmax=250, colors='k', linestyles='dashed')

            if t_intervals is not None:
                self.lst_axis[i].axvline(x=t_intervals[0], c='c')
                self.lst_axis[i].axvline(x=t_intervals[1], c='c')

            delta_y[i], bottom[i], y_max, y_max_int = get_delta_y(np.min(y), np.max(y))
            if i == 0 or i == 1:
                top[i] = y_max_int
            else:
                top[i] = y_max
            
            #self.lst_axis[i].text(0.95, 0.8, text[i], fontsize=12, transform=self.lst_axis[i].transAxes, horizontalalignment='right')
            
            self.lst_axis[i].set_ylim(bottom[i], top[i])
            self.lst_axis[i].set_yticks(np.arange(bottom[i], top[i], delta_y[i]))
            #self.lst_axis[i].yaxis.set_minor_locator(MultipleLocator(delta_y[i]/2))
            self.lst_axis[i].set_yticks(np.arange(bottom[i], top[i], delta_y[i]/2), minor=True)

        
        #dt = int((t_min_max[1] - t_min_max[0])/5.0 * 10) /10
        dt = get_delta_x(t_min_max[0], t_min_max[1])
        self._set_xscale(t_min_max[0], t_min_max[1], dt)
    
    
        self.fig.savefig(f'{self.fig_name}.png', dpi=self._png_dpi)

        if self._make_eps:
            self.fig.savefig(f'{self.fig_name}.eps') 
           
        if doShow:
            self.fig.show()

    def plot_th_2side(self, 
        arr_ti_1,
        arr_counts_1, 
        arr_bg_1, 
        t_min_max=[-1, 250], 
        y_limits=None, 
        t_intervals=None, 
        doShow=False
        ):

        pass

def make_gbm_caption(date, time_sod, time_hhmmss):

    return f'Fermi-GBM GRB {date[2:]}\nT$_0$={time_sod} s UT ({time_hhmmss})'

def make_kw_caption(date, time_sod, time_hhmmss, det):

    return f'Konus-Wind GRB {date[2:]}\nT$_0$={time_sod} s UT ({time_hhmmss})\n{det}'

if __name__ == '__main__':

    import pandas as pd

    kw_inf = '20200415 31681.401 08:48:01.401 S1'
    gbm_inf = '20200415 31685.563746 08:48:05.563746'

    kw_data = pd.read_csv('test_data/GRB20200415_T31685/kw20200415_31681_rc2.thr', sep='\s+', header=None)
    kw_data.columns = ['Ti', 'G1', 'G2', 'G3', 'Sum']

    kw_bg = [0,0,0,0]

    gbm_data = pd.read_csv('test_data/GRB20200415_T31685/GBM/GRB200415_GBM_1ms.thr', sep='\s+', skiprows=4, header=None)
    gbm_data.columns = ['Ti', 'G1', 'G2', 'G3']
    gbm_data['Sum'] = gbm_data['G1'] + gbm_data['G2'] + gbm_data['G3']

    cuts = [21.825, 92.037, 381.327, 1559.318]

    gbm_caption = make_gbm_caption(*(gbm_inf.split()))
    kw_caption = make_kw_caption(*(kw_inf.split()))

    print(gbm_data)
    #print(kw_data[['Sum', 'G1', 'G2', 'G3']].to_numpy())
    #sys.exit()

    kw_fig = 'test_kw_2'
    kw_plot = Plot('kw_2',  kw_caption, kw_fig, cuts)
    kw_plot.plot_th(kw_data['Ti'], kw_data[['Sum', 'G1', 'G2', 'G3']].to_numpy(), kw_bg, [-0.5, 0.5])

    gbm_fig = 'test_gbm_1'
    gbm_plot = Plot('gbm_1', gbm_caption, gbm_fig, cuts)
    gbm_plot.plot_th(gbm_data['Ti'], gbm_data[['Sum', 'G1', 'G2', 'G3']].to_numpy(), kw_bg, [-0.5, 0.5])
