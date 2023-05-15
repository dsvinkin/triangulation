import os
import sys

import numpy as np

from ipn_project import IPN_project
from data import gbm_thr_file, kw_thr_file
from plot import Plot
from correlate import correlate, get_dTcc
from plot_cc import plot_cc

class IPN_project:

    def __init__(self, data_path, grb_name):

        self.inst_list = ['Konus', 'GBM'] + inst.ipn_instrument.inst_list
        self.dic_inst = {}

        self.data_path = os.path.join(data_path, grb_name)

        if not os.path.isdir(self.data_path):
            print(f"No data path: {self.data_path}")
            return

        self.gbm_path = os.path.join(self.data_path, 'GBM')
        self.kw_path = self.data_path
        self.ipn_instr_path = self.data_path

        self.dic_inst['Konus'] = inst.konus_wind(self.kw_path)

        if os.path.isdir(self.gbm_path):
           self.dic_inst['GBM'] = inst.fermi_gbm(self.gbm_path)

        else:
           print(f"No GBM data path: {self.gbm_path}")

        for name in inst.ipn_instrument.inst_list:
            self.dic_inst[name] = inst.ipn_instrument(name, self.ipn_instr_path)

    def get_inst(self, name):
        return self.dic_inst[name]
        

    def __str__(self):

        str_out = f'Data in {self.data_path}:\n'
        str_out +="\n".join([
            '{:10s} {:s}'.format(inst, str(self.dic_inst.get(inst, 'No data'))) for inst in self.inst_list])
        return str_out

def main():

    #data_path ='c:/work/IPN_triangulation/'
    data_path ='test_data/'
    grb_name = 'GRB20200415_T31685'

    prj = IPN_project(data_path, grb_name)
    print(prj)

    res_gmb_us = 100
    res_kw_us = 2000

    gbm = gbm_thr_file(prj.get_inst('GBM').get_lc_file(res_gmb_us))
    gbm_lcs = gbm.get_lcs()
    print (gbm_lcs['G1'])
    
 
    kw = kw_thr_file(prj.get_inst('Konus').get_lc_file(res_kw_us), prj.get_inst('Konus').info_file)
    kw_lcs = kw.get_lcs()
    print (kw_lcs['G1'])

    lst_chan = 'G1 G2 G3'.split()

    kw_bg = np.zeros(4)
    for i, ch in enumerate(lst_chan):
        kw_lcs[ch].set_bg_times(-0.1, -0.05)
        kw_bg[i+1] = kw_lcs[ch].bg_cnt
    kw_bg[0] = np.sum(kw_bg[1:])
    
    gbm_bg = np.zeros(4)
    for i, ch in enumerate(lst_chan):
        gbm_lcs[ch].set_bg_times(-0.512, -0.05)
        gbm_bg[i+1] = gbm_lcs[ch].bg_cnt
    gbm_bg[0] = np.sum(gbm_bg[1:])

    kw_fig = 'test_kw_2'
    kw_caption = kw.make_caption()
    cuts = kw.e_min + [kw.e_max[-1],]

    kw_plot = Plot('kw_2',  kw_caption, kw_fig, cuts)
    kw_plot.plot_th(kw.df['Ti'], kw.df[['Sum', 'G1', 'G2', 'G3']].to_numpy(), kw_bg, [-0.02, 0.02])

    gbm_fig = 'test_gbm_01'
    gbm_caption = gbm.make_caption()
    cuts = gbm.e_min + [gbm.e_max[-1],]

    gbm_plot = Plot('gbm_0.1', gbm_caption, gbm_fig, cuts)
    gbm_plot.plot_th(gbm.df['Ti'], gbm.df[['Sum', 'G1', 'G2', 'G3']].to_numpy(), gbm_bg, [-0.02, 0.02])


    kw_times = {'Ti': -0.002, 'Tf': 0.004}
    gbm_times = {'Ti': -0.005, 'Tf': 0.005}

    kw_ti, kw_tf =  kw_lcs['G3'].set_src_times(-0.002, 0.004)
    gbm_ti, gbm_tf =  gbm_lcs['G3'].set_src_times(-0.0045, 0.0015)

    kw_fig = 'test_kw_2_int'
    kw_plot = Plot('kw_2',  kw_caption, kw_fig, cuts)
    kw_plot.plot_th(
        kw.df['Ti'], kw.df[['Sum', 'G1', 'G2', 'G3']].to_numpy(), 
        kw_bg, [-0.02, 0.02], t_intervals=[kw_ti, kw_tf])

    gbm_fig = 'test_gbm_01_int'
    gbm_plot = Plot('gbm_0.1', gbm_caption, gbm_fig, cuts)
    gbm_plot.plot_th(
        gbm.df['Ti'], gbm.df[['Sum', 'G1', 'G2', 'G3']].to_numpy(), 
        gbm_bg, [-0.02, 0.02], t_intervals=[gbm_ti, gbm_tf])

    #sys.exit()

    scale = gbm_lcs['G3'].src_cnt / kw_lcs['G3'].src_cnt

    i_beg = 920
    i_beg_2, i_end_2 = 254, 256

    n_max_1 = 50

    
    print("gbm 01 ms i_start, t_start, n_max: {:d} {:6.3f} {:d}\n".format(
        i_beg, gbm_lcs['G3'].get_times()[i_beg], n_max_1))

    #out_x_y(lcG23.get_times(), lcG23.get_bg_sub_counts(), 'gbm_2ms.txt')

    
    arr_dt, arr_chi, nDOF, fRijMin, dTmin, iMin, nMin = correlate(
       gbm_lcs['G3'].get_times(), gbm_lcs['G3'].get_bg_sub_counts(), gbm_lcs['G3'].get_counts(),
       kw_lcs['G3'].get_times(), kw_lcs['G3'].get_bg_sub_counts(), kw_lcs['G3'].get_counts(),
       i_beg_2-1, i_end_2+1, i_beg, n_max_1, scale, 0.0001, 0.002)

    print("Cross-correlation:\ndTmin chiMin: {:.4f} {:.2f}".format(dTmin, fRijMin))

    #chi2_file_name = os.path.join(path,'ccGBM1KW2.txt')
    #out_x_y(arr_dt, arr_chi, chi2_file_name)

    dTlower, dTupper, fSigma = get_dTcc(arr_dt, arr_chi, nDOF, fRijMin, nMin, nSigma=3)
    print("dTlower dTupper fSigma: {:.4f} {:.4f} {:.3f}".format(dTlower, dTupper, fSigma))
    print("dTcc dTcc- dTcc+: {:.4f} {:+.4f} {:+.4f}".format(dTmin, dTlower-dTmin, dTupper-dTmin))

    fig_file_name = 'ccGBM{:d}KW{:d}.png'.format(res_gmb_us, res_kw_us)
    plot_cc(arr_dt, arr_chi, (dTlower, dTupper), fSigma, dTmin, fig_file_name)



main()