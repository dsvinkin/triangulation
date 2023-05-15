
"""

"""
import sys
import re
import numpy as np  

import pandas as pd

import clock

class light_curve():

    def __init__(self, sc_name, date, time, dT, E_min, E_max, bg_cnt, lc_data):

        self.sc_name = sc_name
        self.date = date 
        self.time = time 
        self.dT = dT
        self.E_min = E_min
        self.E_max = E_max 
        self.bg_cnt = bg_cnt
        self.lc_data = lc_data # array of bin start times and counts

        self.bg_idx_i = 0
        self.bg_idx_f = 0
        self.src_idx_i = 0
        self.src_idx_f = 0

    def __str__(self,):
        return f"{self.sc_name} {self.date} {self.time} {self.dT}"
        
    def find_index(self, t):
        return np.argwhere(self.lc_data[:,0] <= t).flatten()[-1]

    def get_bg_rate(self):
        if self.bg_cnt:
            return self.bg_cnt/(self.lc_data[1,0] - self.lc_data[0,0])
        else:
            return None

    def get_bg_counts(self):
        return np.ones_like(self.lc_data[:,0]) * self.bg_cnt

    def get_bg_sub_counts(self):
        return self.lc_data[:,1] - self.bg_cnt
        
    def get_times(self):
        return self.lc_data[:,0]
    
    def get_counts(self):
        return self.lc_data[:,1]

    def get_counts_err(self):
        return (self.lc_data[:,1])**0.5

    def set_bg_rate(self, bg_rate):
        res = self.lc_data[1,0] - self.lc_data[0,0]
        self.bg_cnt = bg_rate * res
    
    def set_bg_times(self, Ti, Tf):
        b = np.logical_and(self.lc_data[:,0] >= Ti, self.lc_data[:,0] <= Tf)
        self.bg_cnt = np.sum(self.lc_data[b,1])/self.lc_data[b,1].size
 
        w = np.argwhere(b).flatten()
        #print w
        self.bg_idx_i = w[0]
        self.bg_idx_f = w[-1]
        print ("{:s} background interval: {:8d} {:8d}".format(self.sc_name, self.bg_idx_i, self.bg_idx_f))
        return self.lc_data[self.bg_idx_i, 0], self.lc_data[self.bg_idx_f, 0]

    def set_bg_idx(self, iBeg, iEnd):

        if iBeg < 0 or iEnd < 0 or iEnd >= self.lc_data.shape[0] or iBeg >= self.lc_data.shape[0]:
            print(f"Wrong interval boundaries: {iBeg} {iEnd}")
            print(f"Allowed indexes are between: 0 and {self.lc_data.shape[0]}")
            return

        if iEnd <= iBeg:
            print(f"Wrong interval boundaries: {iBeg} {iEnd}")
            print("The first index should be lower than the second")
            return

        self.bg_cnt = np.sum(self.lc_data[iBeg:iEnd+1,1])/self.lc_data[iBeg:iEnd+1,1].size
 
        self.bg_idx_i = iBeg
        self.bg_idx_f = iEnd
    
        return self.lc_data[self.bg_idx_i, 0], self.lc_data[self.bg_idx_f, 0]

    def set_src_times(self, Ti, Tf):

        b = np.logical_and(self.lc_data[:,0] >= Ti, self.lc_data[:,0] <= Tf)

        w = np.argwhere(b).flatten()
        self.src_idx_i = w[0]
        self.src_idx_f = w[-1] - 1

        print ("{:s} source interval: {:8d} {:8d}".format(
            self.sc_name, self.src_idx_i, self.src_idx_f))

        if(self.bg_cnt):
            self.src_cnt = np.sum(self.lc_data[self.src_idx_i:self.src_idx_f+1,1] - self.bg_cnt)
        else:
            print ("{:s} background is not set!".format(self.sc_name))
            return None, None

        return self.lc_data[self.src_idx_i, 0], self.lc_data[self.src_idx_f + 1, 0]

    def set_src_idx(self, iBeg, iEnd):
        if(self.bg_cnt):
            self.src_cnt = np.sum(self.lc_data[iBeg:iEnd+1,1] - self.bg_cnt)
        else:
            print ("{:s} background is not set!".format(self.sc_name))
            return None, None
 
        self.src_idx_i = iBeg
        self.src_idx_f = iEnd
        print ("Src interval: {:8d} {:8d}".format(self.src_idx_i, self.src_idx_f))
        return self.lc_data[self.src_idx_i, 0], self.lc_data[self.src_idx_f + 1, 0]
        
    def write_log(self):
        file_name  = "{:s}_{:05d}_{:s}_par.txt".format(self.date, int(self.time), self.sc_name)
        with open(file_name, 'w') as f:
            f.write("bg:  {:5d} {:5d}\n".format(self.bg_idx_i,self.bg_idx_f))
            f.write("src: {:5d} {:5d}\n".format(self.src_idx_i,self.src_idx_f))

def text_to_array(text, skip_lines):

    lst_data = []
    for l in text.split('\n')[skip_lines:]:
        if len(l.split()) == 0:
            continue
        lst_data.append([float(s) for s in l.split()])

    return np.array(lst_data)

class ipn_file():
    """
    Class for reading IPN formatted lightcurve file
    """    
    def __init__(self, file_name):

        self.file_name = file_name
        self.lc, self.str_info = self.read(file_name)
        
    def read(self, file_name):
        """Read IPN thr and returns lc object"""

        with open(file_name, 'r') as f:
            text = f.read()
 
        s = " ".join(text.split('\n')[0:3])
        
        flt = r'[-+]?\d*\.?\d*(?:[eE][-+]?\d+)?' # float regexp
        h1 = r"\'\s*(.+?)\s*\'\s+\'(\d{2})/(\d{2})/(\d{2})\'\s+" + r"({:s})\s+".format(flt)
        h2 = r"({:s})\s+({:s})\s*".format(flt, flt)
        h = h1 + h2 + h2

        m = re.search(h, s)
        if m:
            sc_name = m.group(1)
            date = m.group(4) + m.group(3) + m.group(2)
            time = float(m.group(5))
            E1 = float(m.group(6))
            E2 = float(m.group(7))
            bg = float(m.group(8))
            dt = float(m.group(9))
        else:
            #print h, s
            return None
    
        if int(date[:2]) < 90:
            date = '20' + date # 21 century
        else:
            date = '19' + date # 20 century

        str_info = "Header info: {:10s}\nDate time: {:6s} {:9.3f}\n".format(sc_name, date, time)
        str_info +="Emin Emax: {:8.1f} {:8.1f}\nBg lev.:{:5.1f}\nres: {:5.3f}\n".format(E1, E2, bg, dt)
        
        lst_data = []
        for l in text.split('\n')[3:]:
            if len(l.split()) == 0:
                continue
            lst_data.append([float(s) for s in l.split()])

        data = np.array(lst_data)
        
        lc = light_curve(sc_name, date, time, dt, E1, E2, bg, data)

        str_info +="t_min t_max: {:8.3f} {:8.3f}\nlc_zize: {:d}\n".format(
            lc.get_times()[0], lc.get_times()[-1], lc.get_times().size)
            
        return lc, str_info

    def get_lc(self):
        return self.lc

    def get_info(self):
        return self.str_info

class gbm_thr_file():

    def __init__(self, file_name):
        self.file_name = file_name
       
        with open(file_name) as f:
           self.text = f.read()

        self.str_info = self.read(self.text)

    def split_header(self, text):

        lines = text.split('\n')
        m = re.search(r'(\d+\.\d+)', lines[1])
        if not m:
            print("Cannot split header of {:s}".format(self.file_name))
            sys.exit()

        #print(m.group(0), lines[1])
        self.trig_time = float(m.group(1))

        date_utc = clock.fermi2utc(self.trig_time)
        frac_s = date_utc.microsecond/1e6
        time2sec = date_utc.hour*3600 + date_utc.minute*60 + date_utc.second
        self.time_sod = time2sec + frac_s

        self.date = date_utc.strftime("%Y%m%d")
        self.time_hhmmss = date_utc.strftime("%H%M%S.%f")
        #print(self.date, self.sod)

        self.e_min = [float(s) for s in lines[2].split()[-3:]]
        self.e_max = [float(s) for s in lines[3].split()[-3:]]

    def make_caption(self,):

        return f'Fermi-GBM GRB{self.date[2:]}\nT$_0$={self.time_sod} s UT ({self.time_hhmmss})'

    def read(self, text):
        """
        reads GBM thr and returns lc object
        """

        self.split_header(text)
        
        sc_name = 'Fermi-GBM'
        
        bg = 0.0
        
        data = text_to_array(text, 4)
        dt = np.around(data[1,0] - data[0,0], decimals=4)

        self.lc_names = "G1 G2 G3 G23 G123".split()

        self.df = pd.DataFrame(data=data, columns = ['Ti', 'G1', 'G2', 'G3'])
        self.df['Sum'] = self.df['G1'] + self.df['G2'] + self.df['G3'] 

        self.lc = {
           'G1' :light_curve(sc_name, self.date, self.time_sod, dt, 
                self.e_min[0], self.e_max[0], bg, data[:,0:2]),

           'G2': light_curve(sc_name, self.date, self.time_sod, dt, 
                self.e_min[1], self.e_max[1], bg, np.stack((data[:,0], data[:,2]), axis=-1)),

           'G3': light_curve(sc_name, self.date, self.time_sod, dt, 
                self.e_min[2], self.e_max[2], bg, np.stack((data[:,0], data[:,3]), axis=-1)),

           'G23': light_curve(sc_name, self.date, self.time_sod, dt, 
                self.e_min[1], self.e_max[2], bg, 
                np.stack((data[:,0], data[:,2]+data[:,3]), axis=-1)),

           'G123': light_curve(sc_name, self.date, self.time_sod, dt, 
                self.e_min[0], self.e_max[2], bg, 
                np.stack((data[:,0], data[:,1]+data[:,2]+data[:,3]), axis=-1))
        }

        str_info = f"Header info: {sc_name:10s}\nDate time: {self.date:8s} {self.time_sod:9.3f}\n"
        str_info += f"Emin Emax: {self.e_min[0]:5.1f} {self.e_max[2]:8.1f}\n"
        str_info += f"Bg lev.:{bg:5.1f}\nres: {dt:5.3f}\n"
        str_info += "t_min t_max: {:8.3f} {:8.3f}\nlc_zize: {:d}\n".format(
            self.lc['G1'].get_times()[0], self.lc['G1'].get_times()[-1], 
            self.lc['G1'].get_times().size)

        return str_info
     
    def get_lcs(self):
        return dict([ (s, self.lc[s]) for s in self.lc_names ])

    def get_info(self):
        return self.str_info


class kw_thr_file():

    def __init__(self, thr_name, info_name):
        """
        file_name is one of
        kwYYYYMMDD_SSSSS_rc[2,16,64,256].thr
        

        assume that *.info file is in the same path
        """
   
        self.thr_name = thr_name
        self.info_name = info_name

        self.e_min = [float(s) for s in "20 90 38".split()]
        self.e_max = [float(s) for s in "90 38 1600".split()]

        self.info = self.parse_info(info_name)
        self.date, self.time_sod, self.time_hhmmss, self.det = \
            self.info['date'], float(self.info['time']), self.info['time_hhmmss'], self.info['det']

        self.str_info = self.read(thr_name)

    def make_caption(self,):

        return f'Konus-Wind GRB{self.date[2:]}\nT$_0$={self.time_sod} s UT ({self.time_hhmmss})\n{self.det}'

    def parse_info(self, info_file):

        with open(info_file) as f:
            lines = f.read().split('\n')

        dic_info = {}
        for l in lines:
            lst_ = l.split('=')
            if len(lst_) == 2:
                dic_info[lst_[0]] = lst_[1]

        lst_details = dic_info['lzdetails'].split()[2:]

        short_time = '{:05d}'.format(int(float(lst_details[1])))

        return {
            'date':        lst_details[0], 
            'short_date':  lst_details[0][2:],
            'time_hhmmss': lst_details[2].lstrip('(').rstrip(')'),
            'short_time':  short_time,
            'time':        "{:s}".format(lst_details[1]),
            'det':         "{:s}".format(lst_details[3]),
            'name_long' :   "GRB{:s}_T{:s}".format(lst_details[0], short_time)
        }

    def read(self, file_name):
        """
        reads KW thr and returns lc object
        """
        
        sc_name = 'Konus-Wind'
        
        bg = 0.0
        
        with open(file_name, 'r') as f:
           text = f.read()

        data = text_to_array(text, 1)

        dt = np.around(data[1,0] - data[0,0], decimals=3)

        self.lc_names = "G1 G2 G3 G23 G123".split()

        self.df = pd.DataFrame(data=data, columns=['Ti', 'G1', 'G2', 'G3', 'Sum'])

        self.lc = {
           'G1' :light_curve(sc_name, self.date, self.time_sod, dt, 
                self.e_min[0], self.e_max[0], bg, data[:,0:2]),

           'G2': light_curve(sc_name, self.date, self.time_sod, dt, 
                self.e_min[1], self.e_max[1], bg, np.stack((data[:,0], data[:,2]), axis=-1)),

           'G3': light_curve(sc_name, self.date, self.time_sod, dt, 
                self.e_min[2], self.e_max[2], bg, np.stack((data[:,0], data[:,3]), axis=-1)),

           'G23': light_curve(sc_name, self.date,self.time_sod, dt, 
                self.e_min[1], self.e_max[2], bg, 
                np.stack((data[:,0], data[:,2]+data[:,3]), axis=-1)),

           'G123': light_curve(sc_name, self.date,self.time_sod, dt, 
                self.e_min[0], self.e_max[2], bg, 
                np.stack((data[:,0], data[:,1]+data[:,2]+data[:,3]), axis=-1)) 
        }

        str_info = f"Header info: {sc_name:10s}\nDate time: {self.date:8s} {self.time_sod:9.3f}\n"
        str_info += f"Emin Emax: {self.e_min[0]:5.1f} {self.e_max[2]:8.1f}\n"
        str_info += f"Bg lev.:{bg:5.1f}\nres: {dt:5.3f}\n"
        str_info += "t_min t_max: {:8.3f} {:8.3f}\nlc_zize: {:d}\n".format(
            self.lc['G1'].get_times()[0], self.lc['G1'].get_times()[-1], 
            self.lc['G1'].get_times().size)

        return str_info

    def get_lcs(self):
        return dict([ (s, self.lc[s]) for s in self.lc_names ])

    def get_info(self):
        return self.str_info


if __name__ == '__main__':

   hend = ipn_file('test_data/GRB20200415_T31685/20200415_30325_HEND.thr')
   hend_lc = hend.get_lc()
   print(hend_lc)

   hend_lc.set_bg_idx(4350, 4420)
   print(hend_lc.get_bg_sub_counts())

   gbm = gbm_thr_file('test_data/GRB20200415_T31685/GBM/GRB200415_GBM_1ms.thr')
   gbm_lcs = gbm.get_lcs()
   print (gbm_lcs['G1'])
   gbm_lcs['G1'].set_bg_idx(860, 975)
   print(gbm_lcs['G1'].get_bg_sub_counts())

   kw = kw_thr_file('test_data/GRB20200415_T31685/kw20200415_31681_rc2.thr',
       'test_data/GRB20200415_T31685/kw20200415_31681.info')
   kw_lcs = kw.get_lcs()
   print (kw_lcs['G1'])
   kw_lcs['G1'].set_bg_idx(860, 975)
   print(kw_lcs['G1'].get_bg_sub_counts())

