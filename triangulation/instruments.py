"""
Contains main instrument classes
"""

import os
import sys
import re

import utils

class ipn_instrument():

    inst_list = ['BAT', 'SPI-ACS', 'HEND', 'MGNS']

    def __init__(self, name, path):
 
        self.name = name
        self.path = path

        self.thr_name = None
        self.eph_name = None
        
        if name == 'SPI-ACS':
            self.thr_name = utils.search_file(f"{path}/*_INT.thr")
            self.eph_name = utils.search_file(f"{path}/*_INT.eph")

        if name == 'HEND':
            self.thr_name = utils.search_file(f"{path}/*_HEND.thr")
            self.eph_name = utils.search_file(f"{path}/*_HEND_pos.txt")

        if name == 'BAT':
            self.thr_name = utils.search_file(f"{path}/*_BAT64.thr")
            self.eph_name = utils.search_file(f"{path}/*_Swift.eph")

        if name == 'MGNS':
            self.thr_name = utils.search_file(f"{path}/*_MGNS.thr")
            self.eph_name = utils.search_file(f"{path}/*_MGNS_pos_hor.txt")
        

    def __str__(self):

        return f"{self.name}:\nTHR: {self.thr_name}\nEPH: {self.eph_name}"


class konus_wind():

    def __init__(self, path):
        """
        Search for all kw files in path.
 
        Parameters
        ----------
        path            path to thr files, str

        """
        
        self.lst_files = []

        self.info_file = self.get_info_file(path)
        self.info_str, self.info = self.parse_info(self.info_file)

        self.lst_res, self.dic_files = self.get_thr_files(path)

        self.eph_name = utils.search_file(f"{path}/*_Wind.eph")

        #print(self.lst_res)
        #print(self.dic_files)
        #sys.exit()

    def __str__(self):
        str_out = "Konus-Wind files:\n"
        str_out += f"{self.info_file}\nTHR:\n"
        str_out += "\n".join([f'{res} us {self.dic_files[res]}' for res in self.lst_res])
        str_out += f'\nEPH: {self.eph_name}'
   
        return str_out

    def get_lc_file(self, res_ms):
        return self.dic_files[res_ms]

    def get_thr_files(self, path):

        lst_res = []
        dict_files = {}
        for s in os.listdir(path):
            if s.endswith('.thr'):
                
                m = re.search('_rc(\d+)\.thr', s)

                if m is None:
                    #print(f'Cannot split {s}')
                    continue

                res_us = int(m.group(1)) * 1000

                lst_res.append(res_us)
                dict_files[res_us] = os.path.join(path, s)

        return lst_res, dict_files

    def get_info_file(self, path):

        list_files = os.listdir(path)
        lst_ = list(filter(lambda x: x.endswith('.info'), list_files))
        if len(lst_) == 1:
            return os.path.join(path, lst_[0])
        else:
            return None
    
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

        return '\n'.join(lines), {
            'date':        lst_details[0], 
            'short_date':  lst_details[0][2:],
            'time_hhmmss': lst_details[2].lstrip('(').rstrip(')'),
            'short_time':  short_time,
            'time':        "{:s}".format(lst_details[1]),
            'det':         "{:s}".format(lst_details[3]),
            'name_long' :   "GRB{:s}_T{:s}".format(lst_details[0], short_time)
            }


class fermi_gbm():

    def __init__(self, path):
        """
        Search for all thr files in path.
 
        Parameters
        ----------
        path            path to thr files, str

        """

        self.path = path

        self.date_time_file = self.get_date_time(self.path)
        self.lst_res, self.dic_files = self.get_resolutions(self.path)

    def __str__(self):
        str_out = "Fermi-GBM files:\n"
        str_out += f"{self.date_time_file}\n"
        str_out += "\n".join([f'{res} us {self.dic_files[res]}' for res in self.lst_res])
        str_out += '\n'
        return str_out

    def get_date_time(self, path):

        file_name = os.path.join(path, 'fermi_date_time.txt')

        if not os.path.isfile(file_name):
            print(f"File not found {file_name}")
            return None

        return file_name

    def get_resolutions(self, path):

        lst_res = []
        dict_files = {}
        for s in os.listdir(path):
            if s.endswith('ms.thr'):

                m = re.search('_(\d+)ms', s)

                if m is None:
                    print(f'Cannot split {s}')
                    continue

                res = m.group(1)

                if res[0] == '0':
                    res_us = int(res) * 100
                else:
                    res_us = int(res) * 1000

                lst_res.append(res_us)
                dict_files[res_us] = os.path.join(path, s)

        return lst_res, dict_files

    def get_lc_file(self, res_ms):
        return self.dic_files[res_ms]


