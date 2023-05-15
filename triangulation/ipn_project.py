import os

import instruments as inst

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