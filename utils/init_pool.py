'''
Description: 
Autor: caoyh
Date: 2022-11-16 12:59:17
LastEditTime: 2022-11-17 10:41:33
'''
import os
import sys
import json
import pdb
from rdkit import DataStructs
from rdkit import Chem
import re, time
import numpy as np
import pandas as pd
import sys
sys.path.append('.')

from utils import get_config
from utils import read_json
from utils import SingleReaction
from kegg_helper.pykegg import split_equation

# class SingleReaction(object):
#     def __init__(self, S, P, reaction, Tanimoto=0.0, next=None):
#         self.S = S
#         self.P = P
#         self.reaction = reaction
#         self.Tanimoto = Tanimoto
#         self.next = next
   

class MyPool(object):
    def __init__(self, cfg) -> None:
        self.cfg = cfg
        self.rxn_dict = read_json(cfg['file_path']['all_rxn_dict'])
        self.cpd_dict = read_json(cfg['file_path']['all_cpd_dict'])

        self.save_path = self.cfg['file_path']['data_dir']
        if not os.path.exists(self.save_path): os.mkdir(self.save_path)
        
        self.pool = dict()

        self.err_cpds = set()
        self.g_cpds = set()

        self._get_pool()


    def get_Tanimoto(self, s,p):
        # check cpd name
        def check_cpd_name(cpd_name):
            pat = '^C\d{5}'
            res = re.findall(pat, cpd_name)
            if len(res) != 0:
                return True
            else:
                return False
                
        if not(check_cpd_name(s) and check_cpd_name(p)):
            return 0

        try:
            s_smile = self.cpd_dict[s]['smile']
            s_mol = Chem.MolFromSmiles(s_smile)
            p_smile = self.cpd_dict[p]['smile']
            p_mol = Chem.MolFromSmiles(p_smile)
            mols = [s_mol, p_mol] 
            fps = [Chem.RDKFingerprint(x) for x in mols]
            t = DataStructs.FingerprintSimilarity(fps[0], fps[1])
        except Exception as e:
            print(e)
            print("s is %s and p is %s" % (s, p))
            t = 0
        return t
 

    def _get_invalid_cpd_list(self):
        def get_cpd_csv(file_path):
            cpds = pd.read_csv(file_path, header=None)
            cpds = cpds.drop_duplicates()
            cpds_list = cpds[0].values.tolist()
            return cpds_list
        cofactors_list = get_cpd_csv(self.save_path+'cofactors.csv')
        error_cpds_list = get_cpd_csv(self.save_path+'error_cpds.csv')
        return cofactors_list+error_cpds_list 

    def _get_pool(self):   
        invalid_cpds = self._get_invalid_cpd_list()

        KEYS = list(self.rxn_dict.keys())
        len_pool = 0
        for R in KEYS:
            s_list, p_list = split_equation(self.rxn_dict[R]["equation"])
                
            for s in s_list:
                if s in invalid_cpds: continue
                for p in p_list:
                    if p in invalid_cpds: continue
                    srn_front = SingleReaction(s, p, self.rxn_dict[R], Tanimoto=self.get_Tanimoto(s,p))
                    self.pool.setdefault(s, set()).add(srn_front)
                    len_pool += 1

            for p in p_list:
                if p in invalid_cpds: continue
                for s in s_list:
                    if s in invalid_cpds: continue
                    srn_front = SingleReaction(p, s, self.rxn_dict[R], Tanimoto=self.get_Tanimoto(p,s))
                    self.pool.setdefault(p, set()).add(srn_front)
                    len_pool += 1
            print("Reaction %s is done!" %(R))

        pool_all_file = self.save_path + 'MYPOOL_' + time.strftime("%Y%m%d") + '.npy'   
        np.save(pool_all_file, self.pool)

        print('===== Done! ======')
        print("There are %s pairs in the pool" %(len_pool))

def main():
    cfg = get_config() 
    MyPool(cfg)
    


if __name__=='__main__':
    main()
    

