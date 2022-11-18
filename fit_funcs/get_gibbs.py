'''
Description: 
Autor: caoyh
Date: 2022-09-16 17:02:26
LastEditTime: 2022-11-17 17:04:21
'''
# -*- coding: utf-8 -*-
# @Author: Caoyh
# @Date:   2022-06-01 15:37:55
# @Last Modified time: 2022-07-29 16:22:39


import json
import pandas as pd
# Add repository root to the path
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))

sys.path.append('.')

# Import the script to be tested
from kegg_helper.pykegg import split_equation
from fit_funcs.gibbs.mdf import *
from utils import read_json
from utils import get_config
cfg = get_config()


def get_gibbs(rxn_list):
    dG_dict = read_json(cfg['file_path']['dG_dict'])
    rxn_dict = read_json(cfg['file_path']['all_rxn_dict'])
    cpd_dict = read_json(cfg['file_path']['all_cpd_dict'])
    cons_list = pd.read_csv(cfg["file_path"]['gibb_concentration'], sep='\t', header=None).to_dict('list')

    reaction = []
    drGs = []
    con = []
    cpb_list = []

    for rxn in rxn_list:
        reaction += [rxn + '\t' + rxn_dict[rxn]['equation']]
        try: 
            # 如果没有这个值 或者没有这个反应，则认为为0
            drGs += [rxn + '\t' + str(dG_dict[rxn]['dG_Mean'])] 
        except Exception as e:
            print(e)
            drGs += [rxn + '\t' + str(0)] 

        s_list, p_list = split_equation(rxn_dict[rxn]['equation'])

        for s in s_list:
            if s in cpb_list: continue
            else: cpb_list.append(cpd_dict[s])
        for p in p_list:
            if p in cpb_list: continue
            else: cpb_list.append(cpd_dict[p])

    for cpb in cpb_list:
        if cpb in cons_list[0]:
            idx = cons_list[0].index(cpb)
            # min = cons_list[1][idx]
            # max = cons_list[2][idx]
            con += [cpb + '\t' + str(cons_list[1][idx]) + '\t' + str(cons_list[2][idx])]

    reactions_text = "\n".join(reaction) + "\n"
    drGs_text = "\n".join(drGs) + "\n"
    cons_text = '\n'.join(con)


    # Load stoichiometric matrix, drGs and constraints
    S = read_reactions(reactions_text)
    drGs = read_reaction_drGs(drGs_text)
    constraints = read_constraints(cons_text)
    # Calculate MDF inputs
    c = mdf_c(S)
    A = mdf_A(S)
    b = mdf_b(S, drGs, constraints)
    # Perform MDF
    mdf_result = mdf(c, A, b)

    return mdf_result.fun


if __name__=='__main__':
    
    # M00001
    # rxn_list = ['R01786', 'R02740', 'R04779', 'R01070', 'R07159', 'R01518', 'R00658', 'R00200']
    

    # M00002
    rxn_list = ['R07265', 'R11306', 'R04411', 'R04410', 'R09127', 'R06973', 'R00907']
    
    result = get_gibbs(rxn_list)

    print(result)












