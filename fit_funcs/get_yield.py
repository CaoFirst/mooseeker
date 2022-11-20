# -*- coding: utf-8 -*-
# @Author: Caoyh
# @Date:   2022-06-29 11:31:07
# @Last Modified time: 2022-07-29 16:19:38
import sys
sys.path.append('.')

import os
import json
from pathlib import Path
import cobra
# import cobra.test
from cobra.io import read_sbml_model
from cobra import Reaction, Metabolite
import pubchempy as pcp
from kegg_helper.pykegg import *
from utils import get_config
from utils import read_json

cfg = get_config()

dir = cfg['file_path']['yield_dir']
dir_data = Path(dir)

if os.path.exists(dir + 'DatabaseVersions.csv'):
    os.remove(dir + 'DatabaseVersions.csv')
   
rxn_dict = read_json(cfg['file_path']['all_rxn_dict'])
cpd_dict = read_json(cfg['file_path']['all_cpd_dict'])

# def parse_equation(equation): # -1: substrates  1: products
#     """Parses a KEGG reaction string into a tuple of lists of tuples."""
#     # Split the equation into lists
#     # eq_list = re.split(" +", equation)
#     # reactants = eq_list[0:eq_list.index("<=>")]
#     # products = eq_list[eq_list.index("<=>")+1:]
#     reactants, products = split_equation(equation)
#     output = {}
#     for s in reactants:
#         output[s] = -1

#     for p in products:
#         output[p] = 1
#     return output

    # s = 1
    # for segment in reactants:
    #     if re.match("^[0-9]+$", segment):
    #         s = int(segment)
    #         continue
    #     if segment == "+":
    #         continue
    #     else:
    #         output[segment] = -1 * s
    #         s = 1
    # p = 1
    # for segment in products:
    #     if re.match("^[0-9]+$", segment):
    #         p = int(segment)
    #         continue
    #     if segment == "+":
    #         continue
    #     else:
    #         output[segment] = p
    #         # output.append((s,segment))
    #         p = 1
    # return output


def add_cpd2rxn(model, rxn):

    r = rxn_dict[rxn]
    rxn_coff = parser_equation_coef(r['equation'])
    reactants, products = split_equation(r['equation'])
    
    names = locals()
    
    meta_dict = {}
    for cpd in reactants:
        cpd = re.search(r'.*(C\d{5}).*', cpd).group(1)
        # c = MyCompound(cpd)
        c = cpd_dict[cpd]
        
        if c['bigg_id']!=None:
            for bid in c['bigg_id']:
                if Metabolite(bid+'_c') in model.metabolites:
                        # meta_dict[model.metabolites.get_by_id(bid+'_c')] = -1
                        meta_dict[model.metabolites.get_by_id(bid+'_c')] = rxn_coff[cpd]
                        break

        else:
            # print('Something is error. Mannully add compound %s to metabolites ' %(cpd))
            c_name = c['kegg_id'] + '_' + str(c['pubchem'])
            names[c_name] = Metabolite(id = cpd + '_c', 
                                        compartment='c', 
                                        name=c_name, 
                                        formula=c['formula']) 
            # meta_dict[names.get(c_name)] = -1 
            meta_dict[names.get(c_name)] = rxn_coff[cpd]
      
    for cpd in products:
        cpd = re.search(r'.*(C\d{5}).*', cpd).group(1)
        # c = MyCompound(cpd)
        c = cpd_dict[cpd]

        if c['bigg_id']!=None:
            for bid in c['bigg_id']:
                if Metabolite(bid+'_c') in model.metabolites:
                        # meta_dict[model.metabolites.get_by_id(bid+'_c')] = 1
                        meta_dict[model.metabolites.get_by_id(bid+'_c')] = rxn_coff[cpd]
                        break

        else:
            # print('Something is error. Mannully add compound %s to metabolites ' %(cpd))
            c_name = c['kegg_id'] + '_' + str(c['pubchem'])
            names[c_name] = Metabolite(id = cpd + '_c', 
                                        compartment='c', 
                                        name=c_name, 
                                        formula=c['formula']) 
            # meta_dict[names.get(c_name)] = 1  
            meta_dict[names.get(c_name)] = rxn_coff[cpd]

    return meta_dict

def get_yield(rxn_list):
    
    original = read_sbml_model(cfg['file_path']['host_sbml'])
    model = original.copy()

    medium = model.medium
    medium["EX_o2_e"] = 20.0
    medium["EX_glc__D_e"] = 20.0
    model.medium = medium
    
    wt_growth = model.optimize()
    max_growth = wt_growth.objective_value
    min_growth = 0.8*max_growth

    rxn_names = locals()

    obj_rxn_id = None

    for rxn in rxn_list:

        # r = MyReaction(rxn)
        if rxn not in rxn_dict.keys():
            return 0
            
        r = rxn_dict[rxn]

        # First Judge rxn is already in model.reactions
        if r['bigg_id']!=None:
            for rid in r['bigg_id']:
                if Reaction(rid) in model.reactions: 
                    # print('%s is already in model' %(rxn))
                    obj_rxn_id = rid
                    break
                
        else:   
            rxn_names[rxn] = cobra.Reaction(rxn)
            rxn_names[rxn].lower_bound = -1000
            rxn_names[rxn].upper_bound = 1000

            rxn_names[rxn].add_metabolites(add_cpd2rxn(model, rxn))
            
            model.add_reactions([rxn_names[rxn]])
            
            # print('%s is done!' %(rxn))
            obj_rxn_id = rxn

    # model.add_boundary(model.metabolites.get_by_id(), type='demand')

    with model:
        medium = model.medium
        
        model.objective = model.reactions.get_by_id(obj_rxn_id)
        solution = model.optimize()
        max_biomass = solution.objective_value 

        if max_biomass < min_growth: 
            print("This pathway is not feasible!")
            return -999 

        else:
            if model.reactions.get_by_id('EX_glc__D_e').flux==0:
                print("This pathway cannot get EX_glc__D_e, Maximum theoretical yield = 0")
                return 0
            else:
                maximum_yield = max_biomass / (-1*(model.reactions.get_by_id('EX_glc__D_e').flux)) #Target production[mmol/gDW*h]
                # print('Maximum productivity =', max_biomass, 'mmol/gDW*h')
                
                if maximum_yield == 50: 
                    print("This pathway EX_glc__D_e is 50, so the Maximum theoretical yield = 0")
                    return 0.01
                else: 
                    print('Maximum theoretical yield =', maximum_yield, 'mmol-Dgly/mmol-glc')
                    return maximum_yield
    

def main():
    # rxn_list = ['R07265', 'R11306', 'R04411', 'R04410', 'R09127', 'R06973']
    rxn_list = ['R01788', 'R02737', 'R02780', 'R02628', 'R00724']
    print(get_yield(rxn_list))

if __name__ == "__main__":

    main()


    