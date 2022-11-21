'''
Autor: caoyh
Date: 2022-11-17 16:37:08
LastEditTime: 2022-11-20 22:01:19
'''
# -*- coding: utf-8 -*-
# @Author: Caoyh
# @Date:   2022-05-28 20:40:43
# @Last Modified time: 2022-07-04 22:10:23

import os
import sys
import numpy as np
from pymoo.core.problem import ElementwiseProblem

sys.path.append('.')
from utils import SingleLinkList, SingleReaction
from fit_funcs import get_gibbs
from fit_funcs import get_yield

class SingleProblem(ElementwiseProblem):
    def __init__(self, args, cfg, log):
        super().__init__(n_var=1, n_obj=1, n_constr=0)
        self.args = args
        self.cfg = cfg
        self.log = log
        
        self.result_dir = self.cfg['file_path']['result_dir'] + self.args.project + '/'
        if not os.path.exists(self.result_dir): os.mkdir(self.result_dir)

    def _evaluate(self, x, out, *args, **kwargs):
        F1 = self.func_1(x)
        F2 = self.func_2(x)
        F3 = self.func_3(x)

        # 在评估的地方保存个体和适应度值
        with open(self.result_dir+'all_x_f.txt', 'a') as f:
            f.write(str([F1, F2, F3]))
            f.write('\n')
            f.write(str(x[0].chrom.get_rxn_list()))
            f.write('\n')
            travel_list = x[0].chrom.travel()
            for line in travel_list:
                for l in line:
                    f.write(l)
                    f.write('\t')
                f.write('\n')
            f.write('==' * 10)
            f.write('\n')

        fit = float(self.args.weight1) * F1 + float(self.args.weight2) * F2 + float(self.args.weight3) * F3
        out["F"] = np.array([fit], dtype=float)

    def func_1(self, Chrom):
        # length of the pathway
        # less is better

        self.log.logger.info('---get pathway lenght---')
        return Chrom[0].chrom.length()

    def func_2(self, Chrom):
        # Gibbs of the pathway MDF
        # less is better

        self.log.logger.info('---get pathway gibbs---')
        return get_gibbs(Chrom[0].chrom.get_rxn_list())

    def func_3(self, Chrom):
        # yelid of pathway
        # less is better
        rxn_list = Chrom[0].chrom.get_rxn_list()
        self.log.logger.info('---get pathway yield---')
        self.log.logger.info(Chrom[0].chrom.travel())
        self.log.logger.info(rxn_list)
        return -get_yield(rxn_list)
        # return -np.random.rand()


