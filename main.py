'''
Autor: caoyh
Date: 2022-11-18 09:37:01
LastEditTime: 2022-11-20 16:43:19
'''
import os
import sys
import dill
import time
import numpy as np
from tabnanny import verbose
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.optimize import minimize
from argparse import ArgumentParser
from utils import get_config
from model import MultiProblem
from model import SingleProblem
from model.BiOperators import *
sys.path.append('.')

from log.logger import Logger

def parser():
    parser = ArgumentParser('MooSeeker')
    parser.add_argument('-alg', '--algorithm', default='multi', help='Multi or Single objective algorithm to use') ## multi OR single 
    parser.add_argument('-p', '--production', default='vanillin', help='Glycolysis or Vanillin to produce') # glycolysis or vanillin
    parser.add_argument('-w1', '--weight1', default=0.33, help='The weight of length for single objective algorithm')
    parser.add_argument('-w2', '--weight2', default=0.33, help='The weight of gibbs for single objective algorithm')
    parser.add_argument('-w3', '--weight3', default=0.33, help='The weight of yield for single objective algorithm')
    args = parser.parse_args()
    
    args.project = args.algorithm + '_' + args.production

    return args


def train(args, cfg):
    if args.algorithm == 'multi':
        algorithm = NSGA2(pop_size=cfg['alg']['pop_size'],
                    sampling=BioSampling(),
                    crossover=BioCrossover(),
                    mutation=BioMutation(),
                    eliminate_duplicates=BioDuplicateElimination())

        res = minimize(problem=MultiProblem(cfg),
                    algorithm=algorithm,
                    termination=('n_gen', cfg['alg']['gen']),
                    seed=1,
                    callback=BioCallback(),
                    verbose=True,
                    save_history=True)

    if args.algorithm == 'single':
        algorithm = GA(pop_size=cfg['alg']['pop_size'],
                        sampling=BioSampling(),
                        crossover=BioCrossover(),
                        mutation=BioMutation(),
                        eliminate_duplicates=BioDuplicateElimination())

        res = minimize(problem=SingleProblem(cfg),
                    algorithm=algorithm,
                    termination=('n_gen', cfg['alg']['gen']),
                    seed=1,
                    callback=BioCallback(),
                    verbose=True,
                    save_history=True)

    # save checkpoint 
    with open(cfg['file_path']['result_dir']+"checkpoint_"+time.strftime("%Y%m%d"), "wb") as f:
        dill.dump(algorithm, f)

    # print("Best solution found: \nF = %s\nX = %s" % (res.F, res.X[0][0].chrom.travel()))
    # print(res.X[0][0].chrom.get_rxn_list())

    # ------------------ save result ------------------

    # save history
    hist = res.history
    n_evals = []             # corresponding number of function evaluations\
    hist_X = []              # the individuals in each generation
    hist_F = []              # the objective space values in each generation
    hist_cv = []             # constraint violation in each generation
    hist_cv_avg = []         # average constraint violation in the whole population

    for algo in hist:

        # store the number of function evaluations
        n_evals.append(algo.evaluator.n_eval)

        # retrieve the optimum from the algorithm
        opt = algo.opt

        # store the least contraint violation and the average in each population
        hist_cv.append(opt.get("CV").min())
        hist_cv_avg.append(algo.pop.get("CV").mean())

        # get all the individuals 
        hist_X.append(opt.get("X").tolist())

        # filter out only the feasible and append and objective space values
        feas = np.where(opt.get("feasible"))[0]
        hist_F.append(opt.get("F")[feas].tolist()) 
    
    # save to file 
    print('Save result and history !')
    np.savez(cfg['file_path']['result_dir']+'result.npz', X=res.X, F=res.F, CV=res.CV, G=res.G, dtype=object)
    np.savez(cfg['file_path']['result_dir']+'history.npz', n_evals=n_evals, X=hist_X, F=hist_F, cv=hist_cv, cv_avg=hist_cv_avg, dtype=object)
    print('=====Done!=====')
    

if __name__ == '__main__':
    # parser = ArgumentParser('MooSeeker')
    # parser.add_argument('-alg', '--algorithm', default='multi', help='Multi or Single objective algorithm to use') ## multi OR single 
    # parser.add_argument('-p', '--production', default='vanillin', help='Glycolysis or Vanillin to produce') # glycolysis or vanillin
    # args = parser.parse_args()
    # args.project = args.algorithm + '_' + args.production

    # print(args)
    args = parser()
    cfg = get_config()
    log_file = cfg['file_path']['log_dir'] + args.project+ '_' + time.strftime("%Y%m%d")+'.log'
    log = Logger(log_file,level='info')

    log.logger.info(args)

    log.logger.info(type(args.weight3))


    # config_file = 'config/' + args.algorithm + '_' + args.production + '.yaml'
    # assert os.path.exists(config_file)
    
    # cfg = get_config(config_file)
    # train(args, cfg)

    