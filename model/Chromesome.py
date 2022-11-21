# -*- coding: utf-8 -*-
# @Author: Caoyh
# @Date:   2022-05-28 20:45:14
# @Last Modified time: 2022-06-15 12:02:42

import numpy as np
import copy
import sys

sys.path.append('.')
from utils import get_config
from utils import SingleLinkList
from utils import SingleReaction
# from kegg_helper.pykegg import split_equation
from log.logger import Logger

# 1. init the available pool
# 2. get a pair from pool based on the substrate by using the roulette strategy
# 3. update pool based on the products
# 4. 

class Chromesome(object):

    def __init__(self, args, cfg, log):

        self.args = args
        self.cfg = cfg
        self.log = log

        self.ob_sustrate = cfg[args.production]['ob_substrate']
        self.ob_product = cfg[args.production]['ob_product']
        self.abundant = cfg[args.production]['abundant']

        self.NrMax = cfg['init']['NrMax']
        self.all_pool = np.load(cfg['file_path']['mypool'], allow_pickle=True).item()
        
        self.ava_pool = dict()
        self.chrom = SingleLinkList()

    def get_a_chrom(self):
        count = 0
        while (1):
            if count == 10000:
                return False
            elif self.get_a_chrom_roulette():
                return True
            else:
                count += 1

    def init_ava_pool(self):
        """
        重新初始化 
        :return:
        """
        all_pool_cpds = list(self.all_pool.keys())
        for c in self.abundant + [self.ob_sustrate]:
            if c not in all_pool_cpds:
                continue
            else:
                temp = list(self.all_pool[c])
                self.ava_pool[c] = temp

    def update_ava_pool(self, Node):
        """
        根据给定的化合物 更新 pool
        :param compound:
        :return:
        """
        # substrates, products = split_equation(Node.reaction['equation'])
        reactants, products = Node.reaction['reactants'], Node.reaction['products']
        all_pool_cpds = list(self.all_pool.keys())
        for c in reactants + products:
            if c not in all_pool_cpds:
                continue
            temp = list(self.all_pool[c])
            if c in self.ava_pool: 
                continue
            else:
                self.ava_pool[c] = temp

    def update_pool_after_change(self):
        """ 
        After trim, crossover and mutation
        Update pool of chrom
        return self.ava_pool, self.compounds
        """
        self.init_ava_pool()
        cur = self.chrom._head 
        while (cur != None):
            self.update_ava_pool(cur)
            cur = cur.next
        # self.get_all_compounds()

    def roulette(self, NodeList):
        """轮盘赌算法
        :param NodeList: NodeList 是一个链表
        :return: 返回索引值
        """
        sum_val = 0
        for i in range(len(NodeList)):
            if NodeList[i].Tanimoto != 0 and NodeList[i].Tanimoto != None:
                sum_val = sum_val + NodeList[i].Tanimoto

        random_val = np.random.random()
        probability = 0  #累计概率
        for i in range(len(NodeList)):
            if NodeList[i].Tanimoto != 0 and NodeList[i].Tanimoto != None:
                probability = probability + NodeList[i].Tanimoto / sum_val
            if probability >= random_val:
                return i
            else:
                continue

    def get_a_chrom_roulette(self):
        # 按照底物相似性 根据轮盘赌算法 生成个体

        self.init_ava_pool()

        mychrom = SingleLinkList()
        index = self.roulette(list(self.ava_pool[self.ob_sustrate]))
        Node = copy.deepcopy(list(self.ava_pool[self.ob_sustrate])[index])
        self.update_ava_pool(Node)
        mychrom.append(Node)

        while (Node.P != self.ob_product and mychrom.length() < 100):
            index = self.roulette(list(self.ava_pool[Node.P]))
            Node = copy.deepcopy(list(self.ava_pool[Node.P])[index])
            self.update_ava_pool(Node)
            mychrom.append(Node)

        if Node.P == self.ob_product:
            # Get a chromesome sucessfully and trim the chrom
            chrom = self.trim_chrom(mychrom)
            if (chrom.length() < self.NrMax) and self.check_chrom_available(chrom):
                self.chrom = chrom
                return True 
        self.ava_pool.clear()
        return False

    def check_chrom_available(self, chrom):
        """检查是否存在同一个reaction出现两次的情况
        """
        head = chrom._head
        left = head
        if (chrom != None and chrom.length() > 2):
            while(left.next != None):
                if left.reaction['kegg_id'] == left.next.reaction['kegg_id']:
                    return False
                else:
                    left = left.next
        return True
            
    def trim_chrom(self, input_chrom):
        chrom = copy.deepcopy(input_chrom)

        head = chrom._head
        # 修剪情况: a->b b->c c->a a->d ===> a->d 
        if (chrom != None and chrom.length() > 2):
            dummy = SingleReaction(S=None, P=None, reaction={})
            dummy.next = head
            pre = dummy
            left = pre.next
            while(left != None):
                right = left.next
                while(right != None):
                    if left.S == right.S: 
                        pre.next = right
                        pre = dummy
                        left = pre.next
                        break
                    right = right.next
                pre = pre.next
                left = left.next
            return chrom
        else:
            return input_chrom  


if __name__ == '__main__':

    from main import parser
    args = parser()
    cfg = get_config()
    log = Logger('cache/log/chromesome.log')

    cla = Chromesome(args, cfg, log)
    cla.get_a_chrom()

    rxn_list = cla.chrom.travel()
    print(rxn_list)

    # # from fit_funcs import get_yield, get_gibbs

    # y = get_yield(rxn_list, cla.ob_product)
    # g = get_gibbs(rxn_list)

    # print(g)
