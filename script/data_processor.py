# -*- coding: UTF-8 -*-

import argparse
import logging
import os
import sys
from copy import deepcopy
from functools import reduce
from operator import index
from os.path import split
from typing import AnyStr, Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib.pyplot import axis, tick_params, xticks

class DataProcessor:
    def __init__(self, in_config):
        config_dict = dict()
        config_dict = os.path.join(os.path.dirname(os.path.realpath(__file__)), in_config)
        with open(config_dict, 'r') as fin:
            config_dict = yaml.load(fin, Loader=yaml.FullLoader)

        self.exp_config = {'dir': '/home/rdaneel/my_exp_discovery3',
                           'ins_num': 25,
                           'time_limit': 60}

        self.map_config = {'name': 'random-32-32-20',
                            'scen': ['even', 'random'],
                            'sid': 0,
                            'num_of_ag': [20, 40, 60, 80]}

        self.scen_list: List[str] = list(config_dict['scen_dict'].keys())  # even, random
        self.solvers = [{'sol':'EECBS2', 'w':1.05, 'mth': -1, 'plot': ('grey', '+'), 'label': 'Vanilla EECBS'}]
        self.result = dict()

        # Plot parameters
        self.plot_config= {'width': 12, 'height': 9,
                          'ms': 18, 'lw': 2, 'ws': 35, 'mew': 2.5}
        # self.styles = {0: ('grey', '+'),
        #                1: ('purple', '>'),
        #                2: ('deepskyblue', 'o')}

    @staticmethod
    def readFile(f_name: str):
        if os.path.exists(f_name):
            return pd.read_csv(f_name)
        else:
            logging.error('{0} does not exist!'.format(f_name))
            sys.exit()

    def getFileDir(self, cbs_type, num_of_ag=None, in_map=None):
        if num_of_ag is None:
            num_of_ag = self.map_config['num_of_ag']
        if in_map is None:
            in_map = self.map_config['name']

        map_dir = os.path.join(self.exp_config['dir'], in_map)
        out_dir = os.path.join(map_dir, cbs_type)

        return out_dir

    def getFileName(self, scen_name, f_weight, cbs_type, sid=0, merge_th=-1, num_of_ag=None, in_map=None):
        if num_of_ag is None:
            num_of_ag = self.map_config['num_of_ag']
        if in_map is None:
            in_map = self.map_config['name']
        
        out_name =  in_map + '-' + scen_name + '-' + str(num_of_ag) + '-' + \
            '{0:.2f}'.format(f_weight) + '-' + str(sid)
        
        if merge_th == -1:
            out_name = out_name + '-' + cbs_type + '.csv'
        else:
            out_name = out_name + '-' + str(merge_th) + '-' + cbs_type  + '.csv'
        
        return out_name

    def getCSVIns(self, scen_name, f_weight, cbs_type, sid=0, merge_th=-1, num_of_ag=None, in_map=None):
        return self.readFile(os.path.join(
            self.getFileDir(cbs_type, merge_th, num_of_ag, in_map),
            self.getFileName(scen_name, f_weight, cbs_type, sid, merge_th, num_of_ag, in_map)))

    def getIns(self, scen_name, f_weight, cbs_type, merge_th=-1):
        sid_list = self.sid_dict[scen_name]
        df = self.getCSVIns(scen_name, f_weight, cbs_type, min(sid_list), merge_th)
        tmp_key: Tuple[str, float, int, int, str] = (scen_name, f_weight, min(sid_list), merge_th, cbs_type)
        out_dict: Dict[Tuple[str, float, int, int, str], pd.Series] = {tmp_key: df['instance name']}
        
        for sid in sid_list:
            if sid > min(sid_list):
                tmp_df = self.getCSVIns(scen_name, f_weight, cbs_type, sid, merge_th)
                out_dict[(scen_name, f_weight, sid, merge_th, cbs_type)] = tmp_df['instance name']
                df = df.append(tmp_df, ignore_index=True)
        
        return df, out_dict

    def getInsByDF(self, f_ins, scen_name, f_weight, cbs_type, merge_th=-1):
        sid_list = self.sid_dict[scen_name]
        df = self.getCSVIns(scen_name, f_weight, cbs_type, min(sid_list), merge_th)
        cond = df['instance name'].isin(f_ins) == False
        df = df.drop(df[cond].index)
        
        for sid in sid_list:
            if sid > min(sid_list):
                tmp_df = self.getCSVIns(scen_name, f_weight, cbs_type, sid, merge_th)
                cond = tmp_df['instance name'].isin(f_ins) == False
                tmp_df = tmp_df.drop(tmp_df[cond].index)
                df = df.append(tmp_df, ignore_index=True)

        return df
    
    def plot_best_case(self):
        
        type_list = [(-1, 'ECBS_rt30'),  (-1, 'ECBS_flex_rt30')]  # socs2021
        # type_list = [(-1, 'ECBS'), (-1, 'ECBS_flex'), (-1, 'ECBS_rt5'),  (-1, 'ECBS_flex_rt5')]  # socs2021

        if self.plot_common:
            self.common_ins = self.get_common_ins_best()

        self.result = dict()  # Reset the dictionary
        for tp in type_list:
            if tp[1] in self.opt_list:  # optimal cases
                if tp[1] == 'CBS':
                    self.cal_iteration(tp[1], [1.00], [-1], in_mode)
                else:
                    self.cal_iteration(tp[1], [1.00], [tp[0]], in_mode)
                    
            elif tp[1] in self.subopt_list:  # sub-optimal cases
                print(tp[1])
                if tp[1] == 'ECBS' or tp[1] ==  'ECBS_flex' or tp[1] == 'ECBS_flex2' or tp[1] == 'ECBS_rt10' or tp[1] == 'ECBS_rt20' or tp[1] == 'ECBS_rt30' or tp[1] == 'ECBS_rt40':
                    self.cal_iteration(tp[1], self.focal_w, [-1], in_mode)
                else:
                    self.cal_iteration(tp[1], self.focal_w, [tp[0]], in_mode)

        # Plot figure
        f = plt.figure(num=None, figsize=(self.fig_h, self.fig_w), dpi=80, facecolor='w', edgecolor='k')
        ax = f.add_subplot(111)
        ax.tick_params(labelright=False)

        for tp in type_list:
            if tp[1] in self.opt_list:  # optimal cases
                if tp[1] in self.ma_list:
                    self.plot_y_val(tp[1], [1.00], [tp[0]], in_mode)
                else:
                    self.plot_y_val(tp[1], [1.00], [-1], in_mode)

            elif tp[1] in self.subopt_list:
                if tp[1] in self.ma_list:
                    self.plot_y_val(tp[1], self.focal_w, [tp[0]], in_mode)
                else:
                    self.plot_y_val(tp[1], self.focal_w, [-1], in_mode)

        # ax.yaxis.grid(True)
        params = {
            'legend.fontsize': self.legend_size,
            'legend.handlelength': 2,
        }

        out_file: str = "temp.png"
        if in_mode == 0:  # num-success rate
            plt.yticks(np.arange(0, 105, 20))
            plt.ylabel('Success rate (%)', fontsize=self.label_size)
            out_file = self.map_name + '_' + str(self.num_of_agent) + '_num_success.png'

        elif in_mode == 1:  # num-runtime
            plt.ylabel('Average runtime (sec)', fontsize=self.label_size)
            out_file: str = self.map_name + '_' + str(self.num_of_agent) + '_num_runtime.png'

        elif in_mode == 2:  # num-max MA size
            _, end = ax.get_ylim()
            plt.yticks(np.arange(0, np.ceil(end), 1))
            plt.ylabel('Average max size of meta-agent', fontsize=self.label_size)
            out_file = self.map_name + '_' + str(self.num_of_agent) + '_num_maxMA.png'

        elif in_mode == 3:  # avg # generated nodes
            _, end = ax.get_ylim()
            # plt.yticks(np.arange(0, np.ceil(end), 10000))
            plt.ylabel('Average generated nodes', fontsize=self.label_size)
            out_file = self.map_name + '_' + str(self.num_of_agent) + '_num_gnode.png'

        else:
            logging.error("in_mode is required!")
            exit(1)

        plt.rcParams.update(params)
        plt.xticks(self.num_list)
        plt.xlabel('#Agents', fontsize=self.label_size)
        ax.tick_params(labelsize=self.num_size, width=3)
        ax.legend(fontsize=self.legend_size)
        plt.subplots_adjust(left=0.15, bottom=0.15, right=0.98, top=0.98, wspace=0.2, hspace=0.2)

        # ax.legend(fontsize=self.legend_size, loc='lower center', bbox_to_anchor=(0.42, -0.57),
        #     labelspacing=0.1, columnspacing=0.4, ncol=2, borderpad=0.1, borderaxespad=1.1, handletextpad=0.05)
        # plt.subplots_adjust(left=0.165, bottom=0.32, right=0.98, top=0.98, wspace=0.2, hspace=0.2)
        
        plt.savefig(os.path.join(os.path.dirname(os.path.realpath(__file__)), out_file))
        plt.show()
        return

