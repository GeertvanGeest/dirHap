from __future__ import print_function # for compatibilitu with python 2

import os
import shlex # what does this do?
import subprocess
from argparse import ArgumentParser

import numpy as np
import pandas as pd

import itertools
import numpy as np
import pandas as pd
from scipy.stats import multinomial
from scipy.stats import binom
from scipy.stats import chisquare

from math import *

def calc_error(n, total=100, seq_error=0.02):
    return(binom.cdf(n, total, seq_error))

def calc_freq(n, total=100):
    return (int(n) / total)

class multi_nomial_likelihood:
    def __init__(self, haps, counts, ploidy, error=0.02):
        
        self.n = sum(counts)
        self.counts = counts
        self.ploidy = ploidy
        self.haps = np.array(haps)
        self.no_haps = len(self.haps)

        if error < 0.05:
            self.error = error
        else:
            print("Error percentage is too big")
            self.error = 0.05
        
        self.phasings = self.calculate_all_phasings()
        self.prob_array = self.calc_prob_array()
        
        # print(self.n, self.counts, self.ploidy, self.haps, self.no_haps)
        
    def calc_prob_array(self):
        """
        Calculate probabilities array
        """
        
        prob_array = []
        
        for index, phase in enumerate(self.phasings):    

            # hap_indices = list(set(self.phasings_indices[index]))

            # prob = np.array([phase.count(i) for i in self.haps[hap_indices]]) / self.ploidy
            prob = np.array([phase.count(i) for i in self.haps]) / self.ploidy

            if len(list(set(phase))) == 1:
                
                prob[prob.argmax()] =  prob[prob.argmax()] - self.error
                prob[prob.argmin()] =  self.error

            prob_array += [prob.tolist()]
        
        return prob_array


    
    def calculate_all_phasings(self):
        phasings_full = list(itertools.combinations_with_replacement(self.haps, self.ploidy))
        # phasings_indices = list(itertools.combinations_with_replacement(range(self.no_haps), self.ploidy))
        return phasings_full #, phasings_indices
    
    def calculate_likelihood(self):
        """
        Calculate the likelihood in log space. 
        """
        
        # self.probabilities_no_log = np.array([multinomial.pmf(self.counts, n=self.n, p=i) for i in self.prob_array])
        # self.log_probabilities = np.array([multinomial.logpmf(self.counts, n=self.n, p=i) for i in self.prob_array])
        # sum_prob = self.probabilities_no_log.sum()
        # self.likelihoods_no_log = self.probabilities_no_log / sum_prob
        
        # in log space
        self.log_probabilities = np.array([multinomial.logpmf(self.counts, n=self.n, p=i) for i in self.prob_array])
        sum_prob = self.log_probabilities[~np.isinf(self.log_probabilities)].sum()
        self.likelihoods = self.log_probabilities - sum_prob

        # get best solution
        self.best_solution = self.likelihoods.argmax()
        self.selected_phasing = self.phasings[self.best_solution]
        self.mle = self.likelihoods[self.best_solution]

        t = sorted(self.likelihoods)

        if np.isinf(t[-2]) == False:
            # self.mle_ratio = np.log((t[-1] / t[-2]))
            self.mle_ratio = t[-1] - t[-2]

            # self.mle_ratio = t[-1] / t[-2]))
        else:
            self.mle_ratio = 1e6


        if(np.isinf(self.mle_ratio)):
            print(t[-3:])
            self.mle_ratio = 1e6
        
        # chi-square
        # self.chi_square = np.array([chisquare(self.counts, f_exp=np.array(i)*self.n)[1] for i in self.prob_array])
        # self.chi_square_ml = self.chi_square / self.chi_square.sum()

class perform_haplotyping:
    def __init__(self, df, freq=0.05, prob_true=0.99, ploidy=4, error=0.02):

        self.df = df # this dataframe contains the raw haplotype counts. 

        # set thresholds for filtering
        self.freq_thres = freq
        self.prob_true_thres = prob_true
        self.ploidy = ploidy
        self.error = error


    def run(self):
        """
        Run the haplotype analysis over a single group
        """

        self.add_prob_true() # calculate prob of errrounous haplotype
        self.add_freq() # add frequency
        self.add_indices() # add indices
        self.add_single_dosage() # calculate the bi-allelic score per haplotype
        self.add_haplotype_group() # calculate the mnlr

    def add_prob_true(self):
        true_prob = self.df['count'].apply(calc_error, args=(self.df['count'].sum(),))
        self.df.loc[self.df.index.values,'prob_true'] = true_prob

    def add_freq(self):
        freq = self.df['count'].apply(calc_freq, args=(self.df['count'].sum(),))
        self.df.loc[self.df.index.values,'freq'] = freq

    def add_indices(self):
        """
        Defines which haplotypes should be used for dosage estimatiion with the multi-allelic model. 
        """
        indices = (self.df.prob_true > self.prob_true_thres) & (np.array(self.df.freq) > self.freq_thres)
        self.df.loc[self.df.index.values,'indices'] = indices

    def add_haplotype_count(self, phasing):

        def hap(x, phasing):
            if phasing.count(x) == 0:
                return np.nan
            else:
                return phasing.count(x)
            
        dosage = self.df['seq'].apply(hap, args=(phasing,))
        self.df.loc[self.df.index.values, 'MA_dosage'] = dosage

    def add_haplotype_group(self):
        """
        Calculate the multi-nomial model
        """
        
        haps = self.df[self.df.indices]['seq'].tolist()
        cnts = self.df[self.df.indices]['count'].tolist()
        
        if len(haps) > 1:
            ml = multi_nomial_likelihood(haps, cnts, self.ploidy, error=self.error)
            ml.calculate_likelihood()

            self.df.loc[self.df.index.values, 'MA_dosage'] = np.nan    
            self.add_haplotype_count(ml.selected_phasing)

            self.df.loc[self.df.index.values, 'MA_LR_ratio'] = ml.mle_ratio  
        elif len(haps) == 1:
            self.df.loc[self.df.index.values, 'MA_dosage'] = self.ploidy   
            self.df.loc[self.df.index.values, 'MA_LR_ratio'] = np.Inf           
        else:
            self.df.loc[self.df.index.values, 'MA_dosage'] = np.nan    
            self.df.loc[self.df.index.values, 'MA_LR_ratio'] = np.nan

    def calc_dosage(self, n, total):
        """Very simple biniomail model to caluclate dosage per haplotype"""

        expected = np.arange(0,1.1, 1/self.ploidy)
        prob = binom.pmf(n,total,expected)

        selected = prob.argmax()
        ratio = np.log(sorted(prob)[-1] / sorted(prob)[-2])

        return pd.Series({'BA_dosage':selected, 'BA_LR_ratio':ratio})


    def add_single_dosage(self):
        """
        A vs depth, haplotype expected fraction. 
        """
        
        dos = self.df['count'].apply(self.calc_dosage, args=(self.df['count'][self.df['indices']==True].sum(),))
        
        self.df.loc[self.df.index.values,'BA_dosage'] = dos.BA_dosage
        self.df.loc[self.df.index.values,'BA_LR_ratio'] = dos.BA_LR_ratio

        self.df.loc[self.df.indices==False,'BA_dosage'] = np.nan
        self.df.loc[self.df.indices==False,'BA_LR_ratio'] = np.nan
        

def main():

    argparser = ArgumentParser()

    file_path_group = argparser.add_argument_group(title='File structure')
    file_path_group.add_argument('--haplotype_data', type=str, help='File with haplotypes', required=True)
    file_path_group.add_argument('--output_file', type=str, help='Output file name', default="output.csv")
    # file_path_group.add_argument('--output_folder', type=str, help='Output folder',  default="./output/")

    run_group = argparser.add_argument_group(title='Run command args')
    run_group.add_argument('--ploidy', type=str, help='Ploidy', default='4')
    run_group.add_argument('--prob_true', type=str, help='Probability to consider a single haplotype as true haplotype', default='0.99')
    run_group.add_argument('--freq_cutoff', type=str, help='count freq cutoff', default='0.05')
    run_group.add_argument('--error', type=str, help='expected error percentage of sequencing reads', default='0.02')

    args = argparser.parse_args()

    # specify the thresholds
    ploidy = int(args.ploidy)
    prob_true = float(args.prob_true)
    freq_cutoff = float(args.freq_cutoff)
    error = float(args.error)

    # read in the haplotypes hapid	region	sample	count	chrom	start	end	seq
    # THis solution does not scale particulary well. We could read this in chunks and process it parrallel. 
    haplotype_table = pd.read_csv(args.haplotype_data, sep='\t', names=["sample",  "chrom", "start", "end", "seq", "count"])
    print(haplotype_table)
    haplotype_list = haplotype_table.groupby(["sample", "chrom", "start", "end"])

    dfs = []

    # As the input is defined by sample and region we have a single df that is parsed for haplotype analysis. 
    for index, group in enumerate(haplotype_list):
        # print(index, group[0])
        
        if index % 500 == 0:
            print(index)
        
        haplotyper = perform_haplotyping(group[1], freq=freq_cutoff, prob_true=prob_true, ploidy=ploidy, error=error)
        haplotyper.run()

        print(haplotyper.df)

        dfs += [haplotyper.df]


    output_file = args.output_file
    dfs = pd.concat(dfs)
    dfs = dfs.rename(columns={'indices': "used"})   # change to used, more descriptive
    dfs = dfs.round({"freq":2, "prob_true":3, "BA_LR_ratio":1,  "MA_LR_ratio":1}) # more nice output
    dfs.to_csv(output_file, index=None, sep='\t')

if __name__ == '__main__':
    main()