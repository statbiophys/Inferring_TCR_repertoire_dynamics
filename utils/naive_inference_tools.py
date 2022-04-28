import glob
import pandas as pd
import numpy as np
from functools import partial
import re
from full_inference_tools import *


def naive_inference_fthresh(f_thresh, path, filename1, filename2, colnames1, colnames2):
    
    paras_MP_alpha_45_NB = ,^.array([ -1.9282954 ,   0.26787789,   1.02804643, -10.13950223])
    paras_1 = paras_MP_alpha_45_NB #no need of paras to modify
    paras_2 = paras_MP_alpha_45_NB #no need of paras to modify
    Model = Inference_dynamics(path, paras_1, paras_2, filename1, filename2, colnames1, colnames2)
    df = Model.df
    df_copy = df.copy()
    df_bis = df_copy[df_copy['Clone_count_2'] > 0]
    df_bis = df_bis[df_bis['Clone_fraction_1'] >= f_thresh]
    
    df_bis['log_f_1'] = np.log(df_bis['Clone_fraction_1'])
    df_bis['log_f_2'] = np.log(df_bis['Clone_fraction_2'])
    
    s_dist = df_bis['log_f_2']-df_bis['log_f_1']
    
    N_obs = len(s_dist)
    A_t = np.sum(s_dist)/N_obs
    B_t = (1/N_obs)*np.sum((s_dist - (A_t))**2)
    
    #A_t = np.mean(s_dist)
    #B_t = np.std(s_dist)
    
    return s_dist, A_t, B_t, len(s_dist)


def naive_estimate_bin(df, n_threshminus, n_threshplus):
    
    clone_fraction = 'Clone_fraction'
    clone_count = 'Clone_count'
    
    df_bis = df[(df['Clone_count_1'] > n_threshminus) & (df['Clone_count_1'] <= n_threshplus) ]
    df_bis['log_f_1'] = np.log(df_bis[clone_fraction + '_1'])
    df_bis['log_f_2'] = np.log(df_bis[clone_fraction + '_2'])
    
    
    s_dist = df_bis['log_f_2'] - df_bis['log_f_1']
    #f = df_bis[clone_fraction + '_1']
    s_dist = s_dist[s_dist != -np.inf]
    #l = np.argwhere(s_dist != -np.inf)
    #l = l.flatten()
    #f_bis = f.iloc[l]
    #f_bis = np.array(f_bis)
    s_dist = np.array(s_dist)

    N_obs = len(s_dist)
    
    log_fold_sum = np.sum(s_dist)
    A_naif = log_fold_sum/N_obs
    B_naif = (1/N_obs)*np.sum((s_dist - (A_naif))**2)

    return A_naif, B_naif



def naive_estimate_HR(filename_list, individual, path):
    
    
    #path = '../../Data-sets/Data_Harlan_Robins_Zuzia/'
    colnames = ['Clone fraction','Clone count', 'N. Seq. CDR3', 'AA. Seq. CDR3']
    
    for e1, file1 in tqdm(enumerate(filename_list)):
        for e2, file2 in enumerate(filename_list):
            if e1 < e2:
                temp1 = re.findall(r'\d+', file1)
                temp2 = re.findall(r'\d+', file2)
                res1 = list(map(int, temp1))
                res2 = list(map(int, temp2))
                duration = float(res2[1]-res1[1])

                if duration > 150:
                    number_clones, df = import_data(path, file1, file2, colnames,colnames)
                    print('Number of clones ' + str(number_clones) + ' ' + str(duration) + ' days')
                    #df['Clone_fraction_1'] =  df['Clone_fraction_1']/100
                    #df['Clone_fraction_2'] =  df['Clone_fraction_2']/100
                    deltaT = duration
                    optimization_binning(df, deltaT, individual)
                    
                else:
                    pass
            else:
                pass 

def optimization_binning_russian(df, delta_t, individual):
    
    df_copy = df.copy()
    #df_copy['Clone_fraction_1'] = df_copy['Clone_fraction_1']/100
    #df_copy['Clone_fraction_2'] = df_copy['Clone_fraction_2']/100
    df_bis = df_copy[(df_copy['Clone_fraction_1']> 1e-6) & (df_copy['Clone_fraction_2'] > 0)]
    df_bis['Id'] = 1
    N_count = df_bis.groupby(['Clone_count_1'])['Id'].sum()
    nmax = np.max(df_bis['Clone_count_1'])
    n_threshminus = []
    n_threshplus = []
    
    
    #pairs = np.array([[100,200], [200,400], [400,800], [800,2000], [2000, nmax]])
    pairs = np.array([[3,5], [5,7], [7,15], [15,100], [100, 1000], [1000, nmax]])
    
    #n_bins = np.floor(bins*NreadsI)
    #print(bins.astype(int))
    
    paras_db = np.zeros((len(pairs), 2))
    N_bin = np.zeros((len(pairs)))
    median = np.zeros((len(pairs)))
    median_fth = np.zeros((len(pairs)))
    
    for it, (n1,n2) in enumerate(pairs):
        #Bins are build the following way, empirical count is strictly larger than n_1 and smaller than n_2
        paras_db[it, :] = naive_estimate_bin(df_bis, n1, n2)
        paras_db[it, :] = paras_db[it, :]/delta_t
        N_bin[it] = np.sum(N_count[n1+1: n2])
        print(it, ' bin(s) has/(ve) been inferred')
        df_ter = df_bis[(df_bis['Clone_count_1']>n1) & (df_bis['Clone_count_1']<=n2)]
        median[it] = np.median(df_ter['Clone_count_1'])
        median_fth[it] = np.median(df_ter['Clone_fraction_1']) 
        
    d = {'n_thresh_minus': pairs[:,0], 'n_thresh_plus': pairs[:,1], 'A': paras_db[:,0], 'B': paras_db[:,1], '#clones_bin': N_bin, 'median': median, 'median_fth': median_fth}
    df_inference = pd.DataFrame(data =d )
    df_inference.to_csv(individual + '_binning_naive_log' , sep = '\t')

def optimization_binning(df, delta_t, individual):
    
    df_copy = df.copy()
    df_copy['Clone_fraction_1'] = df_copy['Clone_fraction_1']/100
    df_copy['Clone_fraction_2'] = df_copy['Clone_fraction_2']/100
    df_bis = df_copy[(df_copy['Clone_fraction_1']> 5e-6) & (df_copy['Clone_fraction_2'] > 0)]
    df_bis['Id'] = 1
    N_count = df_bis.groupby(['Clone_count_1'])['Id'].sum()
    nmax = np.max(df_bis['Clone_count_1'])
    n_threshminus = []
    n_threshplus = []
    
    #pairs = np.array([[2,5],[3,5], [5,7], [7,15], [15,100], [100, 1000], [1000, nmax]])
    
    #df_bis['bin_qcut'], bins = pd.qcut(df_bis['Clone_count_1'], q =10, duplicates = 'drop', retbins = True)
    
    pairs = np.array([[100,200], [200,400], [400,800], [800,2000], [2000, nmax]])
    #l = len(bins) 
    
    #n_bins = np.floor(bins*NreadsI)
    #print(bins.astype(int))
    
    paras_db = np.zeros((len(pairs), 2))
    N_bin = np.zeros((len(pairs)))
    median = np.zeros((len(pairs)))
    median_fth = np.zeros((len(pairs)))
    
    for it, (n1,n2) in enumerate(pairs):
        #Bins are build the following way, empirical count is strictly larger than n_1 and smaller than n_2
        paras_db[it, :] = naive_estimate_bin(df_bis, n1, n2)
        N_bin[it] = np.sum(N_count[n1+1: n2])
        print(it, ' bin(s) has/(ve) been inferred')
        df_ter = df_bis[(df_bis['Clone_count_1']>n1) & (df_bis['Clone_count_1']<=n2)]
        median[it] = np.median(df_ter['Clone_count_1'])
        median_fth[it] = np.median(df_ter['Clone_fraction_1']) 
        
    d = {'n_thresh_minus': pairs[:,0], 'n_thresh_plus': pairs[:,1], 'A': paras_db[:,0], 'B': paras_db[:,1], '#clones_bin': N_bin, 'median': median, 'median_fth': median_fth}
    df_inference = pd.DataFrame(data =d )
    df_inference.to_csv(individual + '_binning_naive_HR_log_' + str(round(delta_t, 4)), sep = '\t')

def optimization_binning_chu(df, delta_t, individual):
    
    df_copy = df.copy()
    #df_copy['Clone_fraction_1'] = df_copy['Clone_fraction_1']/100
    #df_copy['Clone_fraction_2'] = df_copy['Clone_fraction_2']/100
    df_bis = df_copy[(df_copy['Clone_fraction_1']> 5e-6) & (df_copy['Clone_fraction_2'] > 0)]
    df_bis['Id'] = 1
    N_count = df_bis.groupby(['Clone_count_1'])['Id'].sum()
    nmax = np.max(df_bis['Clone_count_1'])
    n_threshminus = []
    n_threshplus = []
    
    
    #pairs = np.array([[100,200], [200,400], [400,800], [800,2000], [2000, nmax]])
    pairs = np.array([[100,200], [200,400], [400,800], [800,2000], [2000, nmax]])
    
    #n_bins = np.floor(bins*NreadsI)
    #print(bins.astype(int))
    
    paras_db = np.zeros((len(pairs), 2))
    N_bin = np.zeros((len(pairs)))
    median = np.zeros((len(pairs)))
    median_fth = np.zeros((len(pairs)))
    
    for it, (n1,n2) in enumerate(pairs):
        #Bins are build the following way, empirical count is strictly larger than n_1 and smaller than n_2
        paras_db[it, :] = naive_estimate_bin(df_bis, n1, n2)
        paras_db[it, :] = paras_db[it, :]/delta_t
        N_bin[it] = np.sum(N_count[n1+1: n2])
        print(it, ' bin(s) has/(ve) been inferred')
        df_ter = df_bis[(df_bis['Clone_count_1']>n1) & (df_bis['Clone_count_1']<=n2)]
        median[it] = np.median(df_ter['Clone_count_1'])
        median_fth[it] = np.median(df_ter['Clone_fraction_1']) 
        
    d = {'n_thresh_minus': pairs[:,0], 'n_thresh_plus': pairs[:,1], 'A': paras_db[:,0], 'B': paras_db[:,1], '#clones_bin': N_bin, 'median': median, 'median_fth': median_fth}
    df_inference = pd.DataFrame(data =d )
    df_inference.to_csv(individual + '_binning_naive_log' , sep = '\t')

###example for P1:
path = '../../Data-sets/Data_Harlan_Robins_Zuzia/'
filename1_Patient1 = 'S1_0_F1_.txt'
filename2_Patient1 = 'S1_30_F1_.txt'
filename3_Patient1 = 'S1_57_F1_.txt'
filename4_Patient1 = 'S1_85_F1_.txt'
filename5_Patient1 = 'S1_156_F1_.txt'
filename6_Patient1 = 'S1_183_F1_.txt'
filename7_Patient1 = 'S1_212_F1_.txt'
filename8_Patient1 = 'S1_370_F1_.txt'
colnames = ['Clone fraction','Clone count', 'N. Seq. CDR3', 'AA. Seq. CDR3']

filename_list_P1 = [filename1_Patient1, filename2_Patient1, filename3_Patient1, filename4_Patient1, filename5_Patient1, filename6_Patient1, filename7_Patient1, filename8_Patient1]

#filename_list = filename_list_P1
#individual = 'P1'
#naive_estimate_HR(filename_list, individual)
