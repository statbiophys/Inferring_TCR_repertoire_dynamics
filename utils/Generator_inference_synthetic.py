import os
import sys
import math
import time
import numpy as np
import pandas as pd 
from generator_tools import *
from functools import partial
from inference_tools import *
import glob
from multiprocessing import Pool, cpu_count


def optimisation_learning_generations(n_thresh_scal, path, df_list, method, t):

    freq_dtype = float
    nfbins = 1200
    fmin = -10.0
    

    if method == 'poisson':
        acq_model_type = 3

    elif method == 'negative_binomial':
        acq_model_type = 2

    n_thresh_2 = 2
    f_thresh_scal = n_thresh_scal/1e6


    param_estimates = np.zeros((len(df_list), 4))

    for e, file in enumerate(df_list):

        A_index = file.index('A_') 
        print(file[A_index+2:A_index+6])
        A = float(file[A_index+2:A_index+6])
        B_index = file.index('B_')
        print(file[B_index+2:B_index+6])
        B = float(file[B_index+2:B_index+6])

        alpha = -1 + 2*A/B 

        if method == 'poisson':
            paras = [alpha, fmin]
        elif method == 'negative_binomial':
            paras = [alpha, 0.7, 1.1, fmin]


        Model = Inference_dynamics(path + file, paras, paras)
        df = Model.data_thresh_freq_dyn(f_thresh_scal, n_thresh_2)
        
        param_estimates[e, :] = Model.learning_diffusion_1D_two_cond(f_thresh_scal, n_thresh_2)
        
    d = {'n_thresh': n_thresh_scal*np.ones((len(df_list))) , 'A': param_estimates[:,0], 'B': param_estimates[:,1], 'Delta_A': param_estimates[:,2], 'Delta_B': param_estimates[:,3]}
    df_param = pd.DataFrame(data=d)

    print(df_param)


    dirName_para = 'Param_estimates_' + method + '_A_' +  str(A) + '_B_' + str(B) + '_' + str(t) + '_year(s)'     
    os.makedirs(dirName_para, exist_ok=True) 
    

    df_param.to_csv(dirName_para + '/' + 'n_thresh_' + str(n_thresh_scal)  , sep= '\t')


def main():

    np.seterr(divide = 'ignore') 
    np.warnings.filterwarnings('ignore')

    method = 'negative_binomial'
    NreadsI, NreadsII, t_ime = '1e6', '1e6', '2'


    # Synthetic data generation

    print('execution starting...')

    st = time.time()

    # A = -0.1, -0.2, -0.3, -0.4, -0.5, -0.
    A = -0.5
    B = 0.90
    N_0 = 40


    NreadsI = float(NreadsI)
    NreadsII = float(NreadsII)

    #t = float(t_ime)
    t = 2
    N = np.arange(10) #real job on zuzia
    #N = np.arange(2)  #test on local

    pool = Pool(cpu_count())
    synthetic_data_generator_LB_partial = partial(synthetic_data_generator_LB, A = A, B = B, NreadsI = NreadsI, NreadsII=NreadsII, t=t, N_0=N_0, method=method)
    pool.map(synthetic_data_generator_LB_partial, N)
    #print("TIME FOR SYNTHETIC DATA--- %s seconds ---" % (time.time() - st))



    # Inference strategy

    if NreadsI == NreadsII:
        key_sym = '_sym_'

    else:
        key_sym = '_asym_'

    name_B = str(B)
    if len(name_B) == 3:
        name_B = name_B + '0'
    else:
        pass 

    path = 'Synthetic_data_' + method + key_sym + 'A_' + str(A) + '_B_' + name_B + '_' + str(t) + '_year(s)/*'  
    print('repository name' + path)
    df_list = glob.glob(path)

    print(df_list)

    
    n_thresh_vec = np.array([1, 2, 5, 10, 20, 50, 100, 200, 300, 500]) 
    #n_thresh_vec = np.array([1])

    optimisation_nthresh_partial = partial(optimisation_learning_generations, path = '', df_list = df_list, method = method, t=t)
    pool = Pool(cpu_count())
    pool.map(optimisation_nthresh_partial, n_thresh_vec)

    


if __name__ == '__main__': 

    #method =sys.argv[1]
    #print(method)
    #NreadsI=sys.argv[2]
    #print(NreadsI)
    #NreadsII=sys.argv[3]
    #print(NreadsII)
    #t_ime=sys.argv[4]
    #print(t_ime)

    #main(method, NreadsI, NreadsII, t_ime)
    main()