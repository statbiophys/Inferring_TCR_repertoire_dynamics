import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import numdifftools as nd
from numpy import linalg as LA
from scipy.optimize import minimize
from scipy import stats
import math
import time
from copy import deepcopy 
from functools import partial
import numdifftools as nd 
from tqdm.notebook import trange, tqdm
from multiprocessing import Pool, cpu_count, Process
np.seterr(divide = 'ignore') 


#=========================================Parent-object-for-the-inference=====================================

class Inference_tools():

    def __init__(self, path, filename1 = None, filename2 = None, colnames1= None, colnames2 = None, HR = False):

        self.path = path
        if filename1 != None and filename2 != None and colnames1 != None and colnames2 != None:
            self.filename1 = filename1
            self.filename2 = filename2
            self.colnames1 = colnames1
            self.colnames2 = colnames2
            if HR:
                self.N_clones, self.df = self.frequencies_adjustment_HR()
            else:
                self.N_clones , self.df = self.import_data()

        else:
            self.df = pd.read_csv(path, sep = '\t')

        
        # For the integration along the frequencies variable 
        self.nfbins = 1200
        self.N_clones = len(self.df)
            
        self.counts = self.df.loc[:,['Clone_count_1', 'Clone_count_2']]
        self.sparse_rep = self.get_sparserep(self.df)
        self.indn1, self.indn2, self.sparse_rep_counts, self.unicountvals_1, self.unicountvals_2, self.NreadsI, self.NreadsII = self.sparse_rep
        self.N_data_ind = np.vstack((self.indn1, self.indn2)).T




    def frequencies_adjustment_HR(self):

        N_clones , df = self.import_data()
        df_bis = df.copy()
        df_bis['Clone_fraction_1'] = df['Clone_fraction_1']/100
        df_bis['Clone_fraction_2'] = df['Clone_fraction_2']/100

        return N_clones, df_bis



    def cumulative_distributions(self, Clone_fraction):

        """
        This method uses the data-frame and the Clone_fraction name
        output : 
        x_Patient : log(frequencies) sorted from the smallest to the largest values + unique values
        y_Patient : log(P(f>frequencies))
        slope_Patient : slope of the log-cumulative distribution obtained with a fit
        intercept_Patient : intercept of the former linear regression fit
        Y_patient : predicted data from the linear regression

        REECRIRE CETTE FONCTION
        """

        df = self.df
        
        df_true = df[df[Clone_fraction] != 0]
        df_true_new = df_true.copy()
        df_true_new['Id'] = 1
        distribution = df_true_new.groupby([Clone_fraction])['Id'].sum()
        frequencies = df_true_new[Clone_fraction].unique()
        frequencies = sorted(frequencies)
        
        x_Patient = np.log(frequencies)
        y_Patient = np.log(np.max(np.cumsum(distribution)) - np.cumsum(distribution))
        
        slope_Patient, intercept_Patient, r_value_Patient, p_value_Patient, std_err_Patient = stats.linregress(x_Patient[:-1], y_Patient[:-1,])
        Y_Patient = intercept_Patient + slope_Patien1*np.array(x_Patient)
        
        
        return x_Patient, y_Patient, slope_Patient, intercept_Patient, Y_Patient

    def import_data(self):
        '''
        Reads indata from two datasets and merges based on nt sequence
        Outputs dataframe of pair counts for all clones.
        '''
        #print('OK')
        path = self.path
        filename1 = self.filename1
        filename2 = self.filename2
        colnames1 = self.colnames1
        colnames2 = self.colnames2 

        mincount = 0
        maxcount = np.inf
        headerline=0 #line number of headerline
        newnames=['Clone_fraction','Clone_count','ntCDR3','AACDR3']    
        with open(path+filename1, 'r') as f:
            F1Frame_chunk=pd.read_csv(f,delimiter='\t',usecols=colnames1,header=headerline)[colnames1]
        with open(path+filename2, 'r') as f:
            F2Frame_chunk=pd.read_csv(f,delimiter='\t',usecols=colnames2,header=headerline)[colnames2]
        F1Frame_chunk.columns=newnames
        F2Frame_chunk.columns=newnames
        suffixes=('_1','_2')
        mergedFrame=pd.merge(F1Frame_chunk,F2Frame_chunk,on=newnames[2],suffixes=suffixes,how='outer')
        for nameit in [0,1]:
            for labelit in suffixes:
                mergedFrame.loc[:,newnames[nameit]+labelit].fillna(int(0),inplace=True)
                if nameit==1:
                    mergedFrame.loc[:,newnames[nameit]+labelit].astype(int)
        def dummy(x):
            val=x[0]
            if pd.isnull(val):
                val=x[1]    
            return val
        mergedFrame.loc[:,newnames[3]+suffixes[0]]=mergedFrame.loc[:,[newnames[3]+suffixes[0],newnames[3]+suffixes[1]]].apply(dummy,axis=1) #assigns AA sequence to clones, creates duplicates
        mergedFrame.drop(newnames[3]+suffixes[1], 1,inplace=True) #removes duplicates
        mergedFrame.rename(columns = {newnames[3]+suffixes[0]:newnames[3]}, inplace = True)
        mergedFrame=mergedFrame[[newname+suffix for newname in newnames[:2] for suffix in suffixes]+[newnames[2],newnames[3]]]
        filterout=((mergedFrame.Clone_count_1<mincount) & (mergedFrame.Clone_count_2==0)) | ((mergedFrame.Clone_count_2<mincount) & (mergedFrame.Clone_count_1==0)) #has effect only if mincount>0
        number_clones=len(mergedFrame)
        return number_clones,mergedFrame.loc[((mergedFrame.Clone_count_1<=maxcount) & (mergedFrame.Clone_count_2<=maxcount)) & ~filterout]

    def get_sparserep(self, df):


        '''
        Tranforms {(n1,n2)} data stored in pandas dataframe to a sparse 1D representation.
        unicountvals_1(2) are the unique values of n1(2).
        sparse_rep_counts gives the counts of unique pairs.
        indn1(2) is the index of unicountvals_1(2) giving the value of n1(2) in that unique pair.
        len(indn1)=len(indn2)=len(sparse_rep_counts)
        '''

        counts = df.loc[:,['Clone_count_1', 'Clone_count_2']]

        counts['paircount'] = 1  # gives a weight of 1 to each observed clone

        clone_counts = counts.groupby(['Clone_count_1', 'Clone_count_2']).sum()
        sparse_rep_counts = np.asarray(clone_counts.values.flatten(), dtype=int)
        clonecountpair_vals = clone_counts.index.values
        indn1 = np.asarray([clonecountpair_vals[it][0] for it in range(len(sparse_rep_counts))], dtype=int)
        indn2 = np.asarray([clonecountpair_vals[it][1] for it in range(len(sparse_rep_counts))], dtype=int)
        NreadsI = np.sum(counts['Clone_count_1'])
        NreadsII = np.sum(counts['Clone_count_2'])

        unicountvals_1, indn1 = np.unique(indn1, return_inverse=True)
        unicountvals_2, indn2 = np.unique(indn2, return_inverse=True)

        return indn1, indn2, sparse_rep_counts, unicountvals_1, unicountvals_2, NreadsI, NreadsII
    
    def data_thresh(self, n_thresh):
        
        if n_thresh != 0:
            df_thresh= self.df[(self.df['Clone_count_1'] > n_thresh-1) ]
        
        return df_thresh

    def data_thresh_freq_noise(self, f_thresh):
        
        if f_thresh != 0:
            df_thresh= self.df[(self.df['Clone_fraction_1'] >= f_thresh) & (self.df['Clone_fraction_2'] >= f_thresh)]
        
        return df_thresh

    def data_thresh_freq_dyn_all(self, f_thresh):
        
        if f_thresh != 0:
            df_thresh= self.df[(self.df['Clone_fraction_1'] >= f_thresh)]
        return df_thresh

    def data_thresh_freq_dyn(self, f_thresh, n_thresh_2):
        
        if f_thresh != 0:
            df_thresh= self.df[(self.df['Clone_fraction_1'] >= f_thresh) & (self.df['Clone_count_2'] >= n_thresh_2)]
        return df_thresh

    def data_thresh_db(self, n_threshminus, n_threshplus):

        """
        Build the double bouded dataframe such as n_threhsminus <= n_1 <= n_threshplus
        """
        
        df_thresh= self.df[(self.df['Clone_count_1'] > n_threshminus-1) & (self.df['Clone_count_1'] < n_threshplus+1) ]
        
        return df_thresh

    def data_thresh_db_cond_n2(self, n_threshminus, n_threshplus, n2):

        """
        Build the double bouded dataframe such as n_threhsminus < n_1 <= n_threshplus
        """
        
        df_thresh= self.df[(self.df['Clone_count_1'] > n_threshminus) & (self.df['Clone_count_1'] <= n_threshplus) & (self.df['Clone_count_2']>=n2)]
        
        return df_thresh

    def data_thresh_db_cond_fth_n2(self, f_threshminus, f_threshplus, n2):

        """
        Build the double bouded dataframe such as n_threhsminus <= n_1 <= n_threshplus
        """
        
        df_thresh= self.df[(self.df['Clone_fraction_1'] >= f_threshminus) & (self.df['Clone_fraction_1'] < f_threshplus) & (self.df['Clone_count_2']>=n2)]
        
        return df_thresh
        

    
    @staticmethod
    def NegBinParMtr(m,v,nvec): #speed up only insofar as the log and exp are called once on array instead of multiple times on rows
        ''' 
        computes NegBin probabilities over the ordered (but possibly discontiguous) vector (nvec) 
        for mean/variance combinations given by the mean (m) and variance (v) vectors. 
        Note that m<v for negative binomial.
        Output is (len(m),len(nvec)) array
        '''
        nmax = np.max(nvec)
        p = 1-m/v
        r = m*m/v/p
        NBvec=np.arange(nmax+1,dtype=float)
        NBvec=np.log((NBvec+r[:,np.newaxis]-1)*(p[:,np.newaxis]/NBvec))
        NBvec[:,0]=r*np.log(m/v) #handle NBvec[0]=0, treated specially when m[0]=0, see below
        NBvec= np.exp(np.cumsum(NBvec,axis=1)) #save a bit here
        if m[0]==0:
            NBvec[0,:]= 0
            NBvec[0,0]= 1
        NBvec=NBvec[:,nvec]
        return NBvec
    
    @staticmethod
    def PoisPar(Mvec,unicountvals):
        #assert Mvec[0]==0, "first element needs to be zero"
        #nmax=unicountvals[-1]
        nmax = np.max(unicountvals)
        #print(nmax)
        nlen=len(unicountvals)
        mlen=len(Mvec)
        Nvec=unicountvals
        logNvec=-np.insert(np.cumsum(np.log(np.arange(1,nmax+1))),0,0.)[unicountvals] #avoid n=0 nans  
        Nmtr=np.exp(Nvec[np.newaxis,:]*np.log(Mvec)[:,np.newaxis]+logNvec[np.newaxis,:]-Mvec[:,np.newaxis]) # np.log(Mvec) throws warning: since log(0)=-inf
        if Mvec[0]==0:
            Nmtr[0,:]=np.zeros((nlen,)) #when m=0, n=0, and so get rid of nans from log(0)
            Nmtr[0,0]=1. #handled belowacq_model_type
        if unicountvals[0]==0: #if n=0 included get rid of nans from log(0)
            Nmtr[:,0]=np.exp(-Mvec)
        return Nmtr


class Inference_dynamics(Inference_tools):

    def __init__(self, path, paras_1, paras_2, filename1 = None, filename2 = None, colnames1= None, colnames2 = None, t = None, Pers = False):

        super().__init__(path, filename1 , filename2 , colnames1, colnames2 )
        self.paras_1 = paras_1
        self.paras_2 = paras_2
        self.Pers = Pers

        if len(paras_1) == 4:
            self.method = 'negative_binomial'

        elif len(paras_1) == 2:
            self.method = 'poisson'
        self.freq_dtype = float
        self.logrhofvec, self.logfvec = self.get_rhof()
        
        if t != None:
            self.t  = t
        self.svec, self.logfvecwide, self.f2s_step = self.s_vec()
        self.logfvecsecond = self.logfvecwide

        
        
    def s_vec(self):

        """
        s_vec is the log-fold change vector
        """

        #Definition of svec
        logfvec = self.logfvec
        if self.Pers:
            smax = 30.0     #maximum absolute logfold change value

        else:
            smax = 25.0
        s_step = 0.1
        s_step_old= s_step
        logf_step= logfvec[1] - logfvec[0] #use natural log here since f2 increments in increments in exp(). 
        f2s_step= int(round(s_step/logf_step)) #rounded number of f-steps in one s-step
        s_step= float(f2s_step)*logf_step
        smax= s_step*(smax/s_step_old)
        svec= s_step*np.arange(0,int(round(smax/s_step)+1))   
        svec= np.append(-svec[1:][::-1],svec)
        smaxind=(len(svec)-1)/2
        f2s_step=int(round(s_step/logf_step)) #rounded number of f-steps in one s-step
        logfmin=logfvec[0 ]-f2s_step*smaxind*logf_step
        logfmax=logfvec[-1]+f2s_step*smaxind*logf_step
        logfvecwide=np.linspace(logfmin,logfmax,int(len(logfvec)+2*smaxind*f2s_step)) #a wider domain for the second frequency f2=f1*exp(s)

        return svec, logfvecwide, f2s_step

    def get_rhof(self):
        '''
        generates power law (power is alpha_rho) clone frequency distribution over 
        freq_nbins discrete logarithmically spaced frequences between fmin and 1 of dtype freq_dtype
        Outputs log probabilities obtained at log frequencies'''
        
        freq_dtype = self.freq_dtype
        freq_nbins = self.nfbins
        paras_1 = self.paras_1
        alpha_rho = paras_1[0]
        fmin = np.power(10, self.paras_1[-1])
        fmax=1e0
        
        logfvec=np.linspace(np.log10(fmin),np.log10(fmax),freq_nbins)
        logfvec=np.array(np.log(np.power(10,logfvec)) ,dtype=freq_dtype).flatten()  
        logrhovec=logfvec*alpha_rho
        integ=np.exp(logrhovec+logfvec,dtype=freq_dtype)
        normconst=np.log(np.dot(np.diff(logfvec)/2.,integ[1:]+integ[:-1]))
        logrhovec-=normconst 
        return logrhovec, logfvec
            

    def get_logPn_f_db(self, n_threshminus, n_threshplus, day):

        """
        Build the double bouded logPnf such as n_threhsminus <= n_1 <= n_threshplus
        
        """
    
        method = self.method
        paras_1 = self.paras_1
        paras_2 = self.paras_2
        
        df_thresh= super().data_thresh_db(n_threshminus, n_threshplus)
        sparse_rep_thresh = super().get_sparserep(df_thresh)
        
        indn1_thresh, indn2_thresh, sparse_rep_counts_thresh, unicountvals_1_thresh, unicountvals_2_thresh, NreadsI_thresh, NreadsII_thresh = sparse_rep_thresh
        
        NreadsI = self.NreadsI
        NreadsII = self.NreadsII
        
        logrhofvec, logfvec = self.logrhofvec, self.logfvec
        logfvecsecond = self.logfvecsecond
            
        if day == 'initial':
            unicounts, Nreads = unicountvals_1_thresh, NreadsI
            logfvec_tmp=deepcopy(logfvec)
            paras = paras_1
            
        elif day == 'final': 
            unicounts, Nreads = unicountvals_2_thresh, NreadsII
            logfvec_tmp=deepcopy(logfvecsecond)
            paras = paras_2

            
        Pn_f=np.zeros((len(logfvec_tmp),len(unicounts)))
        
        if method == 'poisson':
            mean_n=Nreads*np.exp(logfvec_tmp)
            Pn_f= super().PoisPar(mean_n,unicounts)

        if method == 'negative_binomial':
            beta_mv = paras[1]
            alpha_mv =paras[2]
            mean_n=Nreads*np.exp(logfvec_tmp)
            var_n=mean_n+beta_mv*np.power(mean_n,alpha_mv)
            Pn_f= super().NegBinParMtr(mean_n,var_n,unicounts)
            
        return np.log(Pn_f)

    def get_logPn_f_db_n2_cond(self, n_threshminus, n_threshplus, n2, day):

        """
        Build the double bouded logPnf such as n_threhsminus <= n_1 <= n_threshplus
        
        """
    
        method = self.method
        paras_1 = self.paras_1
        paras_2 = self.paras_2
        
        df_thresh= super().data_thresh_db_cond_n2(n_threshminus, n_threshplus, n2)
        sparse_rep_thresh = super().get_sparserep(df_thresh)
        
        indn1_thresh, indn2_thresh, sparse_rep_counts_thresh, unicountvals_1_thresh, unicountvals_2_thresh, NreadsI_thresh, NreadsII_thresh = sparse_rep_thresh
        
        NreadsI = self.NreadsI
        NreadsII = self.NreadsII
        
        logrhofvec, logfvec = self.logrhofvec, self.logfvec
        logfvecsecond = self.logfvecsecond
            
        if day == 'initial':
            unicounts, Nreads = unicountvals_1_thresh, NreadsI
            logfvec_tmp=deepcopy(logfvec)
            paras = paras_1
            
        elif day == 'final': 
            unicounts, Nreads = unicountvals_2_thresh, NreadsII
            logfvec_tmp=deepcopy(logfvecsecond)
            paras = paras_2

            
        Pn_f=np.zeros((len(logfvec_tmp),len(unicounts)))
        
        if method == 'poisson':
            mean_n=Nreads*np.exp(logfvec_tmp)
            Pn_f= super().PoisPar(mean_n,unicounts)

        if method == 'negative_binomial':
            beta_mv = paras[1]
            alpha_mv =paras[2]
            mean_n=Nreads*np.exp(logfvec_tmp)
            var_n=mean_n+beta_mv*np.power(mean_n,alpha_mv)
            Pn_f= super().NegBinParMtr(mean_n,var_n,unicounts)
            
        return np.log(Pn_f)

    def get_logPn_f_condition_freq(self, f_thresh, replicate, unicountvals1_cond, unicountvals2_cond):


        method = self.method
        paras_1 = self.paras_1
        paras_2 = self.paras_2
        
        NreadsI = self.NreadsI
        NreadsII = self.NreadsII
        
        logrhofvec, logfvec = self.logrhofvec, self.logfvec
        logfvecsecond = self.logfvecsecond
            
        if replicate == 'first':
            unicounts, Nreads = unicountvals1_cond.astype(int), NreadsI
            logfvec_tmp=deepcopy(logfvec)
            paras = paras_1
            
        elif replicate == 'second': 
            unicounts, Nreads = unicountvals2_cond.astype(int), NreadsII
            logfvec_tmp=deepcopy(logfvecsecond)
            paras = paras_2

        Pn_f = np.zeros((len(logfvec_tmp),len(unicounts)))
        
        if method == 'poisson':
            mean_n=Nreads*np.exp(logfvec_tmp)
            Pn_f= super().PoisPar(mean_n,unicounts)

        if method == 'negative_binomial':
            beta_mv = paras[1]
            alpha_mv =paras[2]
            mean_n=Nreads*np.exp(logfvec_tmp)
            var_n=mean_n+beta_mv*np.power(mean_n,alpha_mv)
            Pn_f= super().NegBinParMtr(mean_n,var_n,unicounts)
            
        return np.log(Pn_f)


    def get_logPn_f(self, n_thresh, day):
    
        method = self.method
        paras_1 = self.paras_1
        paras_2 = self.paras_2
        
        df_thresh= super().data_thresh(n_thresh)
        sparse_rep_thresh = super().get_sparserep(df_thresh)
        
        indn1_thresh, indn2_thresh, sparse_rep_counts_thresh, unicountvals_1_thresh, unicountvals_2_thresh, NreadsI_thresh, NreadsII_thresh = sparse_rep_thresh
        
        NreadsI = self.NreadsI
        NreadsII = self.NreadsII
        
        logrhofvec, logfvec = self.logrhofvec, self.logfvec
        logfvecsecond = self.logfvecsecond
            
        if day == 'initial':
            unicounts, Nreads = unicountvals_1_thresh, NreadsI
            logfvec_tmp=deepcopy(logfvec)
            paras = paras_1
            
        elif day == 'final': 
            unicounts, Nreads = unicountvals_2_thresh, NreadsII
            logfvec_tmp=deepcopy(logfvecsecond)
            paras = paras_2

            
        Pn_f=np.zeros((len(logfvec_tmp),len(unicounts)))
        
        if method == 'poisson':
            mean_n=Nreads*np.exp(logfvec_tmp)
            Pn_f= super().PoisPar(mean_n,unicounts)

        if method == 'negative_binomial':
            beta_mv = paras[1]
            alpha_mv =paras[2]
            mean_n=Nreads*np.exp(logfvec_tmp)
            var_n=mean_n+beta_mv*np.power(mean_n,alpha_mv)
            Pn_f= super().NegBinParMtr(mean_n,var_n,unicounts)
            
        return np.log(Pn_f)

    def get_logPn_f_new(self, f_thresh, n_thresh_2, day):
    
        method = self.method
        paras_1 = self.paras_1
        paras_2 = self.paras_2
        
        df_thresh= super().data_thresh_freq_dyn(f_thresh, n_thresh_2)
        sparse_rep_thresh = super().get_sparserep(df_thresh)
        
        indn1_thresh, indn2_thresh, sparse_rep_counts_thresh, unicountvals_1_thresh, unicountvals_2_thresh, NreadsI_thresh, NreadsII_thresh = sparse_rep_thresh
        
        NreadsI = self.NreadsI
        NreadsII = self.NreadsII
        
        logrhofvec, logfvec = self.logrhofvec, self.logfvec
        logfvecsecond = self.logfvecsecond
            
        if day == 'initial':
            unicounts, Nreads = unicountvals_1_thresh.astype(int), NreadsI
            logfvec_tmp=deepcopy(logfvec)
            paras = paras_1
            
        elif day == 'final': 
            unicounts, Nreads = unicountvals_2_thresh.astype(int), NreadsII
            logfvec_tmp=deepcopy(logfvecsecond)
            paras = paras_2

            
        Pn_f=np.zeros((len(logfvec_tmp),len(unicounts)))
        
        if method == 'poisson':
            mean_n=Nreads*np.exp(logfvec_tmp)
            Pn_f= super().PoisPar(mean_n,unicounts)

        if method == 'negative_binomial':
            beta_mv = paras[1]
            alpha_mv =paras[2]
            mean_n=Nreads*np.exp(logfvec_tmp)
            var_n=mean_n+beta_mv*np.power(mean_n,alpha_mv)
            Pn_f= super().NegBinParMtr(mean_n,var_n,unicounts)
            
        return np.log(Pn_f)

    def get_logPn_f_freq_all(self, f_thresh, day):
    
        method = self.method
        paras_1 = self.paras_1
        paras_2 = self.paras_2
        
        df_thresh= super().data_thresh_freq_dyn_all(f_thresh)
        sparse_rep_thresh = super().get_sparserep(df_thresh)
        
        indn1_thresh, indn2_thresh, sparse_rep_counts_thresh, unicountvals_1_thresh, unicountvals_2_thresh, NreadsI_thresh, NreadsII_thresh = sparse_rep_thresh
        
        NreadsI = self.NreadsI
        NreadsII = self.NreadsII
        
        logrhofvec, logfvec = self.logrhofvec, self.logfvec
        logfvecsecond = self.logfvecsecond
            
        if day == 'initial':
            unicounts, Nreads = unicountvals_1_thresh.astype(int), NreadsI
            logfvec_tmp=deepcopy(logfvec)
            paras = paras_1
            
        elif day == 'final': 
            unicounts, Nreads = unicountvals_2_thresh.astype(int), NreadsII
            logfvec_tmp=deepcopy(logfvecsecond)
            paras = paras_2

            
        Pn_f=np.zeros((len(logfvec_tmp),len(unicounts)))
        
        if method == 'poisson':
            mean_n=Nreads*np.exp(logfvec_tmp)
            Pn_f= super().PoisPar(mean_n,unicounts)

        if method == 'negative_binomial':
            beta_mv = paras[1]
            alpha_mv =paras[2]
            mean_n=Nreads*np.exp(logfvec_tmp)
            var_n=mean_n+beta_mv*np.power(mean_n,alpha_mv)
            Pn_f= super().NegBinParMtr(mean_n,var_n,unicounts)
            
        return np.log(Pn_f)

        
    def get_logPn_f_condition(self, f_thresh, day, n_lim = None):
        
        paras_1 = self.paras_1
        paras_2 = self.paras_2
        method = self.method
        logrhofvec, logfvec = self.logrhofvec, self.logfvec
        
        NreadsI = self.NreadsI
        NreadsII = self.NreadsII
        logfvecsecond = self.logfvecsecond

        n_thresh_1 = np.floor(f_thresh*NreadsI)
        n_thresh_2 = 1
        
        #if n_thresh != 0:
            #unicountvals1_cond = np.arange(n_thresh)
            
            #if method == 'negative_binomial':
            #    unicountvals2_cond = np.arange(n_lim)
                
            #elif method == 'poisson':
                #unicountvals2_cond = np.arange(n_lim*10)

        unicountvals1_cond = np.arange(n_thresh_1)
        unicountvals2_cond = np.arange(2)

            
        if day == 'initial':
            unicounts, Nreads = unicountvals1_cond.astype(int), NreadsI
            logfvec_tmp=deepcopy(logfvec)
            paras = paras_1
            
        elif day == 'final': 
            unicounts, Nreads = unicountvals2_cond.astype(int), NreadsII
            logfvec_tmp=deepcopy(logfvecsecond)
            paras = paras_2

        
        Pn_f = np.zeros((len(logfvec_tmp),len(unicounts)))
        
        if method == 'poisson':
            mean_n=Nreads*np.exp(logfvec_tmp)
            Pn_f= super().PoisPar(mean_n,unicounts)

        if method == 'negative_binomial':
            beta_mv = paras[1]
            alpha_mv =paras[2]
            mean_n=Nreads*np.exp(logfvec_tmp)
            var_n=mean_n+beta_mv*np.power(mean_n,alpha_mv)
            Pn_f= super().NegBinParMtr(mean_n,var_n,unicounts)
            
        return np.log(Pn_f)


    def get_logPn_f_condition_db(self, day, unicountvals1_cond, unicountvals2_cond):

        """db is for double bouding 
        the empirical abundance of the clone at initial time is bounded between n_1thresh- and n_1thresh+
        n_1thresh- <= n_1 <= n_1thresh+"""
        
        paras_1 = self.paras_1
        paras_2 = self.paras_2
        method = self.method
        logrhofvec, logfvec = self.logrhofvec, self.logfvec
        
        NreadsI = self.NreadsI
        NreadsII = self.NreadsII
        logfvecsecond = self.logfvecsecond
            
        if day == 'initial':
            unicounts, Nreads = unicountvals1_cond, NreadsI
            logfvec_tmp=deepcopy(logfvec)
            paras = paras_1
            
        elif day == 'final': 
            unicounts, Nreads = unicountvals2_cond, NreadsII
            logfvec_tmp=deepcopy(logfvecsecond)
            paras = paras_2

        Pn_f = np.zeros((len(logfvec_tmp),len(unicounts)))
        
        if method == 'poisson':
            mean_n=Nreads*np.exp(logfvec_tmp)
            Pn_f= super().PoisPar(mean_n,unicounts)

        if method == 'negative_binomial':
            beta_mv = paras[1]
            alpha_mv =paras[2]
            mean_n=Nreads*np.exp(logfvec_tmp)
            var_n=mean_n+beta_mv*np.power(mean_n,alpha_mv)
            Pn_f= super().NegBinParMtr(mean_n,var_n,unicounts)
            
        return np.log(Pn_f)

    def get_Ps(self, A, B):
            svec = self.svec
    
            quad = - 1/(2*B) * np.power((svec - A),2)
            return 1/(math.sqrt(2*(math.pi)*B))*np.exp(quad)

    def get_Ps_for_Post(self, A, B, svec):
            svec = self.svec
    
            quad = - 1/(2*B) * np.power((svec - A),2)
            return 1/(math.sqrt(2*(math.pi)*B))*np.exp(quad)

    
    
    def Plot_surface_diffusion(self, n_thresh):
        
        npoints = 20 #to be chosen by
        Avec = np.linspace(-3, 1, npoints)
        Bvec = np.logspace(-2, 2, npoints)

        LSurface =np.zeros((len(Bvec),len(Avec)))
        for i in range(len(Bvec)):
            for j in range(len(Avec)):
                LSurface[i, j]=  - self.get_log_likelihood_1D([Avec[j], Bvec[i]], n_thresh, n_lim)

        Amesh, Bmesh = np.meshgrid(Avec,Bvec)
        a,b = np.where(LSurface == np.max(LSurface))
        optA = Amesh[a[0],b[0]]
        optB = Bmesh[a[0],b[0]]
        
        return LSurface, Avec, Bvec


    
    
    def _callbackFdiffexpr(self, Xi, n_thresh, n_lim, verbose, A_0, B_0): #case dependent
        '''prints iteration info. called scipy.minimize'''

        global curr_iter 

        if len(Xi) == 2:
            print('{0: d}   {1: 3.6f} {2: 3.6f} {3: 3.6f}'.format(curr_iter, Xi[0], Xi[1], self.get_log_likelihood_dynamic(Xi, n_thresh, n_lim, verbose, A_0, B_0)) +'\n') 
        elif len(Xi) == 4:
            print('{0: d}   {1: 3.6f} {2: 3.6f}   {3: 3.6f} {4: 3.6f} {5: 3.6f}  '.format(curr_iter, Xi[0], Xi[1],Xi[2], Xi[3], self.get_log_likelihood_dynamic(Xi, n_thresh, n_lim, verbose, A_0, B_0)) +'\n') 
        curr_iter += 1

    def _callbackFdiffexpr_simple(self, Xi): #case dependent
        '''prints iteration info. called scipy.minimize'''

        global curr_iter 

        if len(Xi) == 2:
            print('{0: d}   {1: 3.6f} {2: 3.6f} '.format(curr_iter, Xi[0], Xi[1]) +'\n') 
        elif len(Xi) == 4:
            print('{0: d}   {1: 3.6f} {2: 3.6f}   {3: 3.6f} {4: 3.6f} '.format(curr_iter, Xi[0], Xi[1],Xi[2], Xi[3]) +'\n') 
        curr_iter += 1

    def _callbackFdiffexpr_1D(self, Xi, n_thresh, n_lim): #case dependent
        '''prints iteration info. called scipy.minimize'''

        global curr_iter 

        print('{0: d}   {1: 3.6f} {2: 3.6f} {3: 3.6f}'.format(curr_iter, Xi[0], Xi[1], self.get_log_likelihood_1D(Xi, n_thresh, n_lim)) +'\n') 
        curr_iter += 1

    def _callbackFdiffexpr_1D_simple(self, Xi): #case dependent
        '''prints iteration info. called scipy.minimize'''

        global curr_iter 

        print('{0: d}   {1: 3.6f} {2: 3.6f}'.format(curr_iter, Xi[0], Xi[1]) +'\n') 
        curr_iter += 1
        
    
    def get_log_likelihood_1D(self, PARAS, n_thresh, n_lim):

        paras_1 = self.paras_1
        paras_2 = self.paras_2
        freq_dtype = self.freq_dtype
        nfbins = self.nfbins
        method = self.method


        df_thresh= super().data_thresh(n_thresh)
        sparse_rep_thresh = super().get_sparserep(df_thresh)
        indn1_thresh, indn2_thresh, sparse_rep_counts_thresh, unicountvals_1_thresh, unicountvals_2_thresh, NreadsI_thresh, NreadsII_thresh = sparse_rep_thresh


        logrhofvec = self.logrhofvec
        logfvec = self.logfvec
        logfvecsecond = self.logfvecsecond
        svec = self.svec
        f2s_step = self.f2s_step

        unicountvals1_cond = np.arange(n_thresh)

        if method == 'negative_binomial':
            unicountvals2_cond = np.arange(n_lim)

        elif method == 'poisson':
            unicountvals2_cond = np.arange(n_lim*10)

        logPn1_f_cond = self.get_logPn_f_condition(n_thresh, 'initial', n_lim)
        logPn2_f_cond = self.get_logPn_f_condition(n_thresh, 'final', n_lim)

        dlogfby2=np.diff(logfvec)/2

        Pn1n2_s_cond1=np.zeros((len(svec),len(unicountvals1_cond), len(unicountvals2_cond))) 
        #for k in tqdm(range(len(svec))):
        for k in range(len(svec)):
            for i in range (len(unicountvals1_cond)):
                for j in range (len(unicountvals2_cond)):
                    integ = np.exp(logrhofvec+logPn2_f_cond[f2s_step*k:(f2s_step*k+len(logfvec)),j]+logPn1_f_cond[:,i]+ logfvec)
                    Pn1n2_s_cond1[k, i,j] = np.dot(dlogfby2,integ[1:] + integ[:-1])

        logPn1_f = self.get_logPn_f(n_thresh, 'initial')
        logPn2_f = self.get_logPn_f(n_thresh, 'final')

        Pn1n2_s=np.zeros((len(svec),len(sparse_rep_counts_thresh))) 
        for s_it,s in enumerate(svec):
            for it,(n1_it, n2_it) in enumerate(zip(indn1_thresh,indn2_thresh)):
                integ = np.exp(logrhofvec+logPn2_f[f2s_step*s_it:(f2s_step*s_it+len(logfvec)),n2_it]+logPn1_f[:,n1_it]+ logfvec)
                Pn1n2_s[s_it, it] = np.dot(dlogfby2,integ[1:] + integ[:-1])


        ds = np.diff(svec)/2
        N_obs = np.sum(sparse_rep_counts_thresh)

        A = PARAS[0]
        B = PARAS[1]

        integ = np.zeros(len(sparse_rep_counts_thresh))
        Ps = self.get_Ps(A, B)
        for it, (n1_it, n2_it) in enumerate(zip(indn1_thresh, indn2_thresh)):
            Pn1n2_ps = Pn1n2_s[:,it]*Ps
            integ[it] = np.dot(ds, Pn1n2_ps[1:] + Pn1n2_ps[:-1])
            
        cond1 = np.zeros((len(unicountvals1_cond), len(unicountvals2_cond)))
        for i in range (len(unicountvals1_cond)):
            for j in range (len(unicountvals2_cond)):
                Pn1n2_ps = Pn1n2_s_cond1[:,i,j]*Ps
                cond1[i,j] = np.dot(ds, Pn1n2_ps[1:] + Pn1n2_ps[:-1])            
        C = 1 - np.sum(cond1) 
    
        Energy_minus= np.dot(sparse_rep_counts_thresh, np.log(integ)) - N_obs*np.log(C)

        return Energy_minus


    def learning_diffusion_1D(self, f_thresh):

        paras_1 = self.paras_1
        paras_2 = self.paras_2
        freq_dtype = self.freq_dtype
        nfbins = self.nfbins
        method = self.method
        n_lim = 500

        df_thresh= super().data_thresh_freq_dyn_all(f_thresh)
        sparse_rep_thresh = super().get_sparserep(df_thresh)
        indn1_thresh, indn2_thresh, sparse_rep_counts_thresh, unicountvals_1_thresh, unicountvals_2_thresh, NreadsI_thresh, NreadsII_thresh = sparse_rep_thresh

        NreadsI = self.NreadsI
        NreadsII = self.NreadsII

        n_thresh_1 = np.floor(f_thresh*NreadsI)
        #not sure if I should take 
        #n_thresh_2 = 1

        logrhofvec = self.logrhofvec
        logfvec = self.logfvec
        logfvecsecond = self.logfvecsecond
        svec = self.svec
        f2s_step = self.f2s_step

        unicountvals1_cond = np.arange(n_thresh_1).astype(int)
        unicountvals2_cond = np.arange(n_lim)

        #unicountvals1_cond = np.arange(n_thresh)

        #if method == 'negative_binomial':
        #    unicountvals2_cond = np.arange(n_lim)

        #elif method == 'poisson':
        #    unicountvals2_cond = np.arange(n_lim*10)

        
        logPn1_f_cond = self.get_logPn_f_condition_freq(f_thresh, 'first', unicountvals1_cond, unicountvals2_cond)
        logPn2_f_cond = self.get_logPn_f_condition_freq(f_thresh, 'second', unicountvals1_cond, unicountvals2_cond)

        dlogfby2=np.diff(logfvec)/2

        Pn1n2_s_cond1=np.zeros((len(svec),len(unicountvals1_cond), len(unicountvals2_cond))) 
        for k in tqdm(range(len(svec))):
            for i in range (len(unicountvals1_cond)):
                for j in range (len(unicountvals2_cond)):
                    integ = np.exp(logrhofvec+logPn2_f_cond[f2s_step*k:(f2s_step*k+len(logfvec)),j]+logPn1_f_cond[:,i]+ logfvec)
                    Pn1n2_s_cond1[k, i,j] = np.dot(dlogfby2,integ[1:] + integ[:-1])

        logPn1_f = self.get_logPn_f_freq_all(f_thresh, 'initial')
        logPn2_f = self.get_logPn_f_freq_all(f_thresh, 'final')

        Pn1n2_s=np.zeros((len(svec),len(sparse_rep_counts_thresh))) 
        for s_it,s in enumerate(svec):
            for it,(n1_it, n2_it) in enumerate(zip(indn1_thresh,indn2_thresh)):
                integ = np.exp(logrhofvec+logPn2_f[f2s_step*s_it:(f2s_step*s_it+len(logfvec)),n2_it]+logPn1_f[:,n1_it]+ logfvec)
                Pn1n2_s[s_it, it] = np.dot(dlogfby2,integ[1:] + integ[:-1])


        ds = np.diff(svec)/2
        N_obs = np.sum(sparse_rep_counts_thresh)


        def cost(PARAS):

            A = PARAS[0]
            B = PARAS[1]


            integ = np.zeros(len(sparse_rep_counts_thresh))
            Ps = self.get_Ps(A, B)
            for it, (n1_it, n2_it) in enumerate(zip(indn1_thresh, indn2_thresh)):
                Pn1n2_ps = Pn1n2_s[:,it]*Ps
                integ[it] = np.dot(ds, Pn1n2_ps[1:] + Pn1n2_ps[:-1])
            
            cond1 = np.zeros((len(unicountvals1_cond), len(unicountvals2_cond)))
            for i in range (len(unicountvals1_cond)):
                for j in range (len(unicountvals2_cond)):
                    Pn1n2_ps = Pn1n2_s_cond1[:,i,j]*Ps
                    cond1[i,j] = np.dot(ds, Pn1n2_ps[1:] + Pn1n2_ps[:-1])            
            C = 1 - np.sum(cond1) 
    
                
            Energy= -np.dot(sparse_rep_counts_thresh, np.log(integ)) + N_obs*np.log(C)

                
            return Energy

        #--------------------------Compute-the-grid-----------------------------------------
    
        print('Calculation Surface : \n')

        npoints = 20 #to be chosen by 
        Avec = np.linspace(-3, 1, npoints)
        Bvec = np.logspace(-2, 2, npoints)

        LSurface =np.zeros((len(Bvec),len(Avec)))
        for i in range(len(Bvec)):
            for j in range(len(Avec)):
                LSurface[i, j]=  - cost([Avec[j], Bvec[i]])
        
        Amesh, Bmesh = np.meshgrid(Avec,Bvec)
        a,b = np.where(LSurface == np.max(LSurface))
            
        #---------------------------------Optimization-----------------------------------------
        
        optA = Amesh[a[0],b[0]]
        optB = Bmesh[a[0],b[0]]
                      
        print('polish parameter estimate from '+ str(optA)+' '+str(optB))
        initparas=(optA,optB)  
        #initparas =(-0.5, 2)
        
        bnds = ((None, None), (0.001, None))

        global curr_iter 
        curr_iter = int(0)

        #callback_partial = partial(self._callbackFdiffexpr_1D, n_thresh = n_thresh, n_lim = n_lim)
        callback_partial = self._callbackFdiffexpr_1D_simple

        outstruct = minimize(cost, initparas, method='SLSQP', callback=callback_partial, tol=1e-8,options={'ftol':1e-8,'disp': True,'maxiter':300}, bounds=bnds)

        return outstruct.x

    def learning_diffusion_1D_two_cond(self, f_thresh, n_thresh_2):

        paras_1 = self.paras_1
        paras_2 = self.paras_2
        freq_dtype = self.freq_dtype
        nfbins = self.nfbins
        method = self.method
        n_lim = 1000

        df_thresh= super().data_thresh_freq_dyn(f_thresh, n_thresh_2)
        sparse_rep_thresh = super().get_sparserep(df_thresh)
        indn1_thresh, indn2_thresh, sparse_rep_counts_thresh, unicountvals_1_thresh, unicountvals_2_thresh, NreadsI_thresh, NreadsII_thresh = sparse_rep_thresh

        NreadsI = self.NreadsI
        NreadsII = self.NreadsII

        n_thresh_1 = np.floor(f_thresh*NreadsI)
        n_thresh_2 = 1

        logrhofvec = self.logrhofvec
        logfvec = self.logfvec
        logfvecsecond = self.logfvecsecond
        svec = self.svec
        f2s_step = self.f2s_step

        #====CONDITION

        #SUM-1: P(n_1<n_tresh_1)
        unicountvals1_cond_SUM_1 = np.arange(n_thresh_1)
        unicountvals2_cond_SUM_1 = np.arange(n_lim)
        logPn1_f_cond_SUM_1 = self.get_logPn_f_condition_freq(f_thresh, 'first', unicountvals1_cond_SUM_1, unicountvals2_cond_SUM_1)
        logPn2_f_cond_SUM_1 = self.get_logPn_f_condition_freq(f_thresh, 'second', unicountvals1_cond_SUM_1, unicountvals2_cond_SUM_1)

        #SUM-2: P(n_2<n_tresh_2)
        unicountvals1_cond_SUM_2 = np.arange(n_lim)
        unicountvals2_cond_SUM_2 = np.arange(n_thresh_2)
        logPn1_f_cond_SUM_2 = self.get_logPn_f_condition_freq(f_thresh, 'first', unicountvals1_cond_SUM_2, unicountvals2_cond_SUM_2)
        logPn2_f_cond_SUM_2 = self.get_logPn_f_condition_freq(f_thresh, 'second', unicountvals1_cond_SUM_2, unicountvals2_cond_SUM_2)

        #SUM-3 : P(n_1<n_tresh_1,n_2<n_thresh_2)
        unicountvals1_cond_SUM_3 = np.arange(n_thresh_1)
        unicountvals2_cond_SUM_3 = np.arange(n_thresh_2)
        logPn1_f_cond_SUM_3 = self.get_logPn_f_condition_freq(f_thresh, 'first', unicountvals1_cond_SUM_3, unicountvals2_cond_SUM_3)
        logPn2_f_cond_SUM_3 = self.get_logPn_f_condition_freq(f_thresh, 'second', unicountvals1_cond_SUM_3, unicountvals2_cond_SUM_3)


        dlogfby2=np.diff(logfvec)/2


        #SUM-1: P(n_1<n_tresh_1)
        Pn1n2_s_cond_SUM_1=np.zeros((len(svec),len(unicountvals1_cond_SUM_1), len(unicountvals2_cond_SUM_1))) 
        for k in tqdm(range(len(svec))):
            for i in range (len(unicountvals1_cond_SUM_1)):
                for j in range (len(unicountvals2_cond_SUM_1)):
                    integ = np.exp(logrhofvec+logPn2_f_cond_SUM_1[f2s_step*k:(f2s_step*k+len(logfvec)),j]+logPn1_f_cond_SUM_1[:,i]+ logfvec)
                    Pn1n2_s_cond_SUM_1[k, i,j] = np.dot(dlogfby2, integ[1:] + integ[:-1])
        print('SUM-1 HAS BEEN COMPUTED')

        #SUM-2: P(n_2<n_tresh_2)
        Pn1n2_s_cond_SUM_2=np.zeros((len(svec),len(unicountvals1_cond_SUM_2), len(unicountvals2_cond_SUM_2))) 
        for k in tqdm(range(len(svec))):
            for i in range (len(unicountvals1_cond_SUM_2)):
                for j in range (len(unicountvals2_cond_SUM_2)):
                    integ = np.exp(logrhofvec+logPn2_f_cond_SUM_2[f2s_step*k:(f2s_step*k+len(logfvec)),j]+logPn1_f_cond_SUM_2[:,i]+ logfvec)
                    Pn1n2_s_cond_SUM_2[k, i,j] = np.dot(dlogfby2, integ[1:] + integ[:-1])
        print('SUM-2 HAS BEEN COMPUTED')

        #SUM-3 : P(n_1<n_tresh_1,n_2<n_thresh_2)
        Pn1n2_s_cond_SUM_3=np.zeros((len(svec),len(unicountvals1_cond_SUM_3), len(unicountvals2_cond_SUM_3))) 
        for k in tqdm(range(len(svec))):
            for i in range (len(unicountvals1_cond_SUM_3)):
                for j in range (len(unicountvals2_cond_SUM_3)):
                    integ = np.exp(logrhofvec+logPn2_f_cond_SUM_3[f2s_step*k:(f2s_step*k+len(logfvec)),j]+logPn1_f_cond_SUM_3[:,i]+ logfvec)
                    Pn1n2_s_cond_SUM_3[k, i,j] = np.dot(dlogfby2, integ[1:] + integ[:-1])
        print('SUM-3 HAS BEEN COMPUTED')



        logPn1_f = self.get_logPn_f_new(f_thresh, n_thresh_2, 'initial')
        logPn2_f = self.get_logPn_f_new(f_thresh, n_thresh_2, 'final')

        Pn1n2_s=np.zeros((len(svec),len(sparse_rep_counts_thresh))) 
        for s_it,s in enumerate(svec):
            for it,(n1_it, n2_it) in enumerate(zip(indn1_thresh,indn2_thresh)):
                integ = np.exp(logrhofvec+logPn2_f[f2s_step*s_it:(f2s_step*s_it+len(logfvec)),n2_it]+logPn1_f[:,n1_it]+ logfvec)
                Pn1n2_s[s_it, it] = np.dot(dlogfby2,integ[1:] + integ[:-1])
        print('DATA-OK')


        ds = np.diff(svec)/2
        N_obs = np.sum(sparse_rep_counts_thresh)

        def cost(PARAS):

            A = PARAS[0]
            B = PARAS[1]


            integ = np.zeros(len(sparse_rep_counts_thresh))
            Ps = self.get_Ps(A, B)
            for it, (n1_it, n2_it) in enumerate(zip(indn1_thresh, indn2_thresh)):
                Pn1n2_ps = Pn1n2_s[:,it]*Ps
                integ[it] = np.dot(ds, Pn1n2_ps[1:] + Pn1n2_ps[:-1])
            
            cond_SUM_1 = np.zeros((len(unicountvals1_cond_SUM_1), len(unicountvals2_cond_SUM_1)))
            for i in range (len(unicountvals1_cond_SUM_1)):
                for j in range (len(unicountvals2_cond_SUM_1)):
                    Pn1n2_ps = Pn1n2_s_cond_SUM_1[:,i,j]*Ps
                    cond_SUM_1[i,j] = np.dot(ds, Pn1n2_ps[1:] + Pn1n2_ps[:-1]) 

            cond_SUM_2 = np.zeros((len(unicountvals1_cond_SUM_2), len(unicountvals2_cond_SUM_2)))
            for i in range (len(unicountvals1_cond_SUM_2)):
                for j in range (len(unicountvals2_cond_SUM_2)):
                    Pn1n2_ps = Pn1n2_s_cond_SUM_2[:,i,j]*Ps
                    cond_SUM_2[i,j] = np.dot(ds, Pn1n2_ps[1:] + Pn1n2_ps[:-1]) 

            cond_SUM_3 = np.zeros((len(unicountvals1_cond_SUM_3), len(unicountvals2_cond_SUM_3)))
            for i in range (len(unicountvals1_cond_SUM_3)):
                for j in range (len(unicountvals2_cond_SUM_3)):
                    Pn1n2_ps = Pn1n2_s_cond_SUM_3[:,i,j]*Ps
                    cond_SUM_3[i,j] = np.dot(ds, Pn1n2_ps[1:] + Pn1n2_ps[:-1]) 

            SUM_1 = np.sum(cond_SUM_1) #SUM_1 : P(n_1<n_tresh_1)
            SUM_2 = np.sum(cond_SUM_2) #SUM_2 : P(n_2<n_tresh_2)
            SUM_3 = np.sum(cond_SUM_3) #SUM_3 : P(n_1<n_tresh_1,n_2<n_thresh_2)

            C = 1 - SUM_1 - SUM_2 + SUM_3 # C : 1 - P(n_1<n_tresh_1 or n_2<n_tresh_2) : P(n_1>=n_tresh_1 and n_2>=n_tresh_2)

            Energy= -np.dot(sparse_rep_counts_thresh, np.log(integ)) + N_obs*np.log(C)

                
            return Energy

        #--------------------------Compute-the-grid-----------------------------------------
    
        print('NO-Calculation Surface : \n')

        #npoints = 20 #to be chosen by 
        #Avec = np.linspace(-3, 1, npoints)
        #Bvec = np.logspace(-2, 2, npoints)

        #LSurface =np.zeros((len(Bvec),len(Avec)))
        #for i in range(len(Bvec)):
        #    for j in range(len(Avec)):
        #        LSurface[i, j]=  - cost([Avec[j], Bvec[i]])
        
        #Amesh, Bmesh = np.meshgrid(Avec,Bvec)
        #a,b = np.where(LSurface == np.max(LSurface))
            
        #---------------------------------Optimization-----------------------------------------
        
        #optA = Amesh[a[0],b[0]]
        #optB = Bmesh[a[0],b[0]]
                      
        #print('polish parameter estimate from '+ str(optA)+' '+str(optB))
        #initparas=(optA,optB)  
        initparas =(-0.5, 2)
        
        bnds = ((None, None), (0.001, None))

        global curr_iter 
        curr_iter = int(0)

        #callback_partial = partial(self._callbackFdiffexpr_1D, n_thresh = n_thresh, n_lim = n_lim)
        callback_partial = self._callbackFdiffexpr_1D_simple

        outstruct = minimize(cost, initparas, method='SLSQP', callback=callback_partial, tol=1e-8,options={'ftol':1e-8,'disp': True,'maxiter':300}, bounds=bnds)

        return outstruct.x

    def learning_diffusion_1D_db(self, n_threshminus, n_threshplus, n_lim):

        paras_1 = self.paras_1
        paras_2 = self.paras_2
        freq_dtype = self.freq_dtype
        nfbins = self.nfbins
        method = self.method

        df_thresh= self.data_thresh_db(n_threshminus, n_threshplus)
        sparse_rep_thresh = self.get_sparserep(df_thresh)
        indn1_thresh, indn2_thresh, sparse_rep_counts_thresh, unicountvals_1_thresh, unicountvals_2_thresh, NreadsI_thresh, NreadsII_thresh = sparse_rep_thresh

        logrhofvec = self.logrhofvec
        logfvec = self.logfvec
        logfvecsecond = self.logfvecsecond
        svec = self.svec
        f2s_step = self.f2s_step

        unicountvals1_cond = np.arange(n_threshminus, n_threshplus+1)

        if method == 'negative_binomial':
            unicountvals2_cond = np.arange(n_lim)

        elif method == 'poisson':
            unicountvals2_cond = np.arange(n_lim*10)

        logPn1_f_cond = self.get_logPn_f_condition_db('initial', unicountvals1_cond, unicountvals2_cond)
        logPn2_f_cond = self.get_logPn_f_condition_db('final', unicountvals1_cond, unicountvals2_cond)

        dlogfby2=np.diff(logfvec)/2

        Pn1n2_s_cond1=np.zeros((len(svec),len(unicountvals1_cond), len(unicountvals2_cond))) 
        for k in tqdm(range(len(svec))):
            for i in range (len(unicountvals1_cond)):
                for j in range (len(unicountvals2_cond)):
                    integ = np.exp(logrhofvec+logPn2_f_cond[f2s_step*k:(f2s_step*k+len(logfvec)),j]+logPn1_f_cond[:,i]+ logfvec)
                    Pn1n2_s_cond1[k, i, j] = np.dot(dlogfby2,integ[1:] + integ[:-1])

        logPn1_f = self.get_logPn_f_db(n_threshminus, n_threshplus, 'initial')
        logPn2_f = self.get_logPn_f_db(n_threshminus, n_threshplus, 'final')

        Pn1n2_s=np.zeros((len(svec),len(sparse_rep_counts_thresh))) 
        for s_it,s in enumerate(svec):
            for it,(n1_it, n2_it) in enumerate(zip(indn1_thresh,indn2_thresh)):
                integ = np.exp(logrhofvec+logPn2_f[f2s_step*s_it:(f2s_step*s_it+len(logfvec)),n2_it]+logPn1_f[:,n1_it]+ logfvec)
                Pn1n2_s[s_it, it] = np.dot(dlogfby2,integ[1:] + integ[:-1])


        ds = np.diff(svec)/2
        N_obs = np.sum(sparse_rep_counts_thresh)


        def cost(PARAS):

            A = PARAS[0]
            B = PARAS[1]


            integ = np.zeros(len(sparse_rep_counts_thresh))
            Ps = self.get_Ps(A, B)
            for it, (n1_it, n2_it) in enumerate(zip(indn1_thresh, indn2_thresh)):
                Pn1n2_ps = Pn1n2_s[:,it]*Ps
                integ[it] = np.dot(ds, Pn1n2_ps[1:] + Pn1n2_ps[:-1])
                
            cond1 = np.zeros((len(unicountvals1_cond), len(unicountvals2_cond)))
            for i in range (len(unicountvals1_cond)):
                for j in range (len(unicountvals2_cond)):
                    Pn1n2_ps = Pn1n2_s_cond1[:,i,j]*Ps
                    cond1[i,j] = np.dot(ds, Pn1n2_ps[1:] + Pn1n2_ps[:-1])
            #print(cond1[0,0])
            print(cond1[-1,-1], np.shape(cond1))
                
            C = np.sum(cond1) #- np.sum(cond2) + Pn0n0 
                #print("c = " + str(C))
                    
            Energy= -np.dot(sparse_rep_counts_thresh, np.log(integ)) + N_obs*np.log(C)

                    
            return Energy

        #--------------------------Compute-the-grid-----------------------------------------
        
        print('Calculation Surface : \n')

        npoints = 20 #to be chosen by 
        Avec = np.linspace(-3, 1, npoints)
        Bvec = np.logspace(-2, 2, npoints)

        LSurface =np.zeros((len(Bvec),len(Avec)))
        for i in range(len(Bvec)):
            for j in range(len(Avec)):
                LSurface[i, j]=  - cost([Avec[j], Bvec[i]])
            
        Amesh, Bmesh = np.meshgrid(Avec,Bvec)
        a,b = np.where(LSurface == np.max(LSurface))
                
        #---------------------------------Optimization-----------------------------------------
            
        optA = Amesh[a[0],b[0]]
        optB = Bmesh[a[0],b[0]]
                          
        print('polish parameter estimate from '+ str(optA)+' '+str(optB))
        initparas=(optA,optB)  
            
        bnds = ((None, None), (0.001, None))

        global curr_iter 
        curr_iter = int(0)

        callback_partial = self._callbackFdiffexpr_simple

        outstruct = minimize(cost, initparas, method='SLSQP', callback=callback_partial, tol=1e-6,options={'ftol':1e-6 ,'disp': True,'maxiter':10000}, bounds=bnds)

        return outstruct.x

    def learning_diffusion_1D_db_cond_n2(self, n_threshminus, n_threshplus, n_2, n_lim):

        paras_1 = self.paras_1
        paras_2 = self.paras_2
        freq_dtype = self.freq_dtype
        nfbins = self.nfbins
        method = self.method

        df_thresh= self.data_thresh_db_cond_n2(n_threshminus, n_threshplus, n_2)
        sparse_rep_thresh = self.get_sparserep(df_thresh)
        indn1_thresh, indn2_thresh, sparse_rep_counts_thresh, unicountvals_1_thresh, unicountvals_2_thresh, NreadsI_thresh, NreadsII_thresh = sparse_rep_thresh
        print(unicountvals_1_thresh, unicountvals_2_thresh)

        logrhofvec = self.logrhofvec
        logfvec = self.logfvec
        logfvecsecond = self.logfvecsecond
        svec = self.svec
        f2s_step = self.f2s_step

        unicountvals1_cond = np.arange(n_threshminus+1, n_threshplus+1)

        if method == 'negative_binomial':
            unicountvals2_cond = np.arange(n_2,n_lim)

        elif method == 'poisson':
            unicountvals2_cond = np.arange(n_2,n_lim*10)

        logPn1_f_cond = self.get_logPn_f_condition_db('initial', unicountvals1_cond, unicountvals2_cond)
        logPn2_f_cond = self.get_logPn_f_condition_db('final', unicountvals1_cond, unicountvals2_cond)

        dlogfby2=np.diff(logfvec)/2

        Pn1n2_s_cond1=np.zeros((len(svec),len(unicountvals1_cond), len(unicountvals2_cond))) 
        for k in tqdm(range(len(svec))):
            for i in range (len(unicountvals1_cond)):
                for j in range (len(unicountvals2_cond)):
                    integ = np.exp(logrhofvec+logPn2_f_cond[f2s_step*k:(f2s_step*k+len(logfvec)),j]+logPn1_f_cond[:,i]+ logfvec)
                    Pn1n2_s_cond1[k, i,j] = np.dot(dlogfby2,integ[1:] + integ[:-1])

        logPn1_f = self.get_logPn_f_db_n2_cond(n_threshminus, n_threshplus, n_2, 'initial')
        logPn2_f = self.get_logPn_f_db_n2_cond(n_threshminus, n_threshplus, n_2, 'final')

        Pn1n2_s=np.zeros((len(svec),len(sparse_rep_counts_thresh))) 
        print(len(sparse_rep_counts_thresh))
        for s_it,s in enumerate(svec):
            for it,(n1_it, n2_it) in enumerate(zip(indn1_thresh,indn2_thresh)):
                integ = np.exp(logrhofvec+logPn2_f[f2s_step*s_it:(f2s_step*s_it+len(logfvec)),n2_it]+logPn1_f[:,n1_it]+ logfvec)
                Pn1n2_s[s_it, it] = np.dot(dlogfby2,integ[1:] + integ[:-1])


        ds = np.diff(svec)/2
        N_obs = np.sum(sparse_rep_counts_thresh)


        def cost(PARAS):

            A = PARAS[0]
            B = PARAS[1]


            integ = np.zeros(len(sparse_rep_counts_thresh))
            Ps = self.get_Ps(A, B)
            for it, (n1_it, n2_it) in enumerate(zip(indn1_thresh, indn2_thresh)):
                Pn1n2_ps = Pn1n2_s[:,it]*Ps
                integ[it] = np.dot(ds, Pn1n2_ps[1:] + Pn1n2_ps[:-1])
                
            cond1 = np.zeros((len(unicountvals1_cond), len(unicountvals2_cond)))
            for i in range (len(unicountvals1_cond)):
                for j in range (len(unicountvals2_cond)):
                    Pn1n2_ps = Pn1n2_s_cond1[:,i,j]*Ps
                    cond1[i,j] = np.dot(ds, Pn1n2_ps[1:] + Pn1n2_ps[:-1])
            #print(cond1[0,0])
            print(cond1[-1,-1], np.shape(cond1))
                
            C = np.sum(cond1) #- np.sum(cond2) + Pn0n0 
                #print("c = " + str(C))
                    
            Energy= -np.dot(sparse_rep_counts_thresh, np.log(integ)) + N_obs*np.log(C)

                    
            return Energy

        #--------------------------Compute-the-grid-----------------------------------------
        
        print('Calculation Surface : \n')

        npoints = 20 #to be chosen by 
        Avec = np.linspace(-3, 1, npoints)
        Bvec = np.logspace(-2, 2, npoints)

        LSurface =np.zeros((len(Bvec),len(Avec)))
        for i in range(len(Bvec)):
            for j in range(len(Avec)):
                LSurface[i, j]=  - cost([Avec[j], Bvec[i]])
            
        Amesh, Bmesh = np.meshgrid(Avec,Bvec)
        a,b = np.where(LSurface == np.max(LSurface))
                
        #---------------------------------Optimization-----------------------------------------
            
        optA = Amesh[a[0],b[0]]
        optB = Bmesh[a[0],b[0]]
                          
        print('polish parameter estimate from '+ str(optA)+' '+str(optB))
        initparas=(optA,optB)  
            
        bnds = ((None, None), (0.001, None))

        global curr_iter 
        curr_iter = int(0)

        callback_partial = self._callbackFdiffexpr_simple

        outstruct = minimize(cost, initparas, method='SLSQP', callback=callback_partial, tol=1e-6,options={'ftol':1e-6 ,'disp': True,'maxiter':10000}, bounds=bnds)

        return outstruct.x

    def learning_diffusion_1D_db_cond_fth_n2(self, f_threshminus, f_threshplus, n_2, n_lim):

        paras_1 = self.paras_1
        paras_2 = self.paras_2
        freq_dtype = self.freq_dtype
        nfbins = self.nfbins
        method = self.method

        NreadsI = self.NreadsI
        NreadsII = self.NreadsII

        n_threshminus = np.floor(f_threshminus*NreadsI)
        n_threshplus = np.floor(f_threshplus*NreadsI)

        df_thresh= self.data_thresh_db_cond_fth_n2(f_threshminus, f_threshplus, n_2)
        sparse_rep_thresh = self.get_sparserep(df_thresh)
        indn1_thresh, indn2_thresh, sparse_rep_counts_thresh, unicountvals_1_thresh, unicountvals_2_thresh, NreadsI_thresh, NreadsII_thresh = sparse_rep_thresh

        logrhofvec = self.logrhofvec
        logfvec = self.logfvec
        logfvecsecond = self.logfvecsecond
        svec = self.svec
        f2s_step = self.f2s_step

        unicountvals1_cond = np.arange(n_threshminus, n_threshplus+1)
        print(type(unicountvals1_cond))

        unicountvals1_cond = unicountvals1_cond.astype(int)
        print(type(unicountvals1_cond))


        if method == 'negative_binomial':
            unicountvals2_cond = np.arange(n_2, n_lim)

        elif method == 'poisson':
            unicountvals2_cond = np.arange(n_2,n_lim*10)

        logPn1_f_cond = self.get_logPn_f_condition_db('initial', unicountvals1_cond, unicountvals2_cond)
        logPn2_f_cond = self.get_logPn_f_condition_db('final', unicountvals1_cond, unicountvals2_cond)

        dlogfby2=np.diff(logfvec)/2

        Pn1n2_s_cond1=np.zeros((len(svec),len(unicountvals1_cond), len(unicountvals2_cond))) 
        for k in tqdm(range(len(svec))):
            for i in range (len(unicountvals1_cond)):
                for j in range (len(unicountvals2_cond)):
                    integ = np.exp(logrhofvec+logPn2_f_cond[f2s_step*k:(f2s_step*k+len(logfvec)),j]+logPn1_f_cond[:,i]+ logfvec)
                    Pn1n2_s_cond1[k, i,j] = np.dot(dlogfby2,integ[1:] + integ[:-1])

        logPn1_f = self.get_logPn_f_db_n2_cond(n_threshminus, n_threshplus, n_2, 'initial')
        logPn2_f = self.get_logPn_f_db_n2_cond(n_threshminus, n_threshplus, n_2, 'final')

        Pn1n2_s=np.zeros((len(svec),len(sparse_rep_counts_thresh))) 
        for s_it,s in enumerate(svec):
            for it,(n1_it, n2_it) in enumerate(zip(indn1_thresh,indn2_thresh)):
                integ = np.exp(logrhofvec+logPn2_f[f2s_step*s_it:(f2s_step*s_it+len(logfvec)),n2_it]+logPn1_f[:,n1_it]+ logfvec)
                Pn1n2_s[s_it, it] = np.dot(dlogfby2,integ[1:] + integ[:-1])


        ds = np.diff(svec)/2
        N_obs = np.sum(sparse_rep_counts_thresh)


        def cost(PARAS):

            A = PARAS[0]
            B = PARAS[1]


            integ = np.zeros(len(sparse_rep_counts_thresh))
            Ps = self.get_Ps(A, B)
            for it, (n1_it, n2_it) in enumerate(zip(indn1_thresh, indn2_thresh)):
                Pn1n2_ps = Pn1n2_s[:,it]*Ps
                integ[it] = np.dot(ds, Pn1n2_ps[1:] + Pn1n2_ps[:-1])
                
            cond1 = np.zeros((len(unicountvals1_cond), len(unicountvals2_cond)))
            for i in range (len(unicountvals1_cond)):
                for j in range (len(unicountvals2_cond)):
                    Pn1n2_ps = Pn1n2_s_cond1[:,i,j]*Ps
                    cond1[i,j] = np.dot(ds, Pn1n2_ps[1:] + Pn1n2_ps[:-1])
            #print(cond1[0,0])
            print(cond1[-1,-1], np.shape(cond1))
                
            C = np.sum(cond1) #- np.sum(cond2) + Pn0n0 
                #print("c = " + str(C))
                    
            Energy= -np.dot(sparse_rep_counts_thresh, np.log(integ)) + N_obs*np.log(C)

                    
            return Energy

        #--------------------------Compute-the-grid-----------------------------------------
        
        print('Calculation Surface : \n')

        npoints = 20 #to be chosen by 
        Avec = np.linspace(-3, 1, npoints)
        Bvec = np.logspace(-2, 2, npoints)

        LSurface =np.zeros((len(Bvec),len(Avec)))
        for i in range(len(Bvec)):
            for j in range(len(Avec)):
                LSurface[i, j]=  - cost([Avec[j], Bvec[i]])
            
        Amesh, Bmesh = np.meshgrid(Avec,Bvec)
        a,b = np.where(LSurface == np.max(LSurface))
                
        #---------------------------------Optimization-----------------------------------------
            
        optA = Amesh[a[0],b[0]]
        optB = Bmesh[a[0],b[0]]
                          
        print('polish parameter estimate from '+ str(optA)+' '+str(optB))
        initparas=(optA,optB)  
            
        bnds = ((None, None), (0.001, None))

        global curr_iter 
        curr_iter = int(0)

        callback_partial = self._callbackFdiffexpr_simple

        outstruct = minimize(cost, initparas, method='SLSQP', callback=callback_partial, tol=1e-6,options={'ftol':1e-6 ,'disp': True,'maxiter':10000}, bounds=bnds)

        return outstruct.x

    def get_Posterior_utils(self, PARAS, f_thresh):

        paras_1 = self.paras_1
        paras_2 = self.paras_2
        freq_dtype = self.freq_dtype
        nfbins = self.nfbins
        method = self.method
        n_thresh_2 = 1


        df_thresh= super().data_thresh_freq_dyn(f_thresh, n_thresh_2)
        sparse_rep_thresh = super().get_sparserep(df_thresh)
        indn1_thresh, indn2_thresh, sparse_rep_counts_thresh, unicountvals_1_thresh, unicountvals_2_thresh, NreadsI_thresh, NreadsII_thresh = sparse_rep_thresh


        logrhofvec = self.logrhofvec
        logfvec = self.logfvec
        logfvecsecond = self.logfvecsecond
        svec = self.svec
        f2s_step = self.f2s_step


        dlogfby2=np.diff(logfvec)/2 

        logPn1_f = self.get_logPn_f_new(f_thresh, n_thresh_2, 'initial')
        logPn2_f = self.get_logPn_f_new(f_thresh, n_thresh_2, 'final')

        Pn1n2_s=np.zeros((len(svec),len(sparse_rep_counts_thresh))) 
        for s_it,s in enumerate(svec):
            for it,(n1_it, n2_it) in enumerate(zip(indn1_thresh,indn2_thresh)):
                integ = np.exp(logrhofvec+logPn2_f[f2s_step*s_it:(f2s_step*s_it+len(logfvec)),n2_it]+logPn1_f[:,n1_it]+ logfvec)
                Pn1n2_s[s_it, it] = np.dot(dlogfby2,integ[1:] + integ[:-1])


        ds = np.diff(svec)/2
        N_obs = np.sum(sparse_rep_counts_thresh)

        A = PARAS[0]
        B = PARAS[1]

        Ps = self.get_Ps(A, B)


        return Ps, Pn1n2_s

    def get_Posterior_matrix(self, PARAS, f_thresh):

        Ps, Pn1n2S = self.get_Posterior_utils(PARAS, f_thresh)
    
        (a,b) = np.shape(Pn1n2S)
        PsTile = np.tile(Ps, (b,1)).T
        M = np.multiply(PsTile, Pn1n2S)
        svec = self.svec
        ds = np.tile(np.diff(svec)/2,(b,1)).T 

        intPs = np.sum(np.multiply(ds, M[1:, :] + M[:-1,:]), axis = 0)
        MPosterior = np.divide(M, intPs)
    
        return MPosterior


    def get_Ps_full(self, PARAS, f_thresh):

        paras_1 = self.paras_1
        paras_2 = self.paras_2
        freq_dtype = self.freq_dtype
        nfbins = self.nfbins
        method = self.method
        n_lim = 500
        n_thresh_2 = 1


        df_thresh= super().data_thresh_freq_dyn(f_thresh, n_thresh_2)
        sparse_rep_thresh = super().get_sparserep(df_thresh)
        indn1_thresh, indn2_thresh, sparse_rep_counts_thresh, unicountvals_1_thresh, unicountvals_2_thresh, NreadsI_thresh, NreadsII_thresh = sparse_rep_thresh

        NreadsI = self.NreadsI
        NreadsII = self.NreadsII

        n_thresh_1 = np.floor(f_thresh*NreadsI)
        #n_thresh_2 = 1

        logrhofvec = self.logrhofvec
        logfvec = self.logfvec
        logfvecsecond = self.logfvecsecond
        svec = self.svec
        f2s_step = self.f2s_step

        #====CONDITION

        #SUM-1: P(n_1<n_tresh_1)
        unicountvals1_cond_SUM_1 = np.arange(n_thresh_1)
        unicountvals2_cond_SUM_1 = np.arange(n_lim)
        logPn1_f_cond_SUM_1 = self.get_logPn_f_condition_freq(f_thresh, 'first', unicountvals1_cond_SUM_1, unicountvals2_cond_SUM_1)
        logPn2_f_cond_SUM_1 = self.get_logPn_f_condition_freq(f_thresh, 'second', unicountvals1_cond_SUM_1, unicountvals2_cond_SUM_1)

        #SUM-2: P(n_2<n_tresh_2)
        unicountvals1_cond_SUM_2 = np.arange(n_lim)
        unicountvals2_cond_SUM_2 = np.arange(n_thresh_2)
        logPn1_f_cond_SUM_2 = self.get_logPn_f_condition_freq(f_thresh, 'first', unicountvals1_cond_SUM_2, unicountvals2_cond_SUM_2)
        logPn2_f_cond_SUM_2 = self.get_logPn_f_condition_freq(f_thresh, 'second', unicountvals1_cond_SUM_2, unicountvals2_cond_SUM_2)

        #SUM-3 : P(n_1<n_tresh_1,n_2<n_thresh_2)
        unicountvals1_cond_SUM_3 = np.arange(n_thresh_1)
        unicountvals2_cond_SUM_3 = np.arange(n_thresh_2)
        logPn1_f_cond_SUM_3 = self.get_logPn_f_condition_freq(f_thresh, 'first', unicountvals1_cond_SUM_3, unicountvals2_cond_SUM_3)
        logPn2_f_cond_SUM_3 = self.get_logPn_f_condition_freq(f_thresh, 'second', unicountvals1_cond_SUM_3, unicountvals2_cond_SUM_3)


        dlogfby2=np.diff(logfvec)/2


        #SUM-1: P(n_1<n_tresh_1)
        Pn1n2_s_cond_SUM_1=np.zeros((len(svec),len(unicountvals1_cond_SUM_1), len(unicountvals2_cond_SUM_1))) 
        for k in tqdm(range(len(svec))):
            for i in range (len(unicountvals1_cond_SUM_1)):
                for j in range (len(unicountvals2_cond_SUM_1)):
                    integ = np.exp(logrhofvec+logPn2_f_cond_SUM_1[f2s_step*k:(f2s_step*k+len(logfvec)),j]+logPn1_f_cond_SUM_1[:,i]+ logfvec)
                    Pn1n2_s_cond_SUM_1[k, i,j] = np.dot(dlogfby2, integ[1:] + integ[:-1])

        #SUM-2: P(n_2<n_tresh_2)
        Pn1n2_s_cond_SUM_2=np.zeros((len(svec),len(unicountvals1_cond_SUM_2), len(unicountvals2_cond_SUM_2))) 
        for k in tqdm(range(len(svec))):
            for i in range (len(unicountvals1_cond_SUM_2)):
                for j in range (len(unicountvals2_cond_SUM_2)):
                    integ = np.exp(logrhofvec+logPn2_f_cond_SUM_2[f2s_step*k:(f2s_step*k+len(logfvec)),j]+logPn1_f_cond_SUM_2[:,i]+ logfvec)
                    Pn1n2_s_cond_SUM_2[k, i,j] = np.dot(dlogfby2, integ[1:] + integ[:-1])

        #SUM-3 : P(n_1<n_tresh_1,n_2<n_thresh_2)
        Pn1n2_s_cond_SUM_3=np.zeros((len(svec),len(unicountvals1_cond_SUM_3), len(unicountvals2_cond_SUM_3))) 
        for k in tqdm(range(len(svec))):
            for i in range (len(unicountvals1_cond_SUM_3)):
                for j in range (len(unicountvals2_cond_SUM_3)):
                    integ = np.exp(logrhofvec+logPn2_f_cond_SUM_3[f2s_step*k:(f2s_step*k+len(logfvec)),j]+logPn1_f_cond_SUM_3[:,i]+ logfvec)
                    Pn1n2_s_cond_SUM_3[k, i,j] = np.dot(dlogfby2, integ[1:] + integ[:-1])

        logPn1_f = self.get_logPn_f_new(f_thresh, n_thresh_2,  'initial')
        logPn2_f = self.get_logPn_f_new(f_thresh, n_thresh_2, 'final')

        Pn1n2_s=np.zeros((len(svec),len(sparse_rep_counts_thresh))) 
        for s_it,s in enumerate(svec):
            for it,(n1_it, n2_it) in enumerate(zip(indn1_thresh,indn2_thresh)):
                integ = np.exp(logrhofvec+logPn2_f[f2s_step*s_it:(f2s_step*s_it+len(logfvec)),n2_it]+logPn1_f[:,n1_it]+ logfvec)
                Pn1n2_s[s_it, it] = np.dot(dlogfby2,integ[1:] + integ[:-1])


        ds = np.diff(svec)/2
        N_obs = np.sum(sparse_rep_counts_thresh)

        A = PARAS[0]
        B = PARAS[1]

        integ = np.zeros(len(sparse_rep_counts_thresh))
        Ps = self.get_Ps(A, B)
        for it, (n1_it, n2_it) in enumerate(zip(indn1_thresh, indn2_thresh)):
            Pn1n2_ps = Pn1n2_s[:,it]*Ps
            integ[it] = np.dot(ds, Pn1n2_ps[1:] + Pn1n2_ps[:-1])
            
        cond_SUM_1 = np.zeros((len(unicountvals1_cond_SUM_1), len(unicountvals2_cond_SUM_1)))
        for i in range (len(unicountvals1_cond_SUM_1)):
            for j in range (len(unicountvals2_cond_SUM_1)):
                Pn1n2_ps = Pn1n2_s_cond_SUM_1[:,i,j]*Ps
                cond_SUM_1[i,j] = np.dot(ds, Pn1n2_ps[1:] + Pn1n2_ps[:-1]) 

        cond_SUM_2 = np.zeros((len(unicountvals1_cond_SUM_2), len(unicountvals2_cond_SUM_2)))
        for i in range (len(unicountvals1_cond_SUM_2)):
            for j in range (len(unicountvals2_cond_SUM_2)):
                Pn1n2_ps = Pn1n2_s_cond_SUM_2[:,i,j]*Ps
                cond_SUM_2[i,j] = np.dot(ds, Pn1n2_ps[1:] + Pn1n2_ps[:-1]) 

        cond_SUM_3 = np.zeros((len(unicountvals1_cond_SUM_3), len(unicountvals2_cond_SUM_3)))
        for i in range (len(unicountvals1_cond_SUM_3)):
            for j in range (len(unicountvals2_cond_SUM_3)):
                Pn1n2_ps = Pn1n2_s_cond_SUM_3[:,i,j]*Ps
                cond_SUM_3[i,j] = np.dot(ds, Pn1n2_ps[1:] + Pn1n2_ps[:-1]) 

        SUM_1 = np.sum(cond_SUM_1) #SUM_1 : P(n_1<n_tresh_1)
        SUM_2 = np.sum(cond_SUM_2) #SUM_2 : P(n_2<n_tresh_2)
        SUM_3 = np.sum(cond_SUM_3) #SUM_3 : P(n_1<n_tresh_1,n_2<n_thresh_2)

        C = 1 - SUM_1 - SUM_2 + SUM_3 # C : 1 - P(n_1<n_tresh_1 or n_2<n_tresh_2) : P(n_1>=n_tresh_1 and n_2>=n_tresh_2) 

        P_n1n2_thresh = integ /C
        MPosterior = self.get_Posterior_matrix(PARAS, f_thresh)
        
        (a,b) = np.shape(MPosterior)
        print(np.shape(MPosterior))
        Pn1n2ThreshFull = np.tile(P_n1n2_thresh, (a,1))
        print(np.shape(Pn1n2ThreshFull))

        PSTosum = np.multiply(MPosterior, Pn1n2ThreshFull)
        print(PSTosum.shape)

        sparseRepCountsThresh = np.tile(sparse_rep_counts_thresh, (a,1))

        print(sparseRepCountsThresh.shape)

        P_s_full = np.sum(np.multiply(sparseRepCountsThresh, PSTosum), axis = 1)
        np.sum(np.multiply(sparse_rep_counts_thresh,P_n1n2_thresh))

        PSFull = P_s_full/np.sum(np.multiply(sparse_rep_counts_thresh,P_n1n2_thresh))

        #PSFull = (1/N_obs)* np.sum(np.multiply(sparseRepCountsThresh, MPosterior), axis = 1)#np.sum(P_n1n2_thresh)

        return PSFull

    def P_persistency(self, A,B,t, f_thresh):

        """
        This function compute the Probability for a clone initially present in one sample to be present in an 
        other sample experient after t years knowing that initially the empirical clonal frequency is above a threshold
        """

        NreadsI = self.NreadsI #total number of reads at initial time
        paras_1 = self.paras_1
        paras_2 = self.paras_2
        n_lim = 1000
        logrhofvec = self.logrhofvec
        logfvec = self.logfvec
        logfvecsecond = self.logfvecsecond
        svec = self.svec
        f2s_step = self.f2s_step

        n_thresh_1 = np.floor(f_thresh*NreadsI).astype(int)
        print(n_thresh_1)

        #SUM-1: P(n_1>=n_tresh_1, n_2 = 0)
        unicountvals1_cond_SUM_1 = np.arange(n_thresh_1, n_lim)
        unicountvals2_cond_SUM_1 = np.arange(1)

        logPn1_f_cond_SUM_1 = self.get_logPn_f_condition_freq(f_thresh, 'first', unicountvals1_cond_SUM_1, unicountvals2_cond_SUM_1)
        logPn2_f_cond_SUM_1 = self.get_logPn_f_condition_freq(f_thresh, 'second', unicountvals1_cond_SUM_1, unicountvals2_cond_SUM_1)

        #SUM-2: P(n_1<n_tresh_1)
        unicountvals1_cond_SUM_2 = np.arange(n_thresh_1)
        unicountvals2_cond_SUM_2 = np.arange(n_lim)
        logPn1_f_cond_SUM_2 = self.get_logPn_f_condition_freq(f_thresh, 'first', unicountvals1_cond_SUM_2, unicountvals2_cond_SUM_2)
        logPn2_f_cond_SUM_2 = self.get_logPn_f_condition_freq(f_thresh, 'second', unicountvals1_cond_SUM_2, unicountvals2_cond_SUM_2)

        dlogfby2=np.diff(logfvec)/2
        ds = np.diff(svec)/2
        
        #SUM-1: P(n_1<n_tresh_1)
        Pn1n2_s_cond_SUM_1=np.zeros((len(svec),len(unicountvals1_cond_SUM_1), len(unicountvals2_cond_SUM_1))) 
        for k in tqdm(range(len(svec))):
            for i in range (len(unicountvals1_cond_SUM_1)):
                for j in range (len(unicountvals2_cond_SUM_1)):
                    integ = np.exp(logrhofvec+logPn2_f_cond_SUM_1[f2s_step*k:(f2s_step*k+len(logfvec)),j]+logPn1_f_cond_SUM_1[:,i]+ logfvec)
                    Pn1n2_s_cond_SUM_1[k, i,j] = np.dot(dlogfby2, integ[1:] + integ[:-1])

        #SUM-2: P(n_2=0)
        Pn1n2_s_cond_SUM_2=np.zeros((len(svec),len(unicountvals1_cond_SUM_2), len(unicountvals2_cond_SUM_2))) 
        for k in tqdm(range(len(svec))):
            for i in range (len(unicountvals1_cond_SUM_2)):
                for j in range (len(unicountvals2_cond_SUM_2)):
                    integ = np.exp(logrhofvec+logPn2_f_cond_SUM_2[f2s_step*k:(f2s_step*k+len(logfvec)),j]+logPn1_f_cond_SUM_2[:,i]+ logfvec)
                    Pn1n2_s_cond_SUM_2[k, i,j] = np.dot(dlogfby2, integ[1:] + integ[:-1])    

        Ps = self.get_Ps(A*t, B*t)

        cond_SUM_1 = np.zeros((len(unicountvals1_cond_SUM_1), len(unicountvals2_cond_SUM_1)))
        for i in range (len(unicountvals1_cond_SUM_1)):
            for j in range (len(unicountvals2_cond_SUM_1)):
                Pn1n2_ps = Pn1n2_s_cond_SUM_1[:,i,j]*Ps
                cond_SUM_1[i,j] = np.dot(ds, Pn1n2_ps[1:] + Pn1n2_ps[:-1]) 

        cond_SUM_2 = np.zeros((len(unicountvals1_cond_SUM_2), len(unicountvals2_cond_SUM_2)))
        for i in range (len(unicountvals1_cond_SUM_2)):
            for j in range (len(unicountvals2_cond_SUM_2)):
                Pn1n2_ps = Pn1n2_s_cond_SUM_2[:,i,j]*Ps
                cond_SUM_2[i,j] = np.dot(ds, Pn1n2_ps[1:] + Pn1n2_ps[:-1]) 

        #P(n_2=0, f_1>=f_thresh)
        PF1_thN2_0 = np.sum(cond_SUM_1)
        #P_f1_thresh
        P_F1_threh = 1-np.sum(cond_SUM_2)

        #P(n_2 = 0|f_1>f_thresh)
        RES = 1 - PF1_thN2_0/P_F1_threh

        return RES






class Inference_sampling(Inference_tools):


    def __init__(self, path, filename1 = None, filename2 = None, colnames1= None, colnames2 = None, HR = False):

        super().__init__(path, filename1 , filename2 , colnames1, colnames2, HR)
        self.freq_dtype = float



    def get_rhof(self, paras):
        '''
        generates power law (power is alpha_rho) clone frequency distribution over 
        freq_nbins discrete logarithmically spaced frequences between fmin and 1 of dtype freq_dtype
        Outputs log probabilities obtained at log frequencies'''
        
        freq_dtype = self.freq_dtype
        freq_nbins = self.nfbins
        paras_1 = paras
        alpha_rho = paras_1[0]
        fmin = np.power(10, paras_1[-1])
        fmax=1e0
        
        logfvec=np.linspace(np.log10(fmin),np.log10(fmax),freq_nbins)
        logfvec=np.array(np.log(np.power(10,logfvec)) ,dtype=freq_dtype).flatten()  
        logrhovec=logfvec*alpha_rho
        integ=np.exp(logrhovec+logfvec,dtype=freq_dtype)
        normconst=np.log(np.dot(np.diff(logfvec)/2.,integ[1:]+integ[:-1]))
        logrhovec-=normconst 
        return logrhovec, logfvec

    def get_logPn_f_noise(self, paras, f_thresh, replicate):
    
        if len(paras) == 4:
            method = 'negative_binomial'

        elif len(paras) == 2:
            method = 'poisson'
        
        logrhofvec, logfvec= self.get_rhof(paras)

        df_thresh= super().data_thresh_freq_noise(f_thresh)
        sparse_rep_thresh = super().get_sparserep(df_thresh)
        
        indn1_thresh, indn2_thresh, sparse_rep_counts_thresh, unicountvals_1_thresh, unicountvals_2_thresh, NreadsI_thresh, NreadsII_thresh = sparse_rep_thresh
        NreadsI = self.NreadsI
        NreadsII = self.NreadsII
        
        if replicate == 'first':
            unicounts, Nreads = unicountvals_1_thresh.astype(int), NreadsI
            #print(unicounts)
            logfvec_tmp=deepcopy(logfvec)

        elif replicate == 'second': 
            unicounts, Nreads = unicountvals_2_thresh.astype(int), NreadsII
            #print(unicounts)
            logfvec_tmp=deepcopy(logfvec)

        Pn_f=np.zeros((len(logfvec_tmp),len(unicounts)))
        
        if method == 'poisson':
            mean_n=Nreads*np.exp(logfvec_tmp)
            Pn_f= super().PoisPar(mean_n,unicounts)

        if method == 'negative_binomial':
            beta_mv = paras[1]
            alpha_mv =paras[2]
            mean_n=Nreads*np.exp(logfvec_tmp)
            var_n=mean_n+beta_mv*np.power(mean_n,alpha_mv)
            Pn_f= super().NegBinParMtr(mean_n,var_n,unicounts)
            
        return np.log(Pn_f)

        
    def get_logPn_f_condition_noise(self, paras, f_thresh, replicate, unicountvals1_cond, unicountvals2_cond):
        
        if len(paras) == 4:
            method = 'negative_binomial'

        elif len(paras) == 2:
            method = 'poisson'
        
        logrhofvec, logfvec= self.get_rhof(paras)
        
        NreadsI = self.NreadsI
        NreadsII = self.NreadsII
            
        if replicate == 'first':
            unicounts, Nreads = unicountvals1_cond.astype(int), NreadsI
            logfvec_tmp=deepcopy(logfvec)
            
        elif replicate == 'second': 
            unicounts, Nreads = unicountvals2_cond.astype(int), NreadsII
            logfvec_tmp=deepcopy(logfvec)

        Pn_f = np.zeros((len(logfvec_tmp),len(unicounts)))
        
        if method == 'poisson':
            mean_n=Nreads*np.exp(logfvec_tmp)
            Pn_f= super().PoisPar(mean_n,unicounts)

        if method == 'negative_binomial':
            beta_mv = paras[1]
            alpha_mv =paras[2]
            mean_n=Nreads*np.exp(logfvec_tmp)
            var_n=mean_n+beta_mv*np.power(mean_n,alpha_mv)
            Pn_f= super().NegBinParMtr(mean_n,var_n,unicounts)
            
        return np.log(Pn_f)


    def _callback(self, paras):
        '''prints iteration info. called by scipy.minimize'''
        global curr_iter
        print(''.join(['{0:d} ']+['{'+str(it)+':3.6f} ' for it in range(1,len(paras)+1)]).format(*([curr_iter]+list(paras))))
        curr_iter += 1

    def get_log_likelihood_noise(self, paras, f_thresh):    #changed freq_dtype default to float64 
        
        """
        Likelihood of the sampling noise model to compute
        """

        n_lim = 1000

        if len(paras) == 4:
            method = 'negative_binomial'

        elif len(paras) == 2:
            method = 'poisson'
        
        logrhofvec, logfvec= self.get_rhof(paras)
        
        NreadsI = self.NreadsI
        NreadsII = self.NreadsII

        n_thresh_1 = np.floor(f_thresh*NreadsI)
        n_thresh_2 = np.floor(f_thresh*NreadsII)

        nfbins = self.nfbins
        freq_dtype = float

        df_thresh= super().data_thresh_freq_noise(f_thresh)
        sparse_rep_thresh = super().get_sparserep(df_thresh)
        indn1_thresh, indn2_thresh, sparse_rep_counts_thresh, unicountvals_1_thresh, unicountvals_2_thresh, NreadsI_thresh, NreadsII_thresh = sparse_rep_thresh
       
        
        logPn1_f= self.get_logPn_f_noise(paras, f_thresh, 'first')
        logPn2_f= self.get_logPn_f_noise(paras, f_thresh, 'second') 


        #====CONDITION

        #SUM-1: P(n_1<n_tresh_1)
        unicountvals1_cond_SUM_1 = np.arange(n_thresh_1)
        unicountvals2_cond_SUM_1 = np.arange(n_lim)
        logPn1_f_cond_SUM_1 = self.get_logPn_f_condition_noise(paras, f_thresh, 'first', unicountvals1_cond_SUM_1, unicountvals2_cond_SUM_1)
        logPn2_f_cond_SUM_1 = self.get_logPn_f_condition_noise(paras, f_thresh, 'second', unicountvals1_cond_SUM_1, unicountvals2_cond_SUM_1)

        #SUM-2: P(n_2<n_tresh_2)
        unicountvals1_cond_SUM_2 = np.arange(n_lim)
        unicountvals2_cond_SUM_2 = np.arange(n_thresh_2)
        logPn1_f_cond_SUM_2 = self.get_logPn_f_condition_noise(paras, f_thresh, 'first', unicountvals1_cond_SUM_2, unicountvals2_cond_SUM_2)
        logPn2_f_cond_SUM_2 = self.get_logPn_f_condition_noise(paras, f_thresh, 'second', unicountvals1_cond_SUM_2, unicountvals2_cond_SUM_2)

        #SUM-3 : P(n_1<n_tresh_1,n_2<n_thresh_2)
        unicountvals1_cond_SUM_3 = np.arange(n_thresh_1)
        unicountvals2_cond_SUM_3 = np.arange(n_thresh_2)
        logPn1_f_cond_SUM_3 = self.get_logPn_f_condition_noise(paras, f_thresh, 'first', unicountvals1_cond_SUM_3, unicountvals2_cond_SUM_3)
        logPn2_f_cond_SUM_3 = self.get_logPn_f_condition_noise(paras, f_thresh, 'second', unicountvals1_cond_SUM_3, unicountvals2_cond_SUM_3)


      
        dlogfby2=np.diff(logfvec)/2. #for later integration via trapezoid method
        
        #SUM-1: P(n_1<n_tresh_1)  
        Pn1n2_cond_SUM_1=np.zeros((len(unicountvals1_cond_SUM_1), len(unicountvals2_cond_SUM_1))) 
        for i in range (len(unicountvals1_cond_SUM_1)):
            for j in range (len(unicountvals2_cond_SUM_1)):
                integ = np.exp(logrhofvec+logPn2_f_cond_SUM_1[:,j]+logPn1_f_cond_SUM_1[:,i]+ logfvec)
                Pn1n2_cond_SUM_1[i,j] = np.dot(dlogfby2,integ[1:] + integ[:-1])

        #SUM-2: P(n_2<n_tresh_2)  
        Pn1n2_cond_SUM_2=np.zeros((len(unicountvals1_cond_SUM_2), len(unicountvals2_cond_SUM_2))) 
        for i in range (len(unicountvals1_cond_SUM_2)):
            for j in range (len(unicountvals2_cond_SUM_2)):
                integ = np.exp(logrhofvec+logPn2_f_cond_SUM_2[:,j]+logPn1_f_cond_SUM_2[:,i]+ logfvec)
                Pn1n2_cond_SUM_2[i,j] = np.dot(dlogfby2,integ[1:] + integ[:-1])

        #SUM-3 : P(n_1<n_tresh_1,n_2<n_thresh_2)
        Pn1n2_cond_SUM_3=np.zeros((len(unicountvals1_cond_SUM_3), len(unicountvals2_cond_SUM_3))) 
        for i in range (len(unicountvals1_cond_SUM_3)):
            for j in range (len(unicountvals2_cond_SUM_3)):
                integ = np.exp(logrhofvec+logPn2_f_cond_SUM_3[:,j]+logPn1_f_cond_SUM_3[:,i]+ logfvec)
                Pn1n2_cond_SUM_3[i,j] = np.dot(dlogfby2,integ[1:] + integ[:-1])
                
        
        # Computing P(n1,n2)
        Pn1n2=np.zeros((len(sparse_rep_counts_thresh))) 
        for it,(n1_it, n2_it) in enumerate(zip(indn1_thresh,indn2_thresh)):
            integ = np.exp(logrhofvec+logPn2_f[:,n2_it]+logPn1_f[:,n1_it]+ logfvec )
            Pn1n2[it] = np.dot(dlogfby2,integ[1:] + integ[:-1])
            
        N_obs = np.sum(sparse_rep_counts_thresh)

        SUM_1 = np.sum(Pn1n2_cond_SUM_1)
        SUM_2 = np.sum(Pn1n2_cond_SUM_2)
        SUM_3 = np.sum(Pn1n2_cond_SUM_3)

        C = 1 - SUM_1 - SUM_2 + SUM_3 
        print("c = " + str(C))

        #if C >=1:
        #    raise TypeError("P(cond) must be smaller than 1")
        
        Energy = - (1/N_obs)*np.dot(sparse_rep_counts_thresh, np.log(Pn1n2)) + np.log(C)
        
        return N_obs*Energy



    def learning_Null_model(self, init_paras, f_thresh):

        """
        Different uses for this function to explain, explain the use of NreadsItrue and NreadsIItrue
        """
        if len(init_paras) == 4:
            method = 'negative_binomial'
            parameter_labels = ['alph_rho', 'beta', 'alpha', 'f_min']

        elif len(init_paras) == 2:
            method = 'poisson'
            parameter_labels = ['alph_rho', 'f_min']
            
        assert len(parameter_labels) == len(init_paras), "number of model and initial paras differ!"

        partialobjfunc = partial(self.get_log_likelihood_noise, f_thresh = f_thresh)
            
        nullfunctol = 1e-6
        nullmaxiter = 200
        header = ['Iter'] + parameter_labels
        print(''.join(['{' + str(it) + ':9s} ' for it in range(len(init_paras) + 1)]).format(*header))
        
        global curr_iter
        curr_iter = 1
        callbackp = self._callback
            
        #bnds = ((-3, -0.5), (0.001, 10), (0.0001, 10))
            
        outstruct = minimize(partialobjfunc, init_paras,  method= 'SLSQP', callback=callbackp,
                             options={'ftol': nullfunctol, 'disp': True, 'maxiter': nullmaxiter})

        
        return np.hstack((outstruct.x, outstruct.fun))


