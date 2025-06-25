import numpy as np
import scipy.io



# ================================================================================================================ #
# Sets to zero the diagonal of the matrix
#

def dediag(mat):
    return mat-np.diag(np.diag(mat))


# ================================================================================================================ #
# Loads matrix related to EC (pk or CI estimate): it can be the EC matrix, the matrix of p-values, of counts (rela-
#                                                 ted to p-value), of delays, etc. (see the list below)
#

def load_mat(main_path, sim_folder, binsize=0.3, file_prefix='Cult_', folder_prefix='TECN_TE', ECmeas='TE', suffix='pk', mat_type='peakEC'):
    
    '''
         suffix EC       suffix CI       mat_type
      –––––––––––––––––––––––––––––––––––––––––––––––
         -              _CI              ci
        Pk                -              peakEC

        pk_Pcount       _CI_Pcount       count
        pk_Pval         _CI_Pval         P
        pk_sum          _CI_sum          sum
        pk_sum_sq       _CI_sum_sq       sum_sq
        pk_Zscored      _CI_Zscored      Zscored_EC
      ––––––––––––––––––––––––––––––––––––––––––––––– 
      
          suffix                         mat_type
      –––––––––––––––––––––––––––––––––––––––––––––––
        _comm_delays                     delays
        _delay$num                       mat_del
        _DistanceMat                     mat_d          # only in path_results_TE
      –––––––––––––––––––––––––––––––––––––––––––––––
    
    '''
    path  = main_path+sim_folder+folder_prefix+'_Prog_Length_binsize'+str(binsize)+'/'
    mat   = scipy.io.loadmat(path+file_prefix+ECmeas+suffix+'.mat')
    
    return dediag(mat[mat_type]);



# ================================================================================================================ #
# Load EC, EC significative (jittering test), EC zscored,  EC significative zscored
#

def load_complete_measures(main_path, sim_folder, binsize=0.3, alpha=0.001, alpha_max=0.998, file_prefix='Cult_', folder_prefix='TECN_', ECmeas = 'TE', tp = 'pk'):
    
    if tp=='pk':
        M     = load_mat(main_path, sim_folder, binsize=binsize, file_prefix=file_prefix, folder_prefix=folder_prefix+ECmeas, ECmeas=ECmeas, suffix='Pk',           mat_type='peakEC')
    else:
        M     = load_mat(main_path, sim_folder, binsize=binsize, file_prefix=file_prefix, folder_prefix=folder_prefix+ECmeas, ECmeas=ECmeas, suffix=tp,             mat_type='ci')    
    M[np.isnan(M)] = 0
    M[np.isinf(M)] = 0
    
    
    M_Pval    = load_mat(main_path, sim_folder, binsize=binsize, file_prefix=file_prefix, folder_prefix=folder_prefix+ECmeas, ECmeas=ECmeas, suffix=tp+'_Pval',     mat_type='P')   
    M_Pval[np.isnan(M_Pval)] = 0
    M_Pval[np.isinf(M_Pval)] = 0
    
    M_Zscored = load_mat(main_path, sim_folder, binsize=binsize, file_prefix=file_prefix, folder_prefix=folder_prefix+ECmeas, ECmeas=ECmeas, suffix=tp+'_Zscored',  mat_type='Zscored_EC')
    M_Zscored[np.isnan(M_Zscored)] = 0
    M_Zscored[np.isinf(M_Zscored)] = 0
    '''
    M_Pcount  = load_mat(main_path, sim_folder, binsize=binsize, file_prefix=file_prefix, folder_prefix=folder_prefix+ECmeas, ECmeas=ECmeas, suffix=tp+'_Pcount',   mat_type='count')
    
    M_sum     = load_mat(main_path, sim_folder, binsize=binsize, file_prefix=file_prefix, folder_prefix=folder_prefix+ECmeas, ECmeas=ECmeas, suffix=tp+'_sum',      mat_type='sum')
    
    M_sum_sq  = load_mat(main_path, sim_folder, binsize=binsize, file_prefix=file_prefix, folder_prefix=folder_prefix+ECmeas, ECmeas=ECmeas, suffix=tp+'_sum_sq',   mat_type='sum_sq')
    '''
    
    M_sign                = np.zeros_like(M)
    M_sign[M_Pval<=alpha] = M[M_Pval<=alpha]
    if ECmeas != 'TE' and tp=='pk':
        M_sign[M_Pval>=alpha_max] = M[M_Pval>=alpha_max]
    
    M_Zscored_sign                = np.zeros_like(M)
    M_Zscored_sign[M_Pval<=alpha] = M_Zscored[M_Pval<=alpha]
    if ECmeas != 'TE' and tp=='pk':
        M_Zscored_sign[M_Pval>=alpha_max] = M_Zscored[M_Pval>=alpha_max]
    
    M_Zscored_sign[np.isnan(M_Zscored_sign)] = 0
    M_Zscored_sign[np.isinf(M_Zscored_sign)] = 0
    
    return M, M_sign, M_Zscored, M_Zscored_sign, M_Pval





# ================================================================================================================ #
# Loads TE, returns significant TE, significant TE normalized with the entropy of the sender channel/neuron,
# the communication delays and the matrix of distances between channel pairs.
#


def load_TEsign(sim_folder,main_path,delay_C=5,binsize=0.3,alpha=0.001,tstart=15,DeltaT=30):
    
    if binsize!=1:
        path_TE   = main_path+sim_folder+'TE_Prog_Length_binsize'+str(binsize)+'/'
        TE_1 = scipy.io.loadmat(path_TE+'Cult_DistanceMat.mat')    
        dist_mat =  dediag(TE_1['mat_d'])
        del TE_1
    else:
        path_TE   = main_path+sim_folder+'TE_Prog_Length/'
        TE_1 = scipy.io.loadmat(path_TE+'Cult_DistanceMat.mat')    
        dist_mat =  dediag(TE_1['mat_d'])
        del TE_1
    
    #TE_1 = scipy.io.loadmat(path_TE+'Cult_Pval.mat')    
    #Pval = dediag(TE_1['P'])
    #del TE_1

    
    TE_1 = scipy.io.loadmat(path_TE+'Cult_TEPk.mat')
    TE = dediag(TE_1['peakTE'])
    #TE_sign[Pval>alpha]=0
    del TE_1
    
    Pcount = np.zeros_like(TE)
    TE_1 = scipy.io.loadmat(path_TE+'Cult_Pcount_1_10.mat')    
    Pcount += dediag(TE_1['count'])
    tot = 10
    for n in range(10,1000,10):
        TE_1 = scipy.io.loadmat(path_TE+'Cult_Pcount_'+str(n)+'_'+str(n+9)+'.mat')    
        Pcount += dediag(TE_1['count'])
        del TE_1
        tot += 10
        
    Pval = Pcount/tot
    
    TE_1 = scipy.io.loadmat(path_TE+'Cult_TEPk.mat')
    TE_sign = dediag(TE_1['peakTE'])
    TE_sign[Pval>alpha]=0
    del TE_1
    
    
    del_1 = scipy.io.loadmat(path_TE+'Cult_TEdelays.mat')    
    TE_delays_sign = dediag(del_1['TEdelays_interi'])
    #TE_delays_sign[Pval>alpha]=0
    del del_1
    
    TE_sign[TE_sign<10**-7.5]=0
    
    
    path_M_coinc = main_path+sim_folder+'/M_coincidences/'
    delay_max = int(delay_C/binsize)
    dim = len(TE_sign)
    time_sec = DeltaT*60
    fs = 1000/binsize
    if binsize!=1:
        C = np.loadtxt(path_M_coinc+'coincidences_del'+str(delay_max)+'_'+str(tstart)+'_'+str(tstart+DeltaT)+'_mins'+str(np.round(binsize,1))+'.txt', dtype=float, delimiter=' ', skiprows=0)
    else:
        C = np.loadtxt(path_M_coinc+'coincidences_del'+str(delay_max)+'_'+str(tstart)+'_'+str(tstart+DeltaT)+'_mins.txt', dtype=float, delimiter=' ', skiprows=0)
    M = np.zeros((dim,dim))
    for i in range(dim):
        for j in range(dim):
            M[i,j] = C[dim*i+j, 0] # no dediag qui, così hai la conta per i rate
    
    rec_entropy_rate = np.zeros((dim,dim))
    for idx in range(dim):
        p = M[idx,idx]/time_sec/fs; 
        rec_entropy_rate[:,idx] = np.ones(dim)*(-p*np.log2(p)-(1-p)*np.log2(1-p))
    TE_norm     = np.divide(TE,rec_entropy_rate)
    TE_sig_norm = np.divide(TE_sign,rec_entropy_rate)

    del rec_entropy_rate, M, p
    
    return TE, TE_sign, TE_norm, TE_sig_norm, Pval, dist_mat


# ================================================================================================================ #
# Loads SC, XCov and Q matrices (non-significant) and thresholds each link with the n_sigma sigmas of the cross-
# correlogram function. Returns the thresholded EC, the normalised (with sender-rate) EC and the delays matrices.
#

def load_CQ(sim_folder,main_path,n_sigmaC=0,n_sigmaQ=0,delay_C=5,binsize=0.3,tstart=15,DeltaT=30):
    
    path_mats = main_path+sim_folder+'/M_mats/'
    delay_max = int(delay_C/binsize)
    
    Cdels = np.loadtxt(path_mats+'C_delays_del'+str(delay_max)+'_'+str(tstart)+'_'+str(tstart+DeltaT)+'_mins'+str(binsize)+'.txt', dtype=float, delimiter=' ', skiprows=0)
    Qdels = np.loadtxt(path_mats+'Q_delays_del'+str(delay_max)+'_'+str(tstart)+'_'+str(tstart+DeltaT)+'_mins'+str(binsize)+'.txt', dtype=float, delimiter=' ', skiprows=0)


    Cmax = np.loadtxt(path_mats+'Cmax_del'+str(delay_max)+'_'+str(tstart)+'_'+str(tstart+DeltaT)+'_mins'+str(binsize)+'.txt', dtype=float, delimiter=' ', skiprows=0)
    Cstd = np.loadtxt(path_mats+'Cstd_del'+str(delay_max)+'_'+str(tstart)+'_'+str(tstart+DeltaT)+'_mins'+str(binsize)+'.txt', dtype=float, delimiter=' ', skiprows=0)
    Cmax[np.abs(Cmax)<Cstd*n_sigmaC]=0
    
    XCov = np.loadtxt(path_mats+'XCov_max_del'+str(delay_max)+'_'+str(tstart)+'_'+str(tstart+DeltaT)+'_mins'+str(binsize)+'.txt', dtype=float, delimiter=' ', skiprows=0)
    XCov_std = np.loadtxt(path_mats+'XCov_std_del'+str(delay_max)+'_'+str(tstart)+'_'+str(tstart+DeltaT)+'_mins'+str(binsize)+'.txt', dtype=float, delimiter=' ', skiprows=0)
    XCov_ave = np.loadtxt(path_mats+'XCov_ave_del'+str(delay_max)+'_'+str(tstart)+'_'+str(tstart+DeltaT)+'_mins'+str(binsize)+'.txt', dtype=float, delimiter=' ', skiprows=0)
    XCov[np.abs(XCov)<XCov_std*n_sigmaQ+XCov_ave]=0
    
    Qmax = np.loadtxt(path_mats+'Qmax_del'+str(delay_max)+'_'+str(tstart)+'_'+str(tstart+DeltaT)+'_mins'+str(binsize)+'.txt', dtype=float, delimiter=' ', skiprows=0)
    Qmax[XCov==0]=0
    
    
    Cdels[Cmax==0] = 0
    Qdels[XCov==0] = 0
    
    
    Cmax_norm = np.loadtxt(path_mats+'norm_Cmax_del'+str(delay_max)+'_'+str(tstart)+'_'+str(tstart+DeltaT)+'_mins'+str(binsize)+'.txt', dtype=float, delimiter=' ', skiprows=0)
    Cstd_norm = np.loadtxt(path_mats+'norm_Cstd_del'+str(delay_max)+'_'+str(tstart)+'_'+str(tstart+DeltaT)+'_mins'+str(binsize)+'.txt', dtype=float, delimiter=' ', skiprows=0)
    Cmax_norm[np.abs(Cmax)<Cstd*n_sigmaC]=0
    
    XCov_norm = np.loadtxt(path_mats+'norm_XCov_max_del'+str(delay_max)+'_'+str(tstart)+'_'+str(tstart+DeltaT)+'_mins'+str(binsize)+'.txt', dtype=float, delimiter=' ', skiprows=0)
    XCov_std_norm = np.loadtxt(path_mats+'norm_XCov_std_del'+str(delay_max)+'_'+str(tstart)+'_'+str(tstart+DeltaT)+'_mins'+str(binsize)+'.txt', dtype=float, delimiter=' ', skiprows=0)
    XCov_norm[np.abs(XCov)<XCov_std*n_sigmaQ+XCov_ave]=0
    
    Qmax_norm = np.loadtxt(path_mats+'norm_Qmax_del'+str(delay_max)+'_'+str(tstart)+'_'+str(tstart+DeltaT)+'_mins'+str(binsize)+'.txt', dtype=float, delimiter=' ', skiprows=0)
    Qmax_norm[XCov==0]=0
    
        
    return Cmax, Cdels, XCov, Qdels, Qmax, Cmax_norm, XCov_norm, Qmax_norm


# ================================================================================================================ #

