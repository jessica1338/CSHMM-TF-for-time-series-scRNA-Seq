import numpy as np
import argparse
from cvxpy import *
#import random
import progressbar
from collections import defaultdict
from scipy.stats import spearmanr
import time
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.cluster import KMeans,SpectralClustering
from scipy.spatial import distance
from numpy import inf
import pickle
#import multiprocessing as mp
import sys
from scipy.special import erf
from sklearn.linear_model import LogisticRegression
from scipy.stats import ttest_ind, spearmanr,binom_test
import rpy2
import rpy2.robjects as robjects
import pdb
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri

#from statsmodels.sandbox.stats.multicomp import multipletests
import pkg_resources

numpy2ri.activate()


#pip install rpy2==2.8
genlasso = importr('genlasso')
stats = importr('stats')
base = importr('base')
#Matrix = importr('Matrix')
robjects.r("set.seed(1)")


def load_data_2(file_name,max_gene):
    print 'loading data......'
    lines=open(file_name).readlines()
    #print lines
    head=''
    cell_names=lines[0].replace('\n','').split('\t')[1:-1]
    cell_times=np.array(map(float,map(float,lines[1].replace('\n','').split('\t')[1:-1])))
    cell_labels=[]
    gene_exps=[]
    gene_names=[]
    for i,name in enumerate(cell_names):
        splits=name.split('_')
        cell_lab = splits[1]
        if 'Day' in cell_lab:
            cell_lab='NA'
        cell_labels.append(cell_lab)
    for line in lines[2:]:
        line=line.replace('\n','')
        splits=line.split('\t')[:-1]
        gene_names.append(splits[0])
        gene_exp=map(float,splits[1:])
        gene_exps.append(gene_exp)
    cell_exps=np.transpose(np.array(gene_exps))
    gene_names=np.array(gene_names)
    rm_col=np.all(cell_exps<0.1,axis=0)#remove all < 0.1 genes
    n_cell,n_gene = cell_exps.shape
    for j in range(n_gene):
        if np.count_nonzero(cell_exps[:,j])<n_cell/4: #remove the gene that express in less than 25% of cells
            rm_col[j]=True
    cell_exps=cell_exps[:,~rm_col]
    gene_names=gene_names[~rm_col]
    cell_exps=np.log2(cell_exps+1)
    n_cell,n_gene = cell_exps.shape
    if n_gene>max_gene: #select the top [max_gene] genes by variance
        cell_exps_var=np.var(cell_exps,axis=0)
        sort_index = np.argsort(-cell_exps_var)
        select_gene=sort_index[:max_gene]
        cell_exps=cell_exps[:,select_gene]
        gene_names=gene_names[select_gene]
    n_cell,n_gene = cell_exps.shape
    print n_cell, ' cell loaded with ',n_gene,' selected'
    print np.unique(np.array(cell_labels),return_counts=True)
    return cell_names,cell_times,cell_labels,cell_exps,gene_names
def load_data(file_name,max_gene):
    print 'loading data......'
    lines=open(file_name).readlines()
    head=lines[0].replace('\n','')
    cell_names=[]
    cell_day=[]
    cell_labels=[]
    cell_exps=[]
    gene_names=np.array(head.split('\t')[3:])
    for line in lines[1:]:
        line=line.replace('\n','')
        splits=line.split('\t')
        cell_name=splits[0]
        day=float(splits[1])
        label=splits[2]
        gene_exp=splits[3:]
        cell_names.append(cell_name)
        cell_day.append(day)
        cell_labels.append(label)
        cell_exps.append(map(float,gene_exp))
    cell_exps=np.array(cell_exps)
    n_cell,n_gene = cell_exps.shape
    if n_gene>max_gene: #select the top [max_gene] genes by variance
        cell_exps_var=np.var(cell_exps,axis=0)
        sort_index = np.argsort(-cell_exps_var)
        select_gene=sort_index[:max_gene]
        cell_exps=cell_exps[:,select_gene]
        gene_names=gene_names[select_gene]
    n_cell,n_gene = cell_exps.shape
    n_cell,n_gene = cell_exps.shape
    print n_cell, ' cell loaded with ',n_gene,' genes selected'

    gene_names=np.array(map(lambda x: x.lower(),gene_names))
    return cell_names,cell_day,cell_labels,cell_exps,gene_names
def init_var_Jun(init_file,cell_names,cell_times,cell_exps,cell_labels):
    print 'initializing parameters and hidden variable with Juns model structure......'
    st_line=open(init_file).readlines()[0].replace('\n','').split('\t')
    c_line=open(init_file).readlines()[1].replace('\n','').split('\t')
    n_path=len(st_line)+1
    n_state=n_path+1
    n_cell,n_gene = cell_exps.shape
    path_info=[]
    adj_mat=np.zeros((n_state,n_state))
    adj_mat[0,1]=1
    for i in range(n_path):
        path_info.append(defaultdict(lambda:[]))
    path_info[0]['Sp_idx']=0
    path_info[0]['level']=0
    for i in range(n_path):
        path_info[i]['Sc_idx']=i+1
        path_info[i]['ID']=i
        alpha = np.zeros(n_gene)
        #alpha[nz_idx]=probs[:,1]
        path_info[i]["alpha"]=alpha
    for line in st_line:
        pa,pb = map(int,line.split(' '))
        path_info[pa]['child_path'].append(pb)
        path_info[pb]['parent_path']=pa
        path_info[pb]['Sp_idx']=path_info[pa]['Sc_idx']
        path_info[pb]['level']=path_info[pa]['level']+1
        adj_mat[path_info[pb]['Sp_idx'],pb+1]=1
    for i in range(n_state):
        adj_sum=np.sum(adj_mat[i])
        if adj_sum>0:
            adj_mat[i,:]/=adj_sum
    g_param=np.zeros((n_state,n_gene))
    sigma_param=np.ones(n_gene)
    K_param=np.random.sample((n_path,n_gene))*K_param_range
    A=adj_mat
    cell_path=np.zeros(n_cell,dtype=int)
    #cell_names=cell_names.tolist()
    print cell_names
    for line in c_line:
        cn,p=line.split(' ')
        p=int(p)
        if cn in cell_names:
            cell_path[cell_names.index(cn)]=p
    cell_time=np.random.sample((n_cell,))
    model={}
    model['g_param']=g_param
    model['sigma_param']=sigma_param
    model['K_param']=K_param
    model['trans_mat']=adj_mat
    model['path_info']=path_info
    if args.TF_file is not None:
        TF_start_time=np.zeros((n_path,len(dTD.keys())))
        #TF_start_time=np.random.sample((n_path,len(dTD.keys())))#initialize with 0~0.5
        print 'TF start time initizlize 0: ',TF_start_time.shape
        model['TF_start_time']=TF_start_time
    gene_start_time=np.zeros((n_path,n_gene))
    model['gene_start_time']=gene_start_time
    hid_var={}
    hid_var['cell_time']=cell_time
    hid_var['cell_ori_time']=cell_times
    hid_var['cell_path']=cell_path
    hid_var['cell_labels']=np.array(cell_labels)
    optimize_w_nz(model,hid_var,cell_exps)
    path_trans_prob=compute_path_trans_log_prob(adj_mat,path_info)
    model['path_trans_prob']=path_trans_prob
    return model,hid_var
def calculate_gene_start_time(model,gene_names):
    TF_start_time=model['TF_start_time']
    #TF_start_time=np.random.sample(TF_start_time.shape)
    #print TF_start_time
    gene_tf_dict=defaultdict(lambda:[])
    for key,val in dTD.items():
        for v in val:
            gene_tf_dict[v.lower()]+=[key.lower()]
    
    n_path=TF_start_time.shape[0]
#     n_gene = len(gene_names)
#     for key,val in dTD.items():
#         print key,val
#         reg_genes=set(map(lambda x:x.lower(),val))&set(gene_names)
#         gene_idxs=np.array(map(lambda x: gene_names.tolist().index(x),reg_genes))
#         #print gene_idxs
    def _min_time(p,gene):
        reg_tfs=gene_tf_dict[gene]
        if len(reg_tfs)==0:
            #print "zero!"
            #model['gene_start_time'][i,j]=0
            return 0
        else:
            #print "min!"
            tf_idxs=np.array(map(lambda x: TF_names.tolist().index(x),reg_tfs))
            #model['gene_start_time'][i,j]=min(TF_start_time[i,tf_idxs])
            return np.min(TF_start_time[i,tf_idxs])

    for i in range(n_path):
        #reg_gene_start_time=TF_start_time[i,gene_idxs]
        #print reg_gene_start_time
        #print 'path: ',i
        #x=map(lambda x:_min_time(i,x),gene_names)
        #print x
        model['gene_start_time'][i,:]=np.array(map(lambda x:_min_time(i,x),gene_names))
        #for j,gene in enumerate(gene_names):
        #    reg_tfs=gene_tf_dict[gene]
        #    if len(reg_tfs)==0:
        #        #print "zero!"
        #        model['gene_start_time'][i,j]=0
        #    else:
        #        #print "min!"
        #        tf_idxs=np.array(map(lambda x: TF_names.tolist().index(x),reg_tfs))
        #        model['gene_start_time'][i,j]=min(TF_start_time[i,tf_idxs])
        #        #print TF_start_time[i,tf_idxs]
        #        #print model['gene_start_time'][i,j]
def save_model(file_name,model,hid_var):
    print 'saving model to file: ',file_name
    with open(file_name, 'wb') as handle:
        out_dict={}
        out_dict['model']=model
        out_dict['hid_var']=hid_var
        pickle.dump(out_dict, handle)
def load_model(file_name):
    print 'loading model from file: ',file_name
    with open(file_name, 'rb') as handle:
        out_dict = pickle.load(handle)
    return out_dict['model'],out_dict['hid_var']
def optimize_w_nz(model,hid_var,cell_exps):
    print 'M-step: optimizing w param......'
    path_info=model['path_info']
    cell_path=hid_var['cell_path']
    cell_time=hid_var['cell_time']
    K_param=model['K_param']
    g_param=model['g_param']
    n_cell,n_gene=cell_exps.shape
    n_state=g_param.shape[0]
    n_path=n_state-1
    sigma_param=model['sigma_param']
    w_nz=np.ones((n_path,n_gene))  #non zero ratio for each gene in each path (wpj)
    if optimize_w:
        w_nz=np.zeros((n_path,n_gene))  #non zero ratio for each gene in each path (wpj)
        if progress_bar:
            bar = progressbar.ProgressBar(maxval=n_path*n_gene, \
                    widgets=[' [', progressbar.Timer(), '] ',progressbar.Bar('=','[',']'),' ',progressbar.Percentage(),' (', progressbar.ETA(), ') '] )
            bar.start()
        w_split=n_split
        path_gene_w_table=np.zeros((n_path,n_gene,w_split))
        for p in range(n_path):
            Sp_idx=path_info[p]['Sp_idx']
            Sc_idx=path_info[p]['Sc_idx']
            g_a=g_param[Sp_idx]
            g_b=g_param[Sc_idx]
            p_idx=(cell_path==p)
            cell_exps_p=cell_exps[p_idx]
            cell_time_p=cell_time[p_idx]
            for j in range(n_gene):
                x_js=cell_exps_p[:,j]
                mu_x_js=g_b[j]+(g_a[j]-g_b[j])*np.exp(-K_param[p,j]*cell_time_p)
                tmp=(x_js-mu_x_js)**2./(2.*sigma_param[j]**2.)
                prob2 = np.where(x_js!=0.,0.,drop_out_param)
                prob1= np.exp(-tmp)/(sigma_param[j])/np.sqrt(2.*np.pi)
                for ws in range(1,w_split+1):
                    w=1/float(w_split)*ws
                    mix_prob=w*prob1+(1-w)*prob2
                    sum_log_prob=np.sum(np.log(mix_prob))
                    path_gene_w_table[p,j,ws-1]=sum_log_prob
                max_ws=np.argmax(path_gene_w_table[p,j,:])+1
                max_w=1/float(w_split)*max_ws
                w_nz[p,j]=max_w
                #print 'max_w: ',max_w
                if progress_bar:
                    bar.update(p*n_gene+j+1)
        if progress_bar:
            bar.finish()
    model['w_nz']=w_nz
def compute_path_trans_log_prob(trans_mat,path_info):
    ret=[]
    for i,p in enumerate(path_info):
        mult = 1
        now=p
        while(True):
            Sp=now['Sp_idx']
            Sc=now['Sc_idx']
            if Sp==0:
                break
            mult*=trans_mat[Sp,Sc]
            now=path_info[now['parent_path']]
        ret.append(mult)
    with np.errstate(divide='ignore'):
        ret=np.log(np.array(ret))
    return ret
def calc_cell_exp_prob(p,t,model,x_i):
    path_info=model['path_info']
    path_trans_prob=model['path_trans_prob']
    g_param=model['g_param']
    sigma_param=model['sigma_param']
    K_param=model['K_param']
    w_nz=model['w_nz']
    
    gene_start_time=model['gene_start_time']
    
    Sp_idx=path_info[p]['Sp_idx']
    Sc_idx=path_info[p]['Sc_idx']
    g_a=g_param[Sp_idx]
    g_b=g_param[Sc_idx]
    #mu_x_i=g_b+(g_a-g_b)*np.exp(-K_param[p]*t)
    mu_x_i=g_b+(g_a-g_b)*np.exp(-K_param[p]*np.maximum(t-gene_start_time[p],0))
    tmp=(x_i-mu_x_i)**2./(2.*sigma_param**2.)+np.log((sigma_param*np.sqrt(2.*np.pi)) )
    prob2 = np.where(x_i!=0.,0.,drop_out_param)
    mix_prob=w_nz[p]*np.exp(-tmp)+(1-w_nz[p])*prob2
    log_mix_prob=np.log(mix_prob)
    ret=np.sum(log_mix_prob)+path_trans_prob[p]
    return ret
def log_likelihood(model,hid_var,cell_exps):
    ret=0.
    path_info=model['path_info']
    cell_path=hid_var['cell_path']
    cell_time=hid_var['cell_time']
    g_param=model['g_param']
    K_param=model['K_param']
    n_state,n_gene = g_param.shape
    n_path=n_state-1
    n_cell=cell_exps.shape[0]
    for i in range(n_path):
        s_a=path_info[i]['Sp_idx']
        s_b=path_info[i]['Sc_idx']
        delta_g=g_param[s_a]-g_param[s_b]
        ret+=-lamb*np.sum(np.fabs(delta_g))
    for i in range(n_cell):
        p=cell_path[i]
        t=cell_time[i]
        x_i=cell_exps[i,:]
        ret+=calc_cell_exp_prob(p,t,model,x_i)
    return ret
def model_score(model,hid_var,cell_exps,method):
    print 'calculating ',method,' score......'
    n_cell=cell_exps.shape[0]
    g_param=model['g_param']
    n_state,n_gene = g_param.shape
    k=n_gene * n_state * 3 - n_gene # g_param: G*S, K_param: G*P = G*(S-1), sigma_param: G, w_nz: G*(S-1)
    if not optimize_w:
        k=n_gene * n_state * 2 # g_param: G*S, K_param: G*P = G*(S-1), sigma_param: G, w_nz: G*(S-1)
    ll2= 2*log_likelihood(model,hid_var,cell_exps) 
    BIC_score =  ll2 - np.log(n_cell)*k
    AIC_score =  ll2 - 2*k
    GIC2_score = ll2 - k**(1/3.)*k
    GIC3_score = ll2 - 2*np.log(k)*k
    GIC4_score = ll2 - 2*(np.log(k)+np.log(np.log(k)))*k
    GIC5_score = ll2 - np.log(np.log(n_cell))*np.log(k)*k
    GIC6_score = ll2 - np.log(n_cell)*np.log(k)*k
    if method=='BIC':
        return 'BIC= ', BIC_score
    if method=='AIC':
        return 'AIC= ', AIC_score
    if method=='ALL':
        return '(AIC,BIC,G2,G3,G4,G5,G6)=', (AIC_score,BIC_score,GIC2_score,GIC3_score,GIC4_score,GIC5_score,GIC6_score)

def optimize_transition_prob(model,hid_var):
    trans_mat=model['trans_mat']
    path_info=model['path_info']
    cell_path=hid_var['cell_path']
    new_trans_mat=np.zeros(trans_mat.shape)
    for i in range(cell_path.shape[0]):
        p=cell_path[i]
        Sp=path_info[p]['Sp_idx']
        Sc=path_info[p]['Sc_idx']
        new_trans_mat[Sp,Sc]+=1
    sum_vector=np.sum(new_trans_mat,axis=1)
    for i in range(new_trans_mat.shape[0]):
        if sum_vector[i]>0:
            new_trans_mat[i,:]/=sum_vector[i]
    model['trans_mat']=new_trans_mat
    path_trans_prob=compute_path_trans_log_prob(new_trans_mat,path_info)
    model['path_trans_prob']=path_trans_prob
    return
def assign_path_and_time(model,hid_var,cell_exps):
    print 'E-step: assigning new path and time for cell......'
    n_path=model['K_param'].shape[0]
    n_cell=cell_exps.shape[0]
    if progress_bar:
        bar = progressbar.ProgressBar(maxval=n_cell, \
                widgets=[' [', progressbar.Timer(), '] ',progressbar.Bar('=','[',']'),' ',progressbar.Percentage(),' (', progressbar.ETA(), ') '] )
        bar.start()
    time_split=n_split
    path_time_table=np.zeros((n_path,time_split+1))
    cell_path=hid_var['cell_path']
    cell_time=hid_var['cell_time']
    cell_ori_time=hid_var['cell_ori_time']
    if n_anchor:
        anchor=defaultdict(lambda:defaultdict(lambda:[]))
        for i in range(n_cell):
            p=cell_path[i]
            t=cell_time[i]
            prob=calc_cell_exp_prob(p,t,model,cell_exps[i,:])
            anchor[p]['cell_index'].append(i)
            anchor[p]['cell_prob'].append(prob)
        anchor_cell=np.array([-1])
        for p in range(n_path):
            cell_index=np.array(anchor[p]['cell_index'])
            cell_prob=np.array(anchor[p]['cell_prob'])
            anchor_p= cell_index[np.argsort(-cell_prob)[:n_anchor]]
            anchor_cell=np.union1d(anchor_cell,anchor_p)
        print 'anchor cell: ', anchor_cell
    for i in range(n_cell):
        if progress_bar:
            bar.update(i+1)
        if n_anchor and i in anchor_cell:
            continue
        for p in range(n_path):
            for t_sp in range(time_split+1):
                t=t_sp/float(time_split)
                path_time_table[p,t_sp]=calc_cell_exp_prob(p,t,model,cell_exps[i,:])
        max_time=np.argmax(path_time_table,axis=1) #max_time for every path
        max_prob=np.max(path_time_table,axis=1) #prob of every path with max_time 
        new_path = np.argmax(max_prob)
        ori_prob= np.exp(max_prob-np.max(max_prob))
        norm_prob=ori_prob/np.sum(ori_prob)
        valid_idx=np.array(range(n_path))
        sample_prob=norm_prob[valid_idx]/np.sum(norm_prob[valid_idx])
        sample = np.random.multinomial(1,sample_prob)
        for index,s in enumerate(sample):
            if s==1:
                sampled_path=valid_idx[index]
                break
        new_path=valid_idx[np.argmax(max_prob[valid_idx])]
        new_time=max_time[new_path]/float(time_split)
        cell_path[i]=new_path
        if assign_by_prob_sampling:
            cell_path[i]=sampled_path
        cell_time[i]=new_time
    hid_var['cell_time']=cell_time
    hid_var['cell_path']=cell_path
    if progress_bar:
        bar.finish()
    return    
def optimize_K_param(model,hid_var,cell_exps):
    print 'M-step: optimizing K param......'
    K_param=model['K_param']
    new_K_param=np.zeros(K_param.shape)
    w_nz=model['w_nz']
    n_path,n_gene=K_param.shape
    cell_path=hid_var['cell_path']
    cell_time=hid_var['cell_time']
    path_info=model['path_info']
    g_param=model['g_param']
    sigma_param=model['sigma_param']
    gene_start_time=model['gene_start_time']

    k_split=n_split
    path_gene_k_table=np.zeros((n_path,n_gene,k_split))
    if progress_bar:
        bar = progressbar.ProgressBar(maxval=n_path*n_gene, \
                widgets=[' [', progressbar.Timer(), '] ',progressbar.Bar('=','[',']'),' ',progressbar.Percentage(),' (', progressbar.ETA(), ') '] )
        bar.start()
    count=0
    for p in range(n_path):
        Sp_idx=path_info[p]['Sp_idx']
        Sc_idx=path_info[p]['Sc_idx']
        g_a=g_param[Sp_idx]
        g_b=g_param[Sc_idx]
        p_idx=(cell_path==p)
        cell_exps_p=cell_exps[p_idx]
        cell_time_p=cell_time[p_idx]
        for j in range(n_gene):
            x_js=cell_exps_p[:,j]
            for ks in range(1,k_split+1):
                k=K_param_range/float(k_split)*ks
                mu_x_js=g_b[j]+(g_a[j]-g_b[j])*np.exp(-k*np.maximum(cell_time_p-gene_start_time[p,j],0))
                tmp=((x_js-mu_x_js)**2./(2.*sigma_param[j]**2.)+np.log((sigma_param[j]*np.sqrt(2.*np.pi)) ))
                prob2 = np.where(x_js!=0.,0.,drop_out_param)
                mix_prob=w_nz[p,j]*np.exp(-tmp)+(1-w_nz[p,j])*prob2
                sum_log_prob=np.sum(np.log(mix_prob))
                path_gene_k_table[p,j,ks-1]=sum_log_prob
            max_ks=np.argmax(path_gene_k_table[p,j,:])+1
            max_k=K_param_range/float(k_split)*max_ks
            K_param[p,j]=max_k
            count+=1
            if progress_bar:
                bar.update(count)
    if progress_bar:
        bar.finish()
def optimize_sigma_param(model,hid_var,cell_exps):
    print 'M-step: optimizing sigma param......'
    cell_path=hid_var['cell_path']
    cell_time=hid_var['cell_time']
    path_info=model['path_info']
    g_param=model['g_param']
    n_cell,n_gene=cell_exps.shape
    new_sigma_param=np.zeros(n_gene)
    K_param=model['K_param']
    gene_start_time=model['gene_start_time']

    if progress_bar:
        bar = progressbar.ProgressBar(maxval=n_cell, \
                widgets=[' [', progressbar.Timer(), '] ',progressbar.Bar('=','[',']'),' ',progressbar.Percentage(),' (', progressbar.ETA(), ') '] )
        bar.start()
    for i in range(n_cell):
        p=cell_path[i]
        Sp_idx=path_info[p]['Sp_idx']
        Sc_idx=path_info[p]['Sc_idx']
        g_a=g_param[Sp_idx]
        g_b=g_param[Sc_idx]
        t=cell_time[i]
        x_i=cell_exps[i]
        mu_x_i=g_b+(g_a-g_b)*np.exp(-K_param[p]*np.maximum(t-gene_start_time[p],0))
        new_sigma_param+=(x_i-mu_x_i)**2
        if progress_bar:
            bar.update(i)
    new_sigma_param=(new_sigma_param/float(n_cell))**0.5
    new_sigma_param=np.where(new_sigma_param<1,1,new_sigma_param)
    model['sigma_param']=new_sigma_param    
    if progress_bar:
        bar.finish()
def solve_genlasso(y, X, D, lam):
    #out = genlasso.genlasso(y, X=X, D=Matrix.Matrix(D, sparse = True), minlam=1,verbose=False)
    out = genlasso.genlasso(y, X=X, D=D, minlam=lam,verbose=False)
    #return np.zeros(D.shape[1])
    return np.array(stats.coef(out, lam)[0])
    
def optimize_g_param(model,hid_var,cell_exps,gene_names,opt_method="genlasso"):#,alpha=0.3):
    print 'M-step: optimizing g param with genlasso......'
    path_info=model['path_info']
    cell_path=hid_var['cell_path']
    cell_time=hid_var['cell_time']
    K_param=model['K_param']
    g_param=model['g_param']
    n_cell,n_gene=cell_exps.shape
    n_state=g_param.shape[0]
    sigma_param=model['sigma_param']
    gene_start_time=model['gene_start_time']

    A2=np.zeros((n_state-1,n_state))
    _,path_count=np.unique(hid_var['cell_path'],return_counts=True)
    path_nz_g=check_diff_gene(model)
    de=np.ones((n_state-1,n_gene))
    for index in range(n_state-1):
        path=path_info[index]
        Sp_idx=path['Sp_idx']
        Sc_idx=path['Sc_idx']
        target_list=path["TF_targets"]
        #print gene_names[j]
        #print target_list
        #de=1
        de[index,:]+=path["alpha"]
#         for j in range(n_gene):
#             if gene_names[j] in target_list:
#                 #print j,'matched!'
#                 de[index,j]+=alpha
        A2[index,Sp_idx]=1
        A2[index,Sc_idx]=-1
        if lamb_data_mult=='N':
            A2[index,:]*=path_count[index] # multiply by N
        if lamb_data_mult=='sqrtN':
            A2[index,:]*=np.sqrt(path_count[index]) # multiply by sqrt(N)
        if lamb_data_mult=='logN':
            A2[index,:]*=np.log(path_count[index]) # multiply by log(N)
        if lamb_ratio_mult=='sqrtR':
            A2[index,:]*=np.sqrt(path_nz_g[index])
        if lamb_ratio_mult=='R':
            A2[index,:]*=path_nz_g[index]
    A2*=lamb
    #print 'de: ',de
    #print "max de: ", np.max(de)
    #print "min de: ", np.min(de)
    #print "sum de: ", np.sum(de)
    if progress_bar:
        bar = progressbar.ProgressBar(maxval=n_gene, \
                widgets=[' [', progressbar.Timer(), '] ',progressbar.Bar('=','[',']'),' ',progressbar.Percentage(),' (', progressbar.ETA(), ') '] )
        bar.start()
    for j in range(n_gene):
        A1=np.zeros((n_cell,n_state))
        Xjs=np.zeros(n_cell)
        sigma_j=sigma_param[j]
        #de=np.ones(n_state-1)
        for i in range(n_cell):
            p=cell_path[i]
            Sp_idx=path_info[p]['Sp_idx']
            Sc_idx=path_info[p]['Sc_idx']
            t=cell_time[i]
            x_ij=cell_exps[i,j]
            w_ij=np.exp(-K_param[p,j]*np.maximum(t-gene_start_time[p,j],0))
            A1[i,Sp_idx]=w_ij
            A1[i,Sc_idx]=1-w_ij
            Xjs[i]=x_ij
        D=A2/de[:,j][:,None]
        lam=lamb
        y=Xjs/sigma_j
        X=A1/sigma_j
        if opt_method=="genlasso":
            g_param[:,j]=solve_genlasso(y, X, D, lam).flatten()
        elif opt_method=='cvxpy':
            g_js = Variable(n_state)
            objective = Minimize(0.5*sum_squares((A1*g_js-Xjs)/sigma_j)+pnorm(A2*g_js,1))
            constraints=[g_js>=0]
            prob = Problem(objective, constraints)
            result = prob.solve(solver=SCS)
            g_param[:,j]=g_js.value.flatten()
        if progress_bar:
            bar.update(j)
    g_param[g_param<0]=0
    if progress_bar:
        bar.finish()
def check_diff_gene(model):
    path_info=model['path_info']
    g_param=model['g_param']
    n_gene=g_param.shape[1]
    path_nz_g={}
    for index,path in enumerate(path_info):
        Sp_idx=path['Sp_idx']
        Sc_idx=path['Sc_idx']
        g_a=g_param[Sp_idx]
        g_b=g_param[Sc_idx]
        g_abs_diff=np.fabs(g_a-g_b)
        g_abs_diff_nz=np.where(g_abs_diff<1e-1,0,g_abs_diff)
        nz_count= len(g_abs_diff_nz[np.nonzero(g_abs_diff_nz)])
        nz_ratio = nz_count/float(n_gene)
        if verbose:
            print 'path: ', index, ' nz_ratio: ', nz_ratio
        path_nz_g[index]=nz_ratio
    return path_nz_g

def show_cell_time(hid_var):
    cell_path = hid_var['cell_path']
    cell_labels = hid_var['cell_labels']
    cell_time=hid_var['cell_time']
    cell_ori_time=np.array(hid_var['cell_ori_time'])
    paths=np.unique(cell_path)
    for p in paths:
        print '----------path: ',p,'-------------'
        print 'time\tlabel\t assigned_time'
        p_idx=(cell_path==p)
        cell_time_p=np.around(cell_time[p_idx],decimals=2)
        cell_labels_p=cell_labels[p_idx]
        cell_ori_time_p=cell_ori_time[p_idx]
        sort_idx=np.argsort(cell_time_p)
        for lab,ori_t,t in zip(cell_labels_p[sort_idx],cell_ori_time_p[sort_idx], cell_time_p[sort_idx]):
            print ori_t,'\t', lab, '     \t' , t
def compute_ARI_confuss_mat(hid_var,n_path):
    cell_path = hid_var['cell_path']
    cell_labels = hid_var['cell_labels']
    unique_paths=np.unique(cell_path)
    unique_labels={}
    cell_label_num=[]
    ARI_ans=[]
    ARI_pred=[]
    head=[]
    for i,lab in enumerate(cell_labels):
        if lab in unique_labels.keys():
            cell_label_num.append(unique_labels[lab])
        else:
            ID=len(unique_labels.keys())
            unique_labels[lab]=ID
            cell_label_num.append(unique_labels[lab])
            head.append(lab)
        if lab!='NA':
            ARI_ans.append(unique_labels[lab])
            ARI_pred.append(cell_path[i])
    confuss_mat=np.zeros((n_path,len(unique_labels)))
    for i,num in enumerate(cell_label_num):
        confuss_mat[cell_path[i],num]+=1
    print 'confussion matrix:'
    print head
    print confuss_mat
    ARI= adjusted_rand_score(ARI_ans, ARI_pred)
    print 'ARI: ',ARI
    return confuss_mat,ARI
def load_adj_mat(file_name):
    return np.load(file_name)

def path_distance(pa,pb,cell_exps,cell_path):
    pa_center = np.average(cell_exps[cell_path==pa],axis=0)
    pb_center = np.average(cell_exps[cell_path==pb],axis=0)
    return 1-spearmanr(pa_center,pb_center)[0]


def adjust_model_structure(model,hid_var,cell_exps):
    print 'adjusting model structure '
    path_info=model['path_info']
    cell_path=hid_var['cell_path']
    cell_time=hid_var['cell_time']
    n_path=len(path_info)
    valid_parent=sorted(list(set(np.unique(hid_var['cell_path']).tolist())))
    level_path=defaultdict(lambda:[])
    for i,p in enumerate(path_info):
        p['child_path']=[]
        if i not in valid_parent:
            print 'zero path: ', i
            continue
        level_path[p['level']].append(i)
    def getKey(item):
        return item['level']
    childs = sorted([x for x in path_info], key = getKey)
    for p in childs:
        ID=p['ID']
        if ID not in valid_parent or ID ==0:
            continue
        level= p['level']
        new_parent_list=level_path[level-1]
        if len(new_parent_list)==1:
            new_parent = new_parent_list[0]
            print str(new_parent) + ' -> '+str(ID)
        else:
            par_distance = []
            for par in new_parent_list:
                p['Sp_idx']=path_info[par]['Sc_idx']
                p_idx=(cell_path==ID)
                cell_exps_p=cell_exps[p_idx]
                cell_time_p=cell_time[p_idx]
                s=0
                for i in range(cell_exps_p.shape[0]):
                    s+=calc_cell_exp_prob(ID,cell_time_p[i],model,cell_exps_p[i,:])
                par_distance.append(s)
            #par_distance = [path_distance(ID,x,cell_exps,cell_path) for x in new_parent_list]
            new_parent = new_parent_list[np.argmax(par_distance)]
            print str(new_parent) + ' -> '+str(ID)
        p['parent_path']=new_parent
        new_p=path_info[new_parent]
        p['Sp_idx']=new_p['Sc_idx']
        new_p['child_path'].append(ID)
        p['level']=new_p['level']+1

def optimize_likelihood(cell_exps, gene_names, model, hid_var, model_name,store_model=True):
    #assign_path_TF(dTD,model,hid_var,gene_names,cell_exps,assign_by_K=False)
    if store_model:
        save_model(model_name+'_it0.pickle',model, hid_var)
        #model,hid_var=load_model(model_name)
    
    for out_it in range(1,n_iteration+1):
        prev_path=np.array(hid_var['cell_path'],copy=True)
        print 'training iteration: ', out_it
        
        sys.stdout.flush()
        
        print 'cell paths: ',np.unique(hid_var['cell_path'],return_counts=True)

        score = model_score(model,hid_var,cell_exps,method='ALL')
        print 'model score: ',score
        
        optimize_g_param(model,hid_var,cell_exps,gene_names,opt_method=args.opt_method)
        print 'after M-step g_param full log-likelihood: ',log_likelihood(model,hid_var,cell_exps)
       
        sys.stdout.flush()
        #check_diff_gene(model)
        
        optimize_sigma_param(model,hid_var,cell_exps)
        print 'after M-step sigma_param full log-likelihood: ',log_likelihood(model,hid_var,cell_exps)
        optimize_K_param(model,hid_var,cell_exps)
        print 'after M-step K_param full log-likelihood: ',log_likelihood(model,hid_var,cell_exps)
        
        #assign_path_TF(dTD,model,hid_var,gene_names,cell_exps,assign_by_K=False)
        
        #assign_path_and_time(model,hid_var,cell_exps)
        assign_path_and_time(model,hid_var,cell_exps)
        print 'after E-step full log-likelihood: ',log_likelihood(model,hid_var,cell_exps)
        
        
        if args.TF_file is not None:
            #assign_path_TF(dTD,model,hid_var,gene_names,cell_exps,assign_by_K=False)
            assign_path_TF(dTD,model,hid_var,gene_names,cell_exps)
            set_alpha_logistic_regression(model,gene_tf_table)
            calculate_gene_start_time(model,gene_names)
        
        adjust_model_structure(model,hid_var,cell_exps)
        
        sys.stdout.flush()
        optimize_w_nz(model,hid_var,cell_exps)
        print 'after setting w_nz full log-likelihood: ',log_likelihood(model,hid_var,cell_exps)
        
        optimize_transition_prob(model,hid_var)
        #print 'opt trans mat full log-likelihood: ',log_likelihood(model,hid_var,cell_exps)
        path_trans_prob=compute_path_trans_log_prob(model['trans_mat'],model['path_info'])
        model['path_trans_prob']=path_trans_prob
        print 'after setting trans_prob full log-likelihood: ',log_likelihood(model,hid_var,cell_exps)
        n_path=len(model['path_info'])
        compute_ARI_confuss_mat(hid_var,n_path)
        if verbose:
            show_cell_time(hid_var)
        if np.array_equal(prev_path,hid_var['cell_path']):
            print 'path assignment the same as previous iteration, stop training.'
        #if out_it % 10 ==0:
        
        #print model
        #print hid_var
        model['path_info']=map(dict,model['path_info'])
        if store_model:
            save_model(model_name+'_it'+str(out_it)+'.pickle',model, hid_var)
            #model,hid_var=load_model(model_name)
           
        sys.stdout.flush()
    print 'maximum training iteration reached.'
    #print 'after M-step g_param full log-likelihood: ',log_likelihood(model,hid_var,cell_exps)
    
    #optimize_g_param_close(model,hid_var,cell_exps)
    #print 'after M-step g_param full log-likelihood (close): ',log_likelihood(model,hid_var,cell_exps)
def cv_split_idx(cell_day,n_fold=5):
    cell_day = np.array(cell_day)
    unique_day=np.unique(cell_day)
    n_cell=cell_day.shape[0]
    batch=n_cell/n_fold
    fold_idx=np.zeros(n_cell)
    cell_day_dict={}
    for ud in unique_day:
        ud_idx=np.where(cell_day==ud)[0]
        ud_count=ud_idx.shape[0]
        np.random.shuffle(ud_idx)
        batch=ud_count/n_fold
        for i in range(n_fold):
            fold_idx[ud_idx[batch*i:batch*(i+1)]]=i+1
        if batch*(i+1)<ud_count:
            fold_idx[ud_idx[batch*i:]]=i+1
    return np.array(fold_idx,dtype=int)

def assign_path_TF(dTD,model,hid_var,gene_names,cell_exps,assign_by_K=False,pcut=0.1,cutoff=1,gtop=1000,fold_change=1.5,k_ratio=0.9,ttest_pvs=[0.05],FC_cuts=[0.6,1.0,1.5],marker=[],save=None):
#def assign_path_TF(dTD,model,hid_var,gene_names,cell_exps,pcut=0.1,marker=[],save=None):
    print 'assigning TF to each path......'
    marker_idx = [gene_names.tolist().index(x.lower()) for x in marker]
    print marker
    print marker_idx
    gene_names = np.array(map(lambda x:x.lower(),gene_names))

    k_unique=np.unique(model['K_param'])
    
    cell_path=hid_var['cell_path']
    cell_time=hid_var['cell_time']
    sigma_param=model['sigma_param']
    path_info=model["path_info"] 
    g_param=model["g_param"] 
    n_interval = 5
    pcut=0.05
    ntop=0
    tf_max=10
    print gene_names

    def generate_interval(n_interval, k_unique):
        k_size = k_unique.shape[0]/n_interval
        intervals=[]
        for i in range(n_interval):
            intervals+=[(k_size*i, k_size*(i+1))]
            if k_size*(i+1)>=k_unique.shape[0]:
                break
            intervals+=[(k_size/2+k_size*i, k_size/2+k_size*(i+1))]
        return intervals
    intervals=generate_interval(n_interval, k_unique)
    #n_interval = 10
    #intervals_cont_assign=generate_interval(n_interval, k_unique)
    #cutoff=1
    out_file=open('tf_list.txt','w')
    def get_DE_genes_by_group(exps_a,exps_b):
        pvs=[]
        #res=ttest_ind(exps_a,exps_b)[-1]
        for i in range(cell_exps_p.shape[1]):
            #idx_a=exps_a[:,i]>0
            #idx_b=exps_b[:,i]>0
            pv = ttest_ind(exps_a[:,i],exps_b[:,i])[-1]
            pvs+=[pv]
        return np.array(pvs)
    def get_DE_genes_by_ttest(cell_time_p,cell_exps_p):
        exps_a=cell_exps_p[cell_time_p<0.5]
        exps_b=cell_exps_p[cell_time_p>0.5]
        return get_DE_genes_by_group(exps_a,exps_b)
    def get_sibling_path_idx(p_idx,path_info):
        #print path_info[p].keys()
        #sibling_idx=[]
        for i,p in enumerate(path_info):
            if p_idx in p["child_path"]:
                return np.array(p["child_path"]+[i])       
    def tellDifference(cells_A,cells_B,g_idx,pcut=0.05,fcut=0.6,print_X_Y=False):
        if cells_B.shape[0]==0:
            return [0,1,0]
        X=cells_A[:,g_idx]
        Y=cells_B[:,g_idx]
        pv=ttest_ind(X,Y)[1]
        if print_X_Y:
            print "mean X:",np.sum(X)/X.shape[0]
            print "mean Y:",np.sum(Y)/Y.shape[0]
        fc=np.sum(X)/X.shape[0]-np.sum(Y)/Y.shape[0]
        if pv<pcut and fc > fcut:
            return [1, pv, fc]
        if pv<pcut and fc < -fcut:
            return [-1, pv, fc]
        return [0,pv,fc]
    def assign_eTFs(cell_exps,cell_exps_p,cell_path,sib_idx,gene_names,fcut=0.6,tflist=None):
        if sib_idx is None:
            return None
        tflistpath=pkg_resources.resource_filename(__name__,"HumanTFList.txt") if tflist==None else tflist
        try:
            with open(tflistpath,'r') as f:
                TFs=f.readlines()
                TFs=[item.strip().split()[0] for item in TFs]   
        except:
            print("error! Please check your input TF List file")
            sys.exit(0)
        eTFs=[]
        cell_exps_sibs=cell_exps[cell_path==sib_idx[1]]
        if len(sib_idx)>3:
            for sib in sib_idx[2:-1]:
                print "add ",sib, "as sibling path"
                print "adding matrix shape: ",cell_exps[cell_path==sib].shape
                cell_exps_sibs=np.vstack((cell_exps_sibs,cell_exps[cell_path==sib]))
            print "final shape: ",cell_exps_sibs.shape
        cell_exps_par=cell_exps[cell_path==sib_idx[-1]]
        for gene in gene_names:
            if gene.upper() not in TFs:
                continue
            g_idx = gene_names.tolist().index(gene)
            #print "telldifference with parent: "
            [jflag1,pvp,fcp]=tellDifference(cell_exps_p,cell_exps_par,g_idx,fcut=fcut)
            #print "telldifference with siblings: "
            
            [jflag2,pvs,fcs]=tellDifference(cell_exps_p,cell_exps_sibs,g_idx,fcut=fcut)
            #print [jflag1,pvp,fcp]
            #print [jflag2,pvs,fcs]
            #print "jflag1: ",jflag1
            #print "jflag2: ",jflag2
            #print "len sib_idx: ",len(sib_idx)
            #print " sib_idx: ",sib_idx
            if gene.lower() =="cebpd":
                #Testing
                print [pvp,gene,fcp,fcs]
                print "telldifference with parent: "
                [jflag1,pvp,fcp]=tellDifference(cell_exps_p,cell_exps_par,g_idx,fcut=fcut,print_X_Y=True)
                print "telldifference with siblings: "
                [jflag2,pvs,fcs]=tellDifference(cell_exps_p,cell_exps_sibs,g_idx,fcut=fcut,print_X_Y=True)

            if (jflag1*jflag2>0) or (jflag1!=0 and len(sib_idx)==2):
                #print "yes! ",[pvp,j,fcp,fcs]

                eTFs.append([pvp,gene,fcp,fcs])

            #else:
                #print "no!"
        eTFs.sort()
        return eTFs




    mean_path_g={}
    for i,p in enumerate(model['path_info']):
        out_file.write("path: "+str(i)+'\n')
        print 'path: ',i
        g_pa=model['g_param'][p['Sp_idx']]
        g_pb=model['g_param'][p['Sc_idx']]
        k_p=model['K_param'][p['ID']]
        gene_start_time=model['gene_start_time']
        w_nz=model['w_nz']
       
        p_idx=cell_path==p["ID"]
        cell_exps_p=cell_exps[p_idx]
        mean_path_g[i]=np.mean(cell_exps_p,axis=0)
        cell_time_p=cell_time[p_idx]
        cell_exps_p_g_sum = np.sum(cell_exps_p,axis=0)/cell_exps_p.shape[1]
        
        print "cell_exps_p_g_sum.shape: ",cell_exps_p_g_sum.shape
        sib_idx= get_sibling_path_idx(i,path_info)
        eTFs=[]
        #print "path: ",i, "eTFs:"
        if sib_idx is not None:
        #    print "sib_idx: ",sib_idx
            cell_exps_parent = cell_exps[cell_path==sib_idx[-1]]
            eTFs=assign_eTFs(cell_exps,cell_exps_p,cell_path,sib_idx,gene_names)
        #    print eTFs
        #if sib_idx != None:

        #print cell_exps_p_g_sum
        exp_cutoff=0.01 
        TF_path=set()
        TF_pvalue=defaultdict(lambda:[])
        TF_method=defaultdict(lambda:[])
        '''
        if assign_by_K:
            for itv in intervals:
                idxs=np.array([ k_p==k_unique[x] for x in range(itv[0],itv[1])],dtype=int)
                g_idx = np.array(sum(idxs),dtype=bool)
                g_idx=np.logical_and(g_idx,np.fabs(g_pa-g_pb)>exp_cutoff)
                tf_all,_=getEnrichTF(dTD,gene_names[g_idx],gene_names,pcut=0.1)
                for pv,tf in tf_all:
                    TF_pvalue[tf.lower()]+=[pv]
                #TF_path|=set([x[1] for x in tf_all])
            #print "found TF: ",len(TF_path)
            #TF_path = map(lambda x:x.lower(),TF_path)
            #print TF_path
            #TF_time=[]
        else:
        '''
        marker_dict={}
        for mk in marker:
            marker_dict[mk]=1
        tops=[j*100 for j in range(1,11)]
        
        DE_gene_lists=[]
        DE_gene_method=[]
        DE_gene_top=[]
        
        res1=get_DE_genes_by_ttest(cell_time_p,cell_exps_p)
        #print 'path_ttest: ',res
        #print res.shape
        #print res[res<0.1].shape
        
        #print "path ttest"
        #path_ttest_idx=np.argsort(-res)[:gtop]
        path_ttest_DE_gene=np.argsort(res1)
        #path_ttest_idx = [res<0.1]
        DE_gene_lists+=[path_ttest_DE_gene]
        DE_gene_method+=["path_ttest_top"]
        DE_gene_top+=[tops]
        #DE_gene_lists+=[path_ttest_DE_gene]
        DE_gene_lists+=[res1]
        DE_gene_method+=["path_ttest_cutoff"]
        DE_gene_top+=[[0.05,0.1]]


        res2=get_DE_genes_by_group(cell_exps_p,cell_exps[np.invert(p_idx)])
        #print 'path_all_other: ',res
        path_all_other=np.argsort(res2)
        DE_gene_lists+=[path_all_other]
        DE_gene_method+=["path_all_other_top"]
        DE_gene_top+=[tops]
        #DE_gene_lists+=[path_all_other]
        DE_gene_lists+=[res2]
        DE_gene_method+=["path_all_other_cutoff"]
        DE_gene_top+=[[0.05,0.1]]
            

        if sib_idx is not None:
            res3=get_DE_genes_by_group(cell_exps_p,cell_exps[cell_path==sib_idx[-1]])
            #print 'path_parent: ',res
            path_parent=np.argsort(res3)
            DE_gene_lists+=[path_parent]
            DE_gene_method+=["path_parent_top"]
            DE_gene_top+=[tops]
            #DE_gene_lists+=[path_parent]
            DE_gene_lists+=[res3]
            DE_gene_method+=["path_parent_cutoff"]
            DE_gene_top+=[[0.05,0.1]]
        ##path_ttest_tf = getEnrichTF(dTD,gene_names[path_ttest_idx],gene_names,pcut=0.1,marker=marker_dict)

        
        #print "g_diff"
        g_diff=np.fabs(g_pa-g_pb)

        g_diff_idx = g_diff > cutoff
        #g_diff_tf = getEnrichTF(dTD,gene_names[g_diff_idx],gene_names,pcut=0.1,marker=marker_dict)

        #g_diff_tf = getEnrichTF(dTD,gene_names[g_diff_idx],gene_names,pcut=0.1,marker=marker_dict)
        #g_idx=g_fold_idx
        #print "g_diff_top"
        #g_top_idx=np.argsort(-g_diff)[:gtop]
        g_top_DE_gene=np.argsort(-g_diff)
        #g_top_tf = getEnrichTF(dTD,gene_names[g_top_idx],gene_names,pcut=0.1,marker=marker_dict)
        
        DE_gene_lists+=[g_top_DE_gene]
        DE_gene_method+=["g_diff_top"]
        DE_gene_top+=[tops]
        #DE_gene_lists+=[-g_top_DE_gene]
        DE_gene_lists+=[-g_diff]
        DE_gene_method+=["g_diff_cutoff"]
        DE_gene_top+=[[-0.6,-1.0,-1.5]]
        
        
        '''
        #print "g_fold"
        g_fold = np.fabs(g_diff/(g_pa+1e-5))
        #np.nan_to_num(g_fold)
        #g_idx=np.fabs(g_pa-g_pb)>fold_change
        g_fold_idx=np.logical_and(g_fold>fold_change,g_pa>exp_cutoff)
        g_fold_tf = getEnrichTF(dTD,gene_names[g_fold_idx],gene_names,pcut=0.1,marker=marker_dict)
        
        #print "g_fold_top"
        g_fold = np.fabs(g_diff/(g_pa+1e-5))
        #np.nan_to_num(g_fold)
        #g_idx=np.fabs(g_pa-g_pb)>fold_change
        #g_fold_idx=np.logical_and(g_fold>fold_change,g_pa>exp_cutoff)
        g_fold_top_idx=np.argsort(-g_diff)[:gtop]
        g_fold_top_tf = getEnrichTF(dTD,gene_names[g_fold_top_idx],gene_names,pcut=0.1,marker=marker_dict)
        '''

        #intervals
        itv = intervals[-1]
        low = int(k_unique.shape[0]*k_ratio)
        high = k_unique.shape[0]
        idxs=np.array([ k_p==k_unique[x] for x in range(low,high)],dtype=int)
        #print "k_top"
        k_top_idx = np.array(sum(idxs),dtype=bool)
        k_top_idx=np.logical_and(k_top_idx,np.fabs(g_pa-g_pb)>exp_cutoff)
        #k_top_idx=np.logical_and(k_top_idx,cell_exps_p_g_sum>exp_cutoff)
        #k_top_tf=getEnrichTF(dTD,gene_names[k_top_idx],gene_names,pcut=0.1,marker=marker_dict)
        #print k_top_idx
        #print k_top_idx.shape
        DE_gene_lists+=[np.array(range(len(gene_names)))[k_top_idx]]
        DE_gene_method+=["k_top"]
        DE_gene_top+=[[gene_names[k_top_idx].shape[0]]]
        


        #print "sp_cor"
        sp_cor = []
        for j in range(len(gene_names)):
            sp_cor+=[spearmanr(cell_time_p,cell_exps_p[:,j])[0]]
        
        #sp_cor_top_idx=np.argsort(-np.fabs(sp_cor))[:gtop]
        sp_cor_top_idx=np.argsort(-np.fabs(sp_cor))
        #sp_cor_top_tf = getEnrichTF(dTD,gene_names[sp_cor_top_idx],gene_names,pcut=0.1,marker=marker_dict)
        DE_gene_lists+=[sp_cor_top_idx]
        DE_gene_method+=["sp_cor_top"]
        DE_gene_top+=[tops]
        
        #print "g_diff_ratio"
        g_diff_ratio = g_diff/sigma_param
        g_diff_ratio_idx= np.argsort(-np.fabs(g_diff_ratio))
        #g_diff_ratio_idx= np.argsort(-np.fabs(g_diff_ratio))[:gtop]
        #g_diff_ratio_tf = getEnrichTF(dTD,gene_names[g_diff_ratio_idx],gene_names,pcut=0.1,marker=marker_dict)
        
        DE_gene_lists+=[g_diff_ratio_idx]
        DE_gene_method+=["g_diff_ratio_top"]
        DE_gene_top+=[tops]

        #print "parent_path: ", get_sibling_path_idx(i, path_info)
        

        for key,val in marker_dict.items():
            print "marker p-value: ",key, val

        #sibling_idx=get_sibling_path_idx(i,path_info)
        #if sibling_idx is not None:
        #    print np.argmin(g_param[sibling_idx,:],axis = 0)
        
        
        #g_idx=g_fold>cutoff
        #g_diff_idx=srt_idx[:gtop]
        #g_idx = np.logical_and(g_fold_idx,g_diff_idx)
        tf_all=[]
        method_all=[]
        TF_pv_list=[]
        '''
        for DE_genes,DE_method,DE_top in zip( DE_gene_lists,DE_gene_method,DE_gene_top):
            print DE_method , "top/cutoff: ",DE_top
            for top in DE_top:
                if 'top' in DE_method:
                    print "DE_genes_size: ",top
                    tfs,_=getEnrichTF(dTD,gene_names[DE_genes[:top]],gene_names,pcut=0.1,marker=marker_dict)
                if 'cutoff' in DE_method:
                    print "DE_genes_size: ",gene_names[DE_genes<top].shape
                    tfs,_=getEnrichTF(dTD,gene_names[DE_genes<top],gene_names,pcut=0.1,marker=marker_dict)
                TF_pv_list+=_
                print 'tf size: ',len(tfs)
                tf_all+=tfs
                method_all+=[DE_method +' '+ str(top)]*len(tfs)
        '''
        fc_no_map={}
        cnt=0
        if sib_idx is not None:
            res3=get_DE_genes_by_group(cell_exps_p,cell_exps[cell_path==sib_idx[-1]])
            #for pv_cut in [0.1,0.05]:
            #    for fc_cut in [0.6,1.0,1.5]:
            for pv_cut in ttest_pvs:
            #for pv_cut in [0.05]:
                for fc_cut in FC_cuts:
                #for fc_cut in [0.6]:
            #for pv_cut in [0.05]:
            #    for fc_cut in [1.5]:

                    #fc_cut_idx = g_diff > fc_cut
                    #expressed_fc_idx=np.logical_and(fc_cut_idx,g_diff>exp_cutoff)
                     
                    new_g_diff=np.fabs(mean_path_g[sib_idx[-1]]-mean_path_g[i]) 
                    fc_cut_idx = new_g_diff > fc_cut
                    'path_parent res3'
                    print "path_parent" , "pv cutoff: ",pv_cut, "fc cutoff", fc_cut
                    Jun_DE_idx=np.logical_and(fc_cut_idx,res3<pv_cut)
                    print "DE_genes_size: ",gene_names[Jun_DE_idx].shape
                    #pdb.set_trace()
                    #tfs,_=getEnrichTF(dTD,gene_names[Jun_DE_idx],gene_names,pcut=0.1,marker=marker_dict)
                    tfs,_=getEnrichTF(gene_names[Jun_DE_idx],gene_names,pcut=0.1,marker=marker_dict)
                    #pdb.set_trace()
                    print 'tf size: ',len(tfs)
                    
                    ##filter the tf not expressed: expression > 0 for 85% of cells

                    TF_pv_list+=_
                    tf_all+=tfs
                    method_all+=["path_parent_pv_cutoff_"+str(pv_cut)+"_FC_cutoff_"+str(fc_cut)]*len(tfs)
                    fc_no_map["path_parent_pv_cutoff_"+str(pv_cut)+"_FC_cutoff_"+str(fc_cut)]=cnt
                    cnt+=1
                    ## output TF list for 0.6, 1.0, 1.5 cutoff


                    if save is not None and pv_cut ==0.1 and fc_cut ==0.6:
                        out_DE = open(save+"path_"+str(i)+"_pv"+str(pv_cut)+"fc_"+str(fc_cut)+".txt",'w')
                        FC_list= new_g_diff[Jun_DE_idx]
                        fc_srt_idx =np.argsort(-FC_list)
                        DE_list=gene_names[Jun_DE_idx]
                        #out_DE.write("path: "+str(i)+"\n")
                        for idx in fc_srt_idx:
                            #print  DE_list[idx]+"\t "+str(FC_list[idx])+"\n"
                            out_DE.write( DE_list[idx]+"\t "+str(FC_list[idx])+"\n")
                        out_DE.close()
                        #out_DE.write("\n")
                '''
                'path_ttest res1'
                print "path_ttest" , "pv cutoff: ",pv_cut, "fc cutoff", fc_cut
                Jun_DE_idx=np.logical_and(fc_cut_idx,res1<pv_cut)
                print "DE_genes_size: ",gene_names[Jun_DE_idx].shape
                tfs=getEnrichTF(dTD,gene_names[DE_genes[Jun_DE_idx]],gene_names,pcut=0.1,marker=marker_dict)
                print 'tf size: ',len(tfs)
                tf_all+=tfs
                method_all+=["path_test pv cutoff: "+str(pv_cut)+" fc cutoff "+str(fc_cut)]*len(tfs)
                'path_all_other res2'
                print "path_all_other" , "pv cutoff: ",pv_cut, "fc cutoff", fc_cut
                Jun_DE_idx=np.logical_and(fc_cut_idx,res2<pv_cut)
                print "DE_genes_size: ",gene_names[Jun_DE_idx].shape
                tfs=getEnrichTF(dTD,gene_names[DE_genes[Jun_DE_idx]],gene_names,pcut=0.1,marker=marker_dict)
                print 'tf size: ',len(tfs)
                tf_all+=tfs
                method_all+=["path_all_other pv cutoff: "+str(pv_cut)+" fc cutoff "+str(fc_cut)]*len(tfs)
                '''


        '''
        print 'path_ttest_idx size: ',gene_names[path_ttest_idx].shape
        print 'g_fold_idx size: ',gene_names[g_fold_idx].shape
        print 'g_fold_top_idx size: ',gene_names[g_fold_top_idx].shape
        print 'g_diff_idx size: ',gene_names[g_diff_idx].shape
        print 'g_top_idx size: ',gene_names[g_top_idx].shape
        print 'k_top_idx size: ',gene_names[k_top_idx].shape
        print 'sp_cor_top_idx size: ',gene_names[sp_cor_top_idx].shape
        print 'g_diff_ratio_idx size: ',gene_names[g_diff_ratio_idx].shape
        
        print 'path_ttest tf size: ',len(path_ttest_tf)
        print 'g_fold tf size: ',len(g_fold_tf)
        print 'g_fold_top tf size: ',len(g_fold_top_tf)
        print 'g_diff tf size: ',len(g_diff_tf)
        print 'g_top tf size: ',len(g_top_tf)
        print 'k_top tf size: ',len(k_top_tf)
        print 'sp_cor_top tf size: ',len(sp_cor_top_tf)
        print 'g_diff_ratio tf size: ',len(g_diff_ratio_tf)

        #g_idx=res<0.1
        tf_all = g_diff_tf+g_fold_tf+g_fold_top_tf+g_top_tf+k_top_tf+sp_cor_top_tf+g_diff_ratio_tf+path_ttest_tf
        method_all = ["g_diff"]*len(g_diff_tf) \
                +["g_fold"]*len(g_fold_tf) \
                +["g_fold_top"]*len(g_fold_top_tf) \
                +["g_top"]*len(g_top_tf) \
                +["k_top"]*len(k_top_tf) \
                +["sp_cor_top"]*len(sp_cor_top_tf) \
                +["g_diff_ratio"]*len(g_diff_ratio_tf) \
                +["path_ttest_tf"]*len(path_ttest_tf) 
        '''
        print "tf_all: ",tf_all
        #tf_all=getEnrichTF(dTD,gene_names[g_idx],gene_names)

        TF_pv_dict=defaultdict(lambda:[])
        for pval,tf in TF_pv_list:
            TF_pv_dict[tf.lower()]+=[pval]
        tmp=[]
        #if adjust_p and len(TF_pv_list)>0:
        #    for key in TF_pv_dict.keys():
        #        tmp+=[min(TF_pv_dict[key])]
        #    #print 'before adjusted: ',tmp
        #    p_adjusted = multipletests(tmp, method='fdr_bh')
        #    #print "after adjusted: ",p_adjusted[1]
        #    for key,val in zip(TF_pv_dict.keys(),p_adjusted[1]):
        #        TF_pv_dict[key]=[val]
        for (pv,tf),method in zip(tf_all,method_all):
            #print tf, pv ,method
            #if adjust_p:
            #    TF_pvalue[tf.lower()]+=[TF_pv_dict[tf.lower()]]
            #else:
            TF_pvalue[tf.lower()]+=[pv]
            TF_method[tf.lower()]+=[method]

            
        if save is not None and sib_idx is not None:

            tf_names=[]
            tf_pvs=[]
            for key,val in TF_pv_dict.items():
                tf_names+=[key]
                tf_pvs+=[min(val)]
            tf_names=np.array(tf_names)
            tf_pvs=np.array(tf_pvs)
            srt_idx=np.argsort(tf_pvs)
            out_TF = open(save+"path_"+str(i)+"_TFlist.txt",'w')
            out_TF1 = open(save+"path_"+str(i)+"_TFlist_fc_0.6.txt",'w')
            out_TF2 = open(save+"path_"+str(i)+"_TFlist_fc_1.0.txt",'w')
            out_TF3 = open(save+"path_"+str(i)+"_TFlist_fc_1.5.txt",'w')
            FC_TFs=[out_TF1,out_TF2,out_TF3]
            for idx in srt_idx:
                tf_name = tf_names[idx].lower()
                #out_TF.write( tf_name+"\t "+str(tf_pvs[idx])+'\t'+str(len(TF_method[tf_name])) +'\t'+",".join(TF_method[tf_name])+"\n")
                out_TF.write( tf_name+"\t "+str(tf_pvs[idx])+"\n")
                #print TF_method[tf_name]
                #print fc_no_map
                #print fc_no_map[TF_method[tf_name]]
                for fn in TF_method[tf_name]:
                    FC_TFs[fc_no_map[fn]].write( tf_name+"\t "+str(tf_pvs[idx])+"\n")

            out_TF.close()
            out_TF1.close()
            out_TF2.close()
            out_TF3.close()

        TF_filtered=[]
        TF_filtered_pv=[]
        TF_filtered_method=[]
        TF_no_show=[]
        TF_no_show_pv=[]
        TF_no_show_method=[]
        TF_no_gene=[]
        TF_no_gene_pv=[]
        TF_no_gene_method=[]

        tmp=(g_pa + g_pb)/2
        for mk,mk_idx in zip(marker, marker_idx):
            print "marker exp: ",mk, ":",tmp[mk_idx]
    
    
        
        for tf,pv in TF_pvalue.items():
            min_idx=np.argmin(pv)
            #print 'pv:', TF_pv_dict[tf]
            if TF_pv_dict[tf][0]>pcut:
                continue
        
            if tf in gene_names:
                #print "tf_match!"
                tf_idx = gene_names.tolist().index(tf)
                #print "tf_idx:", tf_idx
                #print "gpa_tf_idx:", g_pa[tf_idx]
                #print "gpb_tf_idx:", g_pb[tf_idx]
                #print cell_exps_parent[:,tf_idx]>0
                #print len(cell_exps_parent[cell_exps_parent[:,tf_idx]>0])
                if g_pa[tf_idx] + g_pb[tf_idx] < exp_cutoff:
                    print tf, "not expressed with value: ",g_pa[tf_idx] + g_pb[tf_idx] 
                elif (len(cell_exps_p[cell_exps_p[:,tf_idx]>0])<cell_exps_p.shape[0]*0.1) and (sib_idx is None or sib_idx is not None and len(cell_exps_parent[cell_exps_parent[:,tf_idx]>0])<cell_exps_parent.shape[0]*0.1):
                    if sib_idx is not None:
                        print "len(cell_exps_parent[cell_exps_parent[:,tf_idx]>0])/#cell in path: : ", len(cell_exps_parent[cell_exps_parent[:,tf_idx]>0]),"/",cell_exps_parent.shape[0]
                    print "len(cell_exps_p[cell_exps_p[:,tf_idx]>0])/#cell in path: : ", len(cell_exps_p[cell_exps_p[:,tf_idx]>0]),"/",cell_exps_p.shape[0]
                    print tf, "not expressed with less than 10% of cells expressed in the path"
                else:
                #if cell_exps_p_g_sum[tf_idx] > exp_cutoff:
                    #print "pv length", len(pv)
                    #print "TF_method[tf] length", len(TF_method[tf])
                    #min_idx=np.argmin(pv)
                    TF_filtered+=[tf]
                    #if adjust_p:
                    #    #print 'add pv:', TF_pv_dict[tf]
                    #    #print 'ori pv: ', min[pv]
                    #    TF_filtered_pv+=[TF_pv_dict[tf][0]]
                    #else:
                    TF_filtered_pv+=[min(pv)]
                    #print tf
                    TF_filtered_method+=[TF_method[tf][min_idx]]
            else:
                print tf, " is not in gene names"
                TF_no_gene+=[tf]
                #if adjust_p:
                #    TF_no_gene_pv+=[TF_pv_dict[tf][0]]
                #else:
                TF_no_gene_pv += [min(pv)]
                TF_no_gene_method += [TF_method[tf][min_idx]]
        TF_filtered=np.array(TF_filtered)
        TF_filtered_pv=np.array(TF_filtered_pv)
        TF_filtered_method=np.array(TF_filtered_method)
        #print len(TF_filtered)
        #print len(TF_filtered_pv)
        #tmp = TF_filtered[TF_filtered_pv<pcut]
        if len(TF_filtered) > tf_max:
            srt_idx=np.argsort(TF_filtered_pv)
            TF_no_show=TF_filtered[srt_idx[tf_max:]]
            TF_no_show_pv=TF_filtered_pv[srt_idx[tf_max:]]
            TF_no_show_method=TF_filtered_method[srt_idx[tf_max:]]

            TF_filtered=TF_filtered[srt_idx[:tf_max]]
            TF_filtered_pv=TF_filtered_pv[srt_idx[:tf_max]]
            
            TF_filtered_method=TF_filtered_method[srt_idx[:tf_max]]
            #TF_no_show = TF_filtered
        print "tf no gene:", TF_no_gene
        print "tf no show:", TF_no_show
        #print "before TF_filtered length: ", len(TF_filtered)
        #print 'tmp length: ',len(tmp)
        #if len(tmp)<ntop:
        #    srt_idx=np.argsort(TF_filtered_pv)
        #    TF_filtered=TF_filtered[srt_idx[:ntop]]
        #    TF_filtered_pv=TF_filtered_pv[srt_idx[:ntop]]
        #elif len(tmp)>tf_max:
        #    srt_idx=np.argsort(TF_filtered_pv)
        #    TF_filtered=TF_filtered[srt_idx[:tf_max]]
        #    TF_filtered_pv=TF_filtered_pv[srt_idx[:tf_max]]

        #print "after TF_filtered length: ", len(TF_filtered)
        print "TF_filtered length", len(TF_filtered)
        print '\npath: ',i
        print 'expressed TF: ', len(TF_filtered)
        print zip(TF_filtered,TF_filtered_pv)
        target_list=[]
        #TF_pos_dict={}
        max_tps=[]
        #print dTD.keys()
        for tf_name in TF_filtered:
            val=dTD[tf_name.upper()]
        #for key, val in dTD.items():
            #tf_name=key.lower()
            #if tf_name in TF_filtered:
            target_list+=val
            #add continuous TF assignment here
#                 best_tp=0
#                 best_sum=0
            tf_idx=TF_names.tolist().index(tf_name)
            scores=np.zeros(6)
            for _ in range(0,6):
                tp=_/10.
#                     delta_t=1e-3
                reg_genes=set(map(lambda x:x.lower(),val))&set(gene_names)
                reg_genes=list(reg_genes)
                reg_genes_idxs=np.array(map(lambda x: gene_names.tolist().index(x),reg_genes))
                #print reg_genes_idxs

                #model['TF_start_time']=TF_start_time
                #gene_start_time=np.zeros((n_path,n_gene))
                gene_start_time[i,reg_genes_idxs]=tp
                #js=reg_genes_idxs
                #k=K_param_range/float(k_split)*ks
                for j in reg_genes_idxs:
                    #print g_pb[j]
                    #print g_pa[j]
                    #print k_p[j]
                    #print cell_time_p
                    #print gene_start_time[i,j]
                    mu_x_j=g_pb[j]+(g_pa[j]-g_pb[j])*np.exp(-k_p[j]*np.maximum(cell_time_p-gene_start_time[i,j],0))
                    x_j=cell_exps_p[:,j]
                    tmp=((x_j-mu_x_j)**2./(2.*sigma_param[j]**2.)+np.log((sigma_param[j]*np.sqrt(2.*np.pi)) ))
                    prob2 = np.where(x_j!=0.,0.,drop_out_param)
                    mix_prob=w_nz[i,j]*np.exp(-tmp)+(1-w_nz[i,j])*prob2
                    #print mix_prob.shape
                    sum_log_prob=np.sum(np.log(mix_prob))
                    scores[_]=+sum_log_prob
                #path_gene_k_table[p,j,ks-1]=sum_log_prob
            max_tp=np.argmax(scores)/10.
            max_tps+=[max_tp]
            #print "path: ",i," tf: ",tf_name, ' best tp: ',max_tp
            if "TF_start_time" in model.keys():
                model['TF_start_time'][i,tf_idx]=max_tp
        #sorted_tf_tp
        #print sorted_tf_tp
        #print map(str,sorted_tf_tp)
        #print "\t".join(map(str,sorted_tf_tp))
        print "TF_filtered"
        sorted_tf_tp=sorted(zip(TF_filtered,max_tps,TF_filtered_pv,TF_filtered_method),key= lambda tup: tup[1])
        print '\n'.join(map(str,sorted_tf_tp))
        print "TF_no_show"
        sorted_tf_no_show=sorted(zip(TF_no_show,TF_no_show_pv,TF_no_show_method),key= lambda tup: tup[1])
        print '\n'.join(map(str,sorted_tf_no_show))
        print "TF_no_gene"
        sorted_tf_no_gene=sorted(zip(TF_no_gene,TF_no_gene_pv,TF_no_gene_method),key= lambda tup: tup[1])
        print '\n'.join(map(str,sorted_tf_no_gene))
        print "#################################################"

                    #print len(val),len(reg_genes)
#                     delta_sum=0
#                     for gene in genes:
#                         j=gene_names.tolist().index(gene)
#                         mu_x_js=g_pb[j]+(g_pa[j]-g_pb[j])*np.exp(-k_p[j]*tp)
#                         mu_x_js_delta=g_pb[j]+(g_pa[j]-g_pb[j])*np.exp(-k_p[j]*(tp+delta_t))
#                         delta_y=mu_x_js_delta-mu_x_js
#                         delta_sum+=delta_y
#                     delta_sum=np.fabs(delta_sum)/len(genes)
#                     #print tp,delta_sum
                    
#                     if best_sum<delta_sum:
#                         best_sum=delta_sum
#                         best_tp=tp
#                 print 'best: ',best_tp,best_sum
                        

                    
        tg,tc =  np.unique(target_list,return_counts=True)
        target_list=map(lambda x: x.lower(),tg)
        #print 'target genes: ',tg.shape
        #print 'hahahawhahahahah', TF_filtered, target_list
        p["TF_list"]=TF_filtered
        p["TF_pvalue"]=TF_filtered_pv
        p["TF_method"]=TF_filtered_method
        p["TF_targets"]=target_list
        p["TF_no_show"]=TF_no_show
        p["TF_no_show_pv"]=TF_no_show_pv
        p["TF_no_show_method"]=TF_no_show_method
        p["TF_no_gene"]=TF_no_gene
        p["TF_no_gene_pv"]=TF_no_gene_pv
        p["TF_no_gene_method"]=TF_no_gene_method
        p["sorted_tf_tp"]=sorted_tf_tp
        p["sorted_eTFs"]=eTFs
        #p["TF_pos_dict"]=TF_pos_dict
def set_alpha_logistic_regression(model,gene_tf_table,cutoff=1):
    g_param=model["g_param"]
    path_info=model["path_info"]
    
    #place outside this function
#     gene_tf_table=[]
#     for i,gene in enumerate(gene_names):
#         feat=np.zeros(len(TF_names))
#         if gene in gene_tf_dict.keys():
#             tfs = map(lambda x: TF_names.tolist().index(x),gene_tf_dict[gene])
#             for t in tfs:
#                 feat[t]=1
#         gene_tf_table.append(feat)

    for i,p in enumerate(path_info):
        Sc_idx = p["Sc_idx"]
        Sp_idx = p["Sp_idx"]
        TF_list= p["TF_list"]
        tf_list_idx = map(lambda x: TF_names.tolist().index(x),TF_list)
        filtr=np.zeros(len(TF_names))
        for t in tf_list_idx:
            filtr[t]=1
        gene_feature = np.array(gene_tf_table)
        gene_feature = np.apply_along_axis(np.logical_and, 1,gene_feature,filtr)
        nz_idx=np.sum(gene_feature,axis=1)>0
        alpha = np.zeros(len(gene_names))
        print 'nz_idx shape: ',nz_idx.shape

        g_a=g_param[Sp_idx]
        g_b=g_param[Sc_idx]
        diff = np.fabs(g_a-g_b)>cutoff
        X=gene_feature[nz_idx]
        y=diff[nz_idx]
        if X.shape[0]==0 or np.unique(y).shape[0]==1:
            p["alpha"]=alpha
            continue
        print "X shape:", X.shape
        print "y shape:", y.shape
        logistic = LogisticRegression(class_weight ='balanced', penalty = "l2")
        logistic.fit(X,y)
        probs=logistic.predict_proba(X)
        alpha[nz_idx]=probs[:,1]
        p["alpha"]=alpha
    
    

#=======================================================================
# TF assignment module starts (by jund)

#-------------------------------------------------------------------
# get enriched TFs 
#---------------------------------------------------------------
# dMi: input sequence scanning result
# dMb: background sequence scanning result
# n: number of sequences in input 
# dTD: dictionary TF->DNA
# dMb: TF binding for background
# review erniched TF
# GL: all genes in the dataset

def getEnrichTF(GeneList,GL,pcut=0.1,ntop=3,marker={}):
    #print "assigning TF by gene list......"
    dMb=batchScanPrior([item.upper() for item in GL],dTD)
    dMi=batchScanPrior([item.upper() for item in GeneList],dTD)
    K=[item for item in dMi.keys() if item in dMb.keys()]
    K.sort()
    #print "testing ",len(K)," TFs......"
    n=len(GeneList)                                   # number of diff genes
    N=len(GL)                                         # N: number of sequences in background (all)
    entf=[]
    #print "K length: ",len(K)
    #pcut=0.01
    #print len(K)
    #pvalues=[]
    #all_tf=[]
    #marker_tf=[]
    all_tfs=[]
    for i in K:
        Ti=len(dMi[i])
        Tb=len(dMb[i])
        pr=float(Tb)/N
        pvi=1-pbinom(Ti-1,n,pr)
        #pvalues+=[pvi]
        if pvi<pcut:
            entf.append([pvi,i])
        #print "i: ",i
        all_tfs+=[[pvi,i]]
        if i.lower() in map(lambda x: x.lower(),marker.keys()):
            #print "marker p-value: ",i, pvi
            if pvi < marker[i.lower()]:
                marker[i.lower()]=pvi
            #marker_tf+=[[pvi,i]]
        #all_tf.append([pvi,i])   
    #p_adjusted = multipletests(pvalues, method='fdr_bh')
    #print sorted(pvalues)
    #print p_adjusted
    #for tf,p1,p2 in zip(K,pvalues,p_adjusted):
    #   print tf,p1,p2
    #print K[p_adjusted<pcut]
    #entf = zip(p_adjusted[p_adjusted<pcut],K[p_adjusted<pcut])
    #entf=[]
    #if len(entf)<ntop:
    #    all_tf.sort()
    #    entf=all_tf[:3]
    entf.sort()
    #print "returing ", len(entf)," TFs...."
    return entf,all_tfs
    
    
def batchScanPrior(A,dTD):
    # dTD  -> dictionary of TF-DNA interaction
    # A -> Gene list
    K=dTD.keys()
    K.sort()
    dA={item.upper():0 for item in A}
    dM={}
    for i in K:
        GI=[item.upper() for item in dTD[i]]
        GI=list(set([item for item in GI if item in dA]))
        if len(GI)>0:
            dM[i]=GI
    return dM

# get TF-DNA dictionary 
# TF->DNA
def getdTD(tfDNA):
    dTD={}
    with open(tfDNA,'r') as f:
        tfRows=f.readlines()
        tfRows=[item.strip().split() for item in tfRows]
        for row in tfRows:
            itf=row[0]
            itarget=row[1]
            if itf not in dTD:
                dTD[itf]=[itarget]
            else:
                dTD[itf].append(itarget)
    
    return dTD
    
#---------------binomial test p-value-----------------------------------
def normal_estimate(s, p, n):
    u = n * p
    o = (u * (1-p)) ** 0.5
    #pdb.set_trace()
    return 0.5 * (1 + erf((s-u)/(o*2**0.5)))
    

def pbinom(x,n,p):
        
    # handling special cases
    if x<0:
        return 0
    if n<=0:
        return 0
    if x>n:
        return 1
                                                
    # use scipy.binom_test to calculate binomial test p-value
    pv=binom_test(x,n,p,alternative="less")
    if (1-pv)<=sys.float_info.epsilon/2:
        return 1
    else:
        return pv

def pbinom_old(x,n,p):
    # this is approximation
    # if n is larger (<2000), approximation 1
    if n<2000:
        q=1.0-p
        pdf=cdf=q**n
        f=p/q
        for i in range(1,x+1):
            pdf*=((n-i+1.0)/i*f)
            cdf+=pdf
        return cdf
    else:
    # if n>=2000 (relatively large, approximiation 2
        return normal_estimate(x,p,n)
        
# TF assignment module end! 
#-------------------------------------------------------------

#=======================================================================

def gen_gene_tf_table():
    
    gene_tf_dict=defaultdict(lambda:[])
    for key,val in dTD.items():
        for v in val:
            gene_tf_dict[v.lower()]+=[key.lower()]
    gene_tf_table=[]
    for i,gene in enumerate(gene_names):
        feat=np.zeros(len(TF_names))
        if gene in gene_tf_dict.keys():
            tfs = map(lambda x: TF_names.tolist().index(x),gene_tf_dict[gene])
            for t in tfs:
                feat[t]=1
        gene_tf_table.append(feat)
    return gene_tf_table 
    
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d',"--data_file", help="specify the data file, if not specified then a default training data file will be used")
    parser.add_argument('-dt',"--data_file_testing", help="specify the testing data file and output best interation for testing, if not specifed then the model will not do testing.", default = None)
    parser.add_argument('-tf',"--TF_file", help="specify the tf-dna file, if not specifed then the model will not take TF into consideration.", default = None)
    
    parser.add_argument('-st',"--structure_file", help="specify the structure file, if not specified then a default structure file will be used")
    parser.add_argument('-seed',"--random_seed", help="specify the random seed, default is 0", type=int,default=0) 
    parser.add_argument('-ni',"--n_iteration", help="specify the number of training iteration, default is 10", type=int,default=10)
    parser.add_argument('-k',"--k_param_range", help="specify the range of K parameter, default is 10", type=int,default=10)
    parser.add_argument('-ns',"--n_split", help="specify the number of splits in learning K and assign cell time, default is 100", type=int,default=100)
    parser.add_argument('-na',"--n_anchor", help="specify the number of anchor cells to remain in each path during training, default is 0", type=int,default=0)
    parser.add_argument('-ng',"--n_gene", help="specify the maximum number of genes used in training, default is 1000", type=int,default=1000)
    parser.add_argument('-lamb',"--lamb", help="specify the regularizing parameter for L1 sparsity, default is 1", type=float,default=1)
    #parser.add_argument('-dop',"--drop_out_param", help="specify the drop-out parameter, default is 0.1", type=float,default=0.1, help=argparse.SUPPRESS)
    parser.add_argument('-dop',"--drop_out_param", type=float,default=0.1, help=argparse.SUPPRESS)
    parser.add_argument('-ps',"--assign_by_prob_sampling", help="specify the whether to use multinomial sampling in path assignment, default is 1", type=int,choices=[0,1],default=1)
    #parser.add_argument('-opt_w',"--optimize_w", help="specify the whether to optimize the w parameter in drop-out event, default is 0", type=int,choices=[0,1],default=0, help=argparse.SUPPRESS)
    parser.add_argument('-opt_w',"--optimize_w", type=int,choices=[0,1],default=0, help=argparse.SUPPRESS)
    #parser.add_argument('-ci',"--cluster_init", help="specify the whether to use k-means clustering as initialization of path assignment, default is 0", type=int,choices=[0,1],default=0, help=argparse.SUPPRESS)
    parser.add_argument('-ci',"--cluster_init", type=int,choices=[0,1],default=0, help=argparse.SUPPRESS)
    #parser.add_argument('-ldm',"--lamb_data_mult", help="specify the multiplier of lambda data parameter, default is logN",choices=['1','sqrtN','logN','N'],default='log(N)')
    #parser.add_argument('-lrm',"--lamb_ratio_mult", help="specify the multiplier of lambda ratio parameter, default is sqrtR",choices=['1','sqrtR','R'],default='sqrt(r)')
#     parser.add_argument('-ldm',"--lamb_data_mult", help="specify the multiplier of lambda data parameter, default is 1",choices=['1','sqrtN','logN','N'],default='1', help=argparse.SUPPRESS)
#     parser.add_argument('-lrm',"--lamb_ratio_mult", help="specify the multiplier of lambda ratio parameter, default is 1",choices=['1','sqrtR','R'],default='1', help=argparse.SUPPRESS)
#     parser.add_argument('-pc',"--path_constraint", help="specify the whether to apply path constraint in training, default is 0", type=int,choices=[0,1],default=0, help=argparse.SUPPRESS)
#     parser.add_argument('-pg',"--progress_bar", help="specify the whether to show progress_bar in training, default is 1", type=int,choices=[0,1],default=0, help=argparse.SUPPRESS)
    parser.add_argument('-ldm',"--lamb_data_mult",choices=['1','sqrtN','logN','N'],default='1', help=argparse.SUPPRESS)
    parser.add_argument('-lrm',"--lamb_ratio_mult",choices=['1','sqrtR','R'],default='1', help=argparse.SUPPRESS)
    parser.add_argument('-pc',"--path_constraint", type=int,choices=[0,1],default=0, help=argparse.SUPPRESS)
    parser.add_argument('-pg',"--progress_bar", type=int,choices=[0,1],default=0, help=argparse.SUPPRESS)
    parser.add_argument('-mn',"--model_name", help="specify the model_name",default = None)
    parser.add_argument('-cv',"--cross_validation", help="specify whether to use 5-fold cross_validation, 0 means not, default is 0", type=int, choices=[0,1],default=0)
    parser.add_argument('-opt',"--opt_method", help="specify what optimization method to solve lasso problem, genlasso or cvxpy, default is cvxpy", type=str, choices=["genlasso","cvxpy"],default="cvxpy")
        
    
    args=parser.parse_args()
    print args
    data_file='data/treutlein2014'
    if args.data_file is not None:
        data_file=args.data_file
    splits=data_file.split('/')
    #structure_file=splits[0]+'/init_cluster_'+splits[1]+'.txt'
    if args.structure_file is not None:
        structure_file=args.structure_file
    else:
        structure_file=splits[0]+'/trained_cluster_'+splits[1]+'.txt'
    verbose=1
    np.random.seed(args.random_seed)
    n_iteration=args.n_iteration
    n_split=args.n_split
    K_param_range=args.k_param_range
    n_anchor=args.n_anchor
    lamb=args.lamb
    max_gene=args.n_gene
    drop_out_param=args.drop_out_param
    assign_by_prob_sampling=args.assign_by_prob_sampling 
    optimize_w=args.optimize_w
    path_constraint=args.path_constraint
    cluster_init=args.cluster_init
    lamb_data_mult=args.lamb_data_mult
    lamb_ratio_mult=args.lamb_ratio_mult
    progress_bar=args.progress_bar
    cv = args.cross_validation
    
    if args.TF_file is not None:
        dTD=getdTD(args.TF_file)
        TF_names = np.array(map(lambda x: x.lower(),dTD.keys()))
#     gene_tf_table = gen_gene_tf_table()
#     gene_tf_table=[]
#     for i,gene in enumerate(gene_names):
#         feat=np.zeros(len(TF_names))
#         if gene in gene_tf_dict.keys():
#             tfs = map(lambda x: TF_names.tolist().index(x),gene_tf_dict[gene])
#             for t in tfs:
#                 feat[t]=1
#         gene_tf_table.append(feat)
 

    if args.model_name is not None:
        model_name=args.model_name
    else:
        model_name = 'model/model_'+splits[1]+'_ns_'+str(n_split)+'_lamb_'+str(lamb)+'_ng_'+str(max_gene)+'_cv_'+str(cv)


    verbose=1
    progress_bar=1
    if args.data_file_testing is not None:
        verbose = 0
        cell_names_train,cell_day_train,cell_labels_train,cell_exps_train,gene_names=load_data(data_file,max_gene)
        cell_names_test,cell_day_test,cell_labels_test,cell_exps_test,gene_names=load_data(args.data_file_testing,max_gene)
        model,hid_var_train = init_var_Jun(structure_file,cell_names_train,cell_day_train,cell_exps_train,cell_labels_train)
        _,hid_var_test = init_var_Jun(structure_file,cell_names_test,cell_day_test,cell_exps_test,cell_labels_test)
        max_it=args.n_iteration
        #max_it=2
        max_test_ll=-float('inf')
        for it in range(1,max_it+1):
            n_iteration = 1
            assign_by_prob_sampling=args.assign_by_prob_sampling 
            optimize_likelihood(cell_exps_train, gene_names, model, hid_var_train,model_name,store_model=False)
            assign_by_prob_sampling=False
            assign_path_and_time(model,hid_var_test,cell_exps_test)
            train_ll = log_likelihood(model,hid_var_train,cell_exps_train)
            test_ll = log_likelihood(model,hid_var_test,cell_exps_test)
            print 'iteration:\t ', it, '\t train_LL:\t ', np.around(train_ll,2), '\t test_ll: \t', np.around(test_ll,2)
            if test_ll > max_test_ll:
                max_test_ll = test_ll
                count = 0
                best_it = it
            else:
                count+=1
            if count>1:
                break
        print 'best_test_it: ', best_it, '\t max_test_ll: ', max_test_ll
        print best_it
        #best_its.append(best_it)
        #best_test_lls.append(max_test_ll)
        #print 'best_its: ', best_its
        #print 'best_test_lls: ', best_test_lls
        #print 'mean_best_its: ',np.mean(best_its)
        #print 'mean_best_test_lls: ',np.mean(best_test_lls)
        #print 'training all data with the best it:', int(np.rint(np.mean(best_its)))
        sys.exit(0)
    
    if args.cross_validation:
        cell_names,cell_day,cell_labels,cell_exps,gene_names=load_data(data_file,max_gene)
        if args.TF_file is not None:
            gene_tf_table = gen_gene_tf_table()
        n_fold = 5
        cv_idx = cv_split_idx(cell_day=cell_day,n_fold=n_fold)
        best_its=[]
        best_test_lls=[]
        verbose = 0
        for i in range(1,n_fold+1):
            print 'fold: ', i
            test_idx = cv_idx==i
            train_idx = cv_idx!=i
            cell_day_test = np.array(cell_day)[test_idx]
            cell_exps_test = np.array(cell_exps)[test_idx]
            cell_labels_test = np.array(cell_labels)[test_idx]
            cell_names_test = np.array(cell_names)[test_idx]
            cell_day_train = np.array(cell_day)[train_idx]
            cell_exps_train = np.array(cell_exps)[train_idx]
            cell_labels_train = np.array(cell_labels)[train_idx]
            cell_names_train = np.array(cell_names)[train_idx]
            
            model,hid_var_train = init_var_Jun(structure_file,cell_names_train.tolist(),cell_day_train,cell_exps_train,cell_labels_train)
            _,hid_var_test = init_var_Jun(structure_file,cell_names_test.tolist(),cell_day_test,cell_exps_test,cell_labels_test)
            
            #model,hid_var_train = init_var(adj_mat,cell_day_train,cell_exps_train,cell_labels_train)
            #_,hid_var_test = init_var(adj_mat,cell_day_test,cell_exps_test,cell_labels_test,testing=True)
            
            #train_ll = log_likelihood(model,hid_var_train,cell_exps_train)
            #test_ll = log_likelihood(model,hid_var_test,cell_exps_test)
            #print 'initial training full log-likelihood: ',train_ll
            #print 'initial testing full log-likelihood: ',test_ll
            #print 'iteration:\t ', 0, '\t train_LL:\t ', np.around(train_ll,2), '\t test_ll: \t', np.around(test_ll,2)
            #compute_ARI_confuss_mat(hid_var_train)
            n_iteration = 1
            max_it=args.n_iteration
            max_test_ll=-float('inf')
            count=0
            best_it = 0
            for it in range(1,max_it):
                assign_by_prob_sampling=args.assign_by_prob_sampling 
                optimize_likelihood(cell_exps_train, gene_names, model, hid_var_train,model_name,store_model=False)
                assign_by_prob_sampling=False
                assign_path_and_time(model,hid_var_test,cell_exps_test)
                train_ll = log_likelihood(model,hid_var_train,cell_exps_train)
                test_ll = log_likelihood(model,hid_var_test,cell_exps_test)
                print 'iteration:\t ', it, '\t train_LL:\t ', np.around(train_ll,2), '\t test_ll: \t', np.around(test_ll,2)
                if test_ll > max_test_ll:
                    max_test_ll = test_ll
                    count = 0
                    best_it = it
                else:
                    count+=1
                if count>1:
                    break
            print 'best_test_it: ', best_it, '\t max_test_ll: ', max_test_ll
            best_its.append(best_it)
            best_test_lls.append(max_test_ll)
        print 'best_its: ', best_its
        print 'best_test_lls: ', best_test_lls
        print 'mean_best_its: ',np.mean(best_its)
        print 'mean_best_test_lls: ',np.mean(best_test_lls)
        print 'training all data with the best it:', int(np.rint(np.mean(best_its)))
        verbose = 1
        n_iteration= int(np.rint(np.mean(best_its)))
        #model,hid_var = init_var(adj_mat,cell_day,cell_exps,cell_labels)
        model,hid_var = init_var_Jun(structure_file,cell_names,cell_day,cell_exps,cell_labels)
        #print 'initial full log-likelihood: ',log_likelihood(model,hid_var,cell_exps)
        n_path=len(model['path_info'] )
        compute_ARI_confuss_mat(hid_var,n_path)
        optimize_likelihood(cell_exps, gene_names, model, hid_var,model_name)
    else:
        cell_names,cell_day,cell_labels,cell_exps,gene_names=load_data(data_file,max_gene)
        if args.TF_file is not None:
            gene_tf_table = gen_gene_tf_table()
        n_cell=len(cell_names)
        model,hid_var = init_var_Jun(structure_file,cell_names,cell_day,cell_exps,cell_labels)
        n_path=len(model['path_info'] )
        compute_ARI_confuss_mat(hid_var,n_path)
        optimize_likelihood(cell_exps, gene_names, model, hid_var,model_name)
