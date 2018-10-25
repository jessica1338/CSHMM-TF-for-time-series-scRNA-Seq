import CSHMM_TF_train_release as ML
import numpy as np
import matplotlib
matplotlib.use('Agg')
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt
from collections import defaultdict
import os
import pandas as pd
import scipy
import pygraphviz as PG
import matplotlib.image as mpimg
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import CubicSpline,splrep,BSpline
from math import acos
from math import sqrt
from math import pi
from matplotlib.ticker import MaxNLocator


import logging
reload(logging)
debug=False
if debug:
    logging.basicConfig(level=logging.DEBUG)
print_top_limit=20

def init_marker_gene_dict():

    marker_genes_dict=defaultdict(lambda: ("",""))

    #treutlein 2014
    marker_genes_dict['Abca3'] = ('AT2','treutlein2014: alveolar/bronchiolar lineages')
    marker_genes_dict['Sftpb'] = ('AT2','treutlein2014: alveolar/bronchiolar lineages')
    marker_genes_dict['Muc1'] = ('AT2','treutlein2014: alveolar/bronchiolar lineages')
    marker_genes_dict['Kyz2'] = ('AT2','treutlein2014: alveolar/bronchiolar lineages')
    marker_genes_dict['Sftpc'] = ('AT2','treutlein2014: alveolar/bronchiolar lineages')
    marker_genes_dict['Aqp5'] = ('AT1','treutlein2014: alveolar/bronchiolar lineages')
    marker_genes_dict['Pdpn'] = ('AT1','treutlein2014: alveolar/bronchiolar lineages')
    marker_genes_dict['Ager'] = ('AT1','treutlein2014: alveolar/bronchiolar lineages')
    marker_genes_dict['Foxj1'] = ('ciliated','treutlein2014: alveolar/bronchiolar lineages')
    marker_genes_dict['Scgb1a1'] = ('Clara','treutlein2014: alveolar/bronchiolar lineages')
    marker_genes_dict['Cftr'] = ('AT2','treutlein2014: novel')
    marker_genes_dict['Cebpa'] = ('AT2','treutlein2014: novel')
    marker_genes_dict['Sftpd'] = ('AT2','treutlein2014: novel')
    marker_genes_dict['Id2'] = ('AT2','treutlein2014: novel')
    marker_genes_dict['Begfa'] = ('AT1','treutlein2014: novel')
    marker_genes_dict['Itgb4'] = ('ciliated','treutlein2014: novel')
    marker_genes_dict['Top2a'] = ('ciliated','treutlein2014: novel')

    #sabrina TASIC
    marker_genes_dict['Gabpb1'] = ('DEgenes','sabrina TASIC: t-test/Scell')
    marker_genes_dict['Rars2'] = ('DEgenes','sabrina TASIC: t-test/Scell')
    marker_genes_dict['Naa50'] = ('DEgenes','sabrina TASIC: t-test/Scell')
    marker_genes_dict['Jund'] = ('DEgenes','sabrina TASIC: t-test/Scell')
    marker_genes_dict['Myo1b'] = ('DEgenes','sabrina TASIC: t-test/Scell')
    marker_genes_dict['Adamfs10'] = ('DEgenes','sabrina TASIC: t-test')
    marker_genes_dict['Foxa2'] = ('DEgenes','sabrina TASIC: t-test/Scell')
    marker_genes_dict['Ppp3r1'] = ('DEgenes','sabrina TASIC: t-test')
    marker_genes_dict['Cdk4'] = ('DEgenes','sabrina TASIC: t-test/Scell')
    marker_genes_dict['Nasp'] = ('DEgenes','sabrina TASIC: t-test/')
    marker_genes_dict['Dlk1'] = ('DEgenes','sabrina TASIC: Scell')
    marker_genes_dict['Hmgb2'] = ('DEgenes','sabrina TASIC: Scell')
    marker_genes_dict['Cdc6'] = ('DEgenes','sabrina TASIC: Scell')

    marker_genes_dict['Birc5'] = ('mitosis','treutlein2016: , main text')
    marker_genes_dict['Ube2c'] = ('mitosis','treutlein2016: , main text')
    # marker_genes_dict['Hmga2'] = ('mitosis','treutlein2016: , main text)
    marker_genes_dict['Cadm1'] = ('neural_projections','treutlein2016: , main text')
    # marker_genes_dict['Dner'] = ('neural projections','treutlein2016: , main text)
    marker_genes_dict['Klhl24'] = ('neural_projections','treutlein2016: , main text')
    # marker_genes_dict['Tubb3'] = ('neural projections','treutlein2016: , main text)
    marker_genes_dict['Mapt'] = ('neural_projections','treutlein2016: , main text')                             
    # marker_genes_dict['Snca'] = ('synaptic transmission','treutlein2016: , main text)
    marker_genes_dict['Stxbp1'] = ('synaptic_transmission','treutlein2016: , main text')
    marker_genes_dict['Vamp2'] = ('synaptic_transmission','treutlein2016: , main text')
    marker_genes_dict['Dmpk'] = ('synaptic_transmission','treutlein2016: , main text')
    marker_genes_dict['Ppp3ca'] = ('synaptic_transmission','treutlein2016: , main text')                              
    marker_genes_dict['Sept3'] = ('cytoskeletal_reorganization','treutlein2016: , main text')
    marker_genes_dict['Sept4'] = ('cytoskeletal_reorganization','treutlein2016: , main text')
    marker_genes_dict['Coro2b'] = ('cytoskeletal_reorganization','treutlein2016: , main text')
    marker_genes_dict['Ank2'] = ('cytoskeletal_reorganization','treutlein2016: , main text')
    marker_genes_dict['Mtap1a'] = ('cytoskeletal_reorganization','treutlein2016: , main text')
    marker_genes_dict['Homer2'] = ('cytoskeletal_reorganization','treutlein2016: , main text')
    marker_genes_dict['Akap9'] = ('cytoskeletal_reorganization','treutlein2016: , main text')                             
    marker_genes_dict['Ascl1'] = ('Ascl1_target|initial_factor','treutlein2016: Extend Data Figure 3b, 8c')
    marker_genes_dict['Hes6'] = ('Ascl1_target|initial_factor','treutlein2016: Extend Data Figure 3b, 8c')
    marker_genes_dict['Zfp238'] = ('Ascl1_target|initial_factor','treutlein2016: Extend Data Figure 3b, 8c')
    marker_genes_dict['Snca'] = ('Ascl1_target|synaptic_transmission','treutlein2016: Extend Data Figure 3b , main text')
    marker_genes_dict['Cox8b'] = ('Ascl1_target','treutlein2016: Extend Data Figure 3b')
    marker_genes_dict['Bex1'] = ('Ascl1_target','treutlein2016: Extend Data Figure 3b')
    marker_genes_dict['Dner'] = ('Ascl1_target|neural_projections','treutlein2016: Extend Data Figure 3b, main text')
    marker_genes_dict['Atoh8'] = ('initial_factor','treutlein2016: Extend Data Figure 8c')
    marker_genes_dict['Sox9'] = ('initial_factor|NPC','treutlein2016: Extend Data Figure 8c,main text')
    marker_genes_dict['Tcf4'] = ('initial_factor','treutlein2016: Extend Data Figure 8c')
    marker_genes_dict['Sox11'] = ('initial_factor','treutlein2016: Extend Data Figure 8c')
    marker_genes_dict['Tcf12'] = ('initial_factor','treutlein2016: Extend Data Figure 8c')
    marker_genes_dict['Dlx3'] = ('initial_factor','treutlein2016: Extend Data Figure 8c')                        
    marker_genes_dict['Ecm1'] = ('MEF','treutlein2016: Extend Data Figure 3b')
    marker_genes_dict['Scd1'] = ('MEF','treutlein2016: Extend Data Figure 3b')
    marker_genes_dict['Dcn'] = ('MEF|Fibroblast','treutlein2016: Extend Data Figure 3b,6i')                         
    marker_genes_dict['Hmga2'] = ('MEF_factors|mitosis','treutlein2016: Extend Data Figure 8e, main text')
    marker_genes_dict['Id3'] = ('MEF_factors','treutlein2016: Extend Data Figure 8e')                            
    marker_genes_dict['Myh3'] = ('Myocyte','treutlein2016: Extend Data Figure 6g')
    marker_genes_dict['Myo18b'] = ('Myocyte','treutlein2016: Extend Data Figure 6g')
    marker_genes_dict['Tnnc2'] = ('Myocyte','treutlein2016: Extend Data Figure 6g, main text')
    marker_genes_dict['Acta1'] = ('Myocyte','treutlein2016: Extend Data Figure 6g')                            
    marker_genes_dict['Map2'] = ('Neuron|pan-neuronal','treutlein2016: Extend Data Figure 6h,main text')                           
    marker_genes_dict['Syp'] = ('Neuron|synaptic_maturation','treutlein2016: main text')
    marker_genes_dict['Gria2'] = ('Neuron|synaptic|synaptic_maturation','treutlein2016: Extend Data Figure 6h, main text')
    marker_genes_dict['Snap25'] = ('Neuron|synaptic|synaptic_maturation','treutlein2016: Extend Data Figure 6h, main text')
    marker_genes_dict['Tubb3'] = ('Neuron|pan-neuronal|neural_projections','treutlein2016: Extend Data Figure 6h, main text')
    marker_genes_dict['Nrxn3'] = ('Neuron|synaptic|synaptic_maturation','treutlein2016: Extend Data Figure 6h, main text')
    marker_genes_dict['Stmn3'] = ('Neuron|synaptic','treutlein2016: Extend Data Figure 6h, main text')
    marker_genes_dict['Rab3c'] = ('synaptic_maturation','treutlein2016: main text')
    marker_genes_dict['Syt4'] = ('synaptic_maturation','treutlein2016: main text')
    marker_genes_dict['Sv2a'] = ('synaptic_maturation','treutlein2016: main text')
    # marker_genes_dict['Map2'] = ('pan-neuronal','treutlein2016: main text)
    # genes involved in synaptic maturation are turned on (Syp, Rab3c, Gria2, Syt4, Nrxn3, Snap25, Sv2a)                          
    marker_genes_dict['Camta1'] = ('Neuron','treutlein2016: Extend Data Figure 8d')
    marker_genes_dict['Insm1'] = ('Neuron','treutlein2016: Extend Data Figure 8d')
    marker_genes_dict['Myt1l'] = ('Neuron','treutlein2016: Extend Data Figure 8d')
    marker_genes_dict['St18'] = ('Neuron','treutlein2016: Extend Data Figure 8d')
    marker_genes_dict['Peg3'] = ('Neuron','treutlein2016: Extend Data Figure 8d')
    marker_genes_dict['Gli3'] = ('NPC','treutlein2016: main text')
    # marker_genes_dict['Sox9'] = ('NPC','treutlein2016: main text)
    marker_genes_dict['Nestin'] = ('NPC','treutlein2016: main text')
    marker_genes_dict['Fabp7'] = ('NPC','treutlein2016: main text')
    marker_genes_dict['Hes1'] = ('NPC','treutlein2016: main text')
    marker_genes_dict['Sox2'] = ('canonical_NPC','treutlein2016: main text')
    marker_genes_dict['Pax6'] = ('canonical_NPC','treutlein2016: main text')                            
    marker_genes_dict['Col1a2'] = ('Fibroblast','treutlein2016: Extend Data Figure 6i')
    # marker_genes_dict['Dcn'] = ('Fibroblast','treutlein2016: Extend Data Figure 6i)                           
    marker_genes_dict['Eln'] = ('Fibroblast','treutlein2016: main text')
    
    
    marker_gene_list = marker_genes_dict.keys()
    treutlein2014_mkgene = []
    treutlein2014_mkgene_dict=defaultdict(lambda: ("",""))
    treutlein2016_mkgene = []
    treutlein2016_mkgene_dict=defaultdict(lambda: ("",""))
    for key, val in marker_genes_dict.items():
        if val[1].split(':')[0]=='treutlein2014':
            treutlein2014_mkgene+=[key]
            treutlein2014_mkgene_dict[key] = val
        if val[1].split(':')[0]=='treutlein2016':
            treutlein2016_mkgene+=[key]
            treutlein2016_mkgene_dict[key] = val
    print treutlein2014_mkgene
    print treutlein2016_mkgene
    print len(treutlein2016_mkgene)
    return marker_gene_list,marker_genes_dict,treutlein2014_mkgene,treutlein2014_mkgene_dict,treutlein2016_mkgene,treutlein2016_mkgene_dict


def get_graphviz_plain_node_xy(filename):
    # model,hid_var=ML.load_model('model/'+model_file)
    # path_info=model['path_info']
    # print path_info
    #     pars=[0]
    #     cell_names,cell_day,cell_labels,cell_exps,gene_names=ML.load_data(data_file,ng)
    #     cell_path = hid_var['cell_path']
    #     cell_labels = hid_var['cell_labels']
    #     cell_time=hid_var['cell_time']
    #     cell_ori_time=hid_var['cell_ori_time']
    #     cell_ori_time=np.array(hid_var['cell_ori_time'])
    lines= open(filename).readlines()
    #print lines
    plt.clf()
    node_pos={}
    #edge={}
    for i,line in enumerate(lines):
        sp=line.split('\t')
        #print sp
        if len(sp)>1 and sp[1].startswith("D"):
            if len(sp[1].split(' '))>1:
    #             print lines[i]
    #             print lines[i+1]
    #             print lines[i+2]
    # #             pos = lines[i+1].split('\t')[2].split('"')[1].split(',')
    # #             edge_x+=[float(pos[0])]
    # #             edge_y+=[float(pos[1])]
    #             pos = lines[i+2].split('\t')[2].split('"')[1].replace(","," ").split(" ")
    #             print pos
    #             for j in range(1,len(pos),2):
    #                 edge_x+=[float(pos[j])]
    #                 edge_y+=[float(pos[j+1])]
                continue

            #print sp
            #print lines[i+1]
            pos = map(float,lines[i+1].split('\t')[2].split('"')[1].split(','))
            node_pos[sp[1]]=pos
            #x+=[float(pos[0])]
            #y+=[float(pos[1])]
    return node_pos       
    #plt.scatter(x,y)
    #plt.scatter(edge_x,edge_y)
#node_pos=get_graphviz_plain_node_xy('plain.txt')
#node_pos

def plot_path_fig(model_file,data_file,circle_size=40,label_add_time=False, change_label=False, time_as_label=False,fontsize = 'x-large', use_data_labels=False):
    print 'plotting path figure for model file: ',model_file
    y_vals=[]
    y_axis=[]
    markers=['.',',','o','s','p','*','+','D','x','h']
    splits=model_file.split('_')
    splits2=('model_'+'_'.join(splits[1:])).split('.')
    #print splits2
    out_fig_name=splits2[0]+'.'+splits2[1]
    #print out_fig_name
    data_file=data_file
    #ng=int(splits[-4])
    model,hid_var=ML.load_model(model_file)
    path_info=model['path_info']
    pars=[0]
    cell_names,cell_day,cell_labels,cell_exps,gene_names=ML.load_data(data_file,16000)
    #print type(cell_labels)
    cell_path = hid_var['cell_path']
    if not use_data_labels:
        cell_labels = hid_var['cell_labels']
    else:
        cell_labels=np.array(cell_labels)
    #print cell_labels
    
    cell_time=hid_var['cell_time']
    cell_ori_time=hid_var['cell_ori_time']
    cell_ori_time=np.array(hid_var['cell_ori_time'])
    #print hid_var.keys()
    n_cell, n_gene = cell_exps.shape
    #print cell_labels.shape
    if change_label:
        new_labels=[]
        for i in range(n_cell):
            if cell_labels[i]=='NA' or label_add_time:
                #print "change_label1"
                #print cell_labels[i] , cell_ori_time[i] , str(cell_ori_time[i])
                #cell_labels[i]=str(cell_ori_time[i])+'_'+cell_labels[i].strip()
                #print cell_labels[i].split(",")[0]+'_'+str(cell_ori_time[i])
                new_labels+=[cell_labels[i].split(",")[0]+'_'+str(cell_ori_time[i])]
                #cell_labels[i]+='_'+str(cell_ori_time[i])
                #print cell_labels[i]
                #print len(new_labels)
            else:
                #print "change_label2"
                #print cell_labels[i]
                new_labels+=[cell_labels[i]]
                
                #print len(new_labels)
        #print new_labels
        cell_labels=np.array(new_labels)
        #print cell_labels.shape
    if time_as_label:
        cell_labels=cell_ori_time
    paths=np.unique(cell_path)
    n_path=len(paths)
    if n_cell != cell_path.shape[0]:
        print '#cell mismatch, failed to produce figure'
        return
    edges=[]
    
    A = PG.AGraph(directed=True, strict=True)
    A.graph_attr['rankdir']='LR'
    while(len(pars)>0):
        par = pars[0]
        path=path_info[par]
        childs = path['child_path']
        a= str(path['Sp_idx'])
        b= str(path['Sc_idx'])
        A.add_edge("D"+a, "D"+b, label='path'+str(par))
        pars+=childs
        pars=pars[1:]
    plt.figure(figsize=(10,7))
    A.layout(prog='dot')
    A.write('structure.txt')
    n_pos=get_graphviz_plain_node_xy('structure.txt')
    #A.draw('figure/'+out_fig_name+'_structure.png')
    A.draw('structure.png')
    print 'structure.png'
    img=mpimg.imread('structure.png')
    imgplot = plt.imshow(img)
    #plt.show()
    
    plt.clf()
    xs=[]
    ys=[]
    plt.figure(figsize=(15,7),tight_layout=True)
    font = {'family': 'serif',
    'color':  'black',
    'weight': 'normal',
    'size': 16,
    }
    def length(v):
        return sqrt(v[0]**2+v[1]**2)
    def dot_product(v,w):
        return v[0]*w[0]+v[1]*w[1]
    def determinant(v,w):
        return v[0]*w[1]-v[1]*w[0]
    def inner_angle(v,w):
        cosx=dot_product(v,w)/(length(v)*length(w))
        rad=acos(cosx) # in radians
        return rad*180/pi # returns degrees
    def angle_clockwise(A, B):
        inner=inner_angle(A,B)
        if inner!=0:
            inner+=15
        det = determinant(A,B)
        if det<0: #this is a property of the det. If the det < 0 then B is clockwise of A
            return inner
        else: # if the det > 0 then A is immediately clockwise of B
            return 360-inner
    for p in range(len(path_info)):
        #p=path_info[i]
        #print i, 
        p_idx=(cell_path==p)
        if cell_path[p_idx].shape[0]==0:
            continue
        sp=path_info[p]['Sp_idx']
        sc=path_info[p]['Sc_idx']
        if 'D'+str(sp) in n_pos.keys():
            sp_pos = np.array(n_pos['D'+str(sp)])
        if 'D'+str(sc) in n_pos.keys():
            sc_pos = np.array(n_pos['D'+str(sc)])
        plt.plot([sp_pos[0],sc_pos[0]],[sp_pos[1],sc_pos[1]],'black')
        plt.text(sp_pos[0]-10, sp_pos[1]+5, 'D'+str(sp), fontdict=font)
        plt.text(sc_pos[0]-10, sc_pos[1]+5, 'D'+str(sc), fontdict=font)
        middle = (sp_pos+ sc_pos)/2
        angle = angle_clockwise(sc_pos-sp_pos,[1,0])
        #print angle
        #angle = 0
        x_change=0
        y_change=0
        #if angle > 180:
        #    y_change*=-1
        #if angle ==360:
        #    y_change =5
        #plt.text(middle[0]+x_change, middle[1]+y_change, 'path '+str(p), fontdict=font,rotation=angle)
        plt.text(middle[0]+x_change, middle[1]+y_change, 'P'+str(p), fontdict=font)
        

    for i in range(n_cell):
        p=cell_path[i]
        t=cell_time[i]
        sp=path_info[p]['Sp_idx']
        sc=path_info[p]['Sc_idx']
        sp_pos = np.array(n_pos['D'+str(sp)])
        sc_pos = np.array(n_pos['D'+str(sc)])
        #weighted_pos=np.round(sp_pos*(1-t)+sc_pos*t,2)
        #weighted_pos=map(int,sp_pos*(1-t)+sc_pos*t)
        weighted_pos=sp_pos*(1-t)+sc_pos*t
        #weighted_pos=map(int,sc_pos*(1-t)+sp_pos*t)
        ys.append(weighted_pos[1])
        xs.append(weighted_pos[0])
        #print weighted_pos
    xs=np.array(xs)
    ys=np.array(ys)
    labs=np.unique(cell_labels)
    lab_xy_pair={}
    #pcolormesh(x, y, Z, vmin=-1., vmax=1., cmap='RdBu_r')
    for i,lab in enumerate(labs):
        #cell_time_lab=cell_time[cell_labels==lab]
        #ys_lab=ys[cell_labels==lab]
        lab_xs=[]
        lab_ys=[]
        lab_ss=[]
        #lab_cs=[]
        lab_xy_pair[lab]=defaultdict(lambda:0)
        #print lab
        tmp=cell_labels==lab
        #print tmp.shape
        #print xs.shape
        #print ys.shape
        #print cell_labels.shape
        #print xs[cell_labels==lab].shape
        #print ys[cell_labels==lab].shape
        for a,b in zip(xs[cell_labels==lab],ys[cell_labels==lab]):
            lab_xy_pair[lab][(a,b)]+=1
        #print lab_xy_pair[lab]
        for key,val in lab_xy_pair[lab].items():
            lab_xs.append(key[0])
            lab_ys.append(key[1])
            #lab_ss.append(val/float(n_cell)*9000)
            lab_ss.append(np.sqrt(val)*circle_size)
            #lab_cs.append(i/float(len(labs)))
    #plt.scatter(xs,ys)
        #print len(lab_xs)
        #print len(lab_ys)
        #print len(lab_ss)
        #print lab_cs
        lab_xs+=[lab_xs[-1]]
        lab_ys+=[lab_ys[-1]]
        lab_ss+=[circle_size]
        plt.scatter(lab_xs,lab_ys,s=lab_ss,marker=markers[2+i/10],label=lab)
        #plt.scatter(lab_xs,lab_ys,s=lab_ss,c=lab_cs,marker=markers[2],label=lab,cmap=plt.get_cmap("tab20"))
        #plt.scatter(lab_xs,lab_ys,s=lab_ss,marker=markers[2],label=lab)
    lgd=plt.legend(loc='upper left',ncol=2,fancybox=True,shadow=True, fontsize = fontsize)
    for lgd_hand in lgd.legendHandles:
        lgd_hand._sizes = [circle_size]
    
    plt.axis('off')

    plt.savefig(model_file+'_cell_tree.png',dpi=300)

        
    plt.show()
    
def GO_ana_gene_list(gene_list,out_file):
    head=["query", "significant", "p_value", "T", "Q", "Q&T", "precision", "recall", "term_id","domain", "group", "description", "depth", "intersection", "evcodes"]
    cmd = 'printf "'+'\\t'.join(head)+'\\n" > '+out_file
    #print cmd
    logging.debug(cmd)
    os.system(cmd)
    cmd = 'python ~/repos/bio_packages/gprofiler-official-0.2.3/gprofiler.py -o mmusculus "'+ ' '.join(gene_list) + '" >> '+out_file
    #print cmd
    logging.debug(cmd)
    os.system(cmd)
    
def analyze_gene(model_file,data_file=None,ng=16000,model_ana_temp=None,out_folder=None):
    print 'analyzing for model file: ',model_file
    #y_vals=[]
    #y_axis=[]
    #markers=['.',',','o','s','p','*','+','D','x','h']
    splits=model_file.split('_')
    splits2=('model_'+'_'.join(splits[1:])).split('.')
    #print splits2
    out_fig_name=splits2[0]+'.'+splits2[1]
    #print out_fig_name
    if out_folder is None:
        out_folder = 'figure/'+out_fig_name
    if data_file is None:
        data_file='data/'+splits[1]
    #ng=int(splits[-4])
    model,hid_var=ML.load_model(model_file)
    path_info=model['path_info']
    pars=[0]
    #par=0
    cell_names,cell_day,cell_labels,cell_exps,gene_names=ML.load_data(data_file,ng)
    #ML.show_cell_time(hid_var)
    cell_path = hid_var['cell_path']
    cell_labels = hid_var['cell_labels']
    cell_time=hid_var['cell_time']
    cell_ori_time=np.array(hid_var['cell_ori_time'])
    n_cell, n_gene = cell_exps.shape
    #n_gene=10
    paths=np.unique(cell_path)
    n_top=300
    n_path=len(paths)
    #print 'gene_names: ',' '.join(gene_names.tolist())
    p_scor_dict={}
    p_gt_dict={}
    p_ge_dict={}
    full_path=defaultdict(lambda:[])
    for p in paths[:]:
        #print '\tpath: ',p
        p_idx=(cell_path==p)
        cell_time_p=np.around(cell_time[p_idx],decimals=2)
        cell_labels_p=cell_labels[p_idx]
        cell_ori_time_p=cell_ori_time[p_idx]
        cell_exps_p=cell_exps[p_idx]
        sort_idx=np.argsort(cell_time_p)
        s_cors=[]
        gt=[]
        ge=[]
        for g_idx in range(n_gene):
            #print 'gene: ',gene_names[g_idx]
            tmp_t=map(lambda x: x+path_info[p]['level']*1.00,cell_time_p[sort_idx])
            tmp_exp = [x for x in cell_exps_p[sort_idx,g_idx]]
            gt.append(tmp_t)
            ge.append(tmp_exp)
            tmp_exp=np.around(tmp_exp,decimals=2)
            tmp_t=np.around(tmp_t,decimals=2)
            scor=scipy.stats.spearmanr(tmp_t,tmp_exp)
            s_cors.append(scor[0])
        p_scor_dict[p]=s_cors
        p_gt_dict[p]=gt
        p_ge_dict[p]=ge

    #print 'paths:', paths
    ori_to_leaf=[]
    for i,p in enumerate(path_info):
        childs=p['child_path']
        if len(childs)>0:
            #print p['ID'],' is not leaf'
            continue
        if p['ID'] not in paths:
            #print p['ID'],' is not in paths'
            continue
        parent=p['parent_path']
        pars=[p['ID']]
        while(parent!=0):
            pars.append(parent)
            path = path_info[parent]
            parent = path['parent_path']
        pars.append(parent)
        ori_to_leaf.append(pars)
    #print ori_to_leaf
    ori_to_leaf_2=[x[:2] for x in ori_to_leaf]
    #print ori_to_leaf_2
    cmd = 'rm -rf '+ out_folder
    print cmd
    os.system(cmd)
    cmd = 'mkdir '+ out_folder
    print cmd
    os.system(cmd)
    
    path_scor_dict={}
    full_ts_es_dict={}

    gene_names_list = map(lambda x: x.lower(), gene_names.tolist())
    #print len(gene_names_list)
    #ana_path = ori_to_leaf
    
    def calculate_scor(fullpath,level):
        fp = fullpath[:level]
        #print fp
        ts=[[] for i in range(n_gene)]
        es=[[] for i in range(n_gene)]
        for i,p in enumerate(fp):
            for i in range(n_gene):
                ts[i]+=p_gt_dict[p][i]
                es[i]+=p_ge_dict[p][i]
        fp_g_scor=[]
        for i in range(n_gene):
            scor=scipy.stats.spearmanr(ts[i],es[i])[0]
            fp_g_scor.append(scor)
        fp_g_scor=np.array(fp_g_scor)
        return fp_g_scor,ts,es
    
    for fp in ori_to_leaf:
        out_pa = '_'.join(map(str,list(reversed(fp))))
        print 'path: ',out_pa                
        fp_g_scor,full_ts,full_es = calculate_scor(fp,len(fp))
        path_scor_dict[tuple(fp)]=fp_g_scor
        full_ts_es_dict[tuple(fp)]=(full_ts,full_es)
        fp_g_scor,_,_ = calculate_scor(fp,1)
        path_scor_dict[tuple(fp[:1])]=fp_g_scor
        fp_g_scor,_,_ = calculate_scor(fp,2)
        path_scor_dict[tuple(fp[:2])]=fp_g_scor
        pos_s_sort_idx=np.argsort(-fp_g_scor)
        neg_s_sort_idx=np.argsort(fp_g_scor)
        pos_genes=gene_names[pos_s_sort_idx[:n_top]]
        neg_genes=gene_names[neg_s_sort_idx[:n_top]]
        mk_gene_dict = marker_genes_dict
        mk_gene_list = marker_gene_list
        if splits[1].split('-')[0]=='treutlein2014':
            mk_gene_dict = treutlein2014_mkgene_dict
            mk_gene_list = treutlein2014_mkgene
        if splits[1].split('-')[0]=='treutlein2016':
            mk_gene_dict = treutlein2016_mkgene_dict
            mk_gene_list = treutlein2016_mkgene
        for index, srt_idx in enumerate(pos_s_sort_idx[:n_top]):
            logging.debug('positive rank: ',index, gene_names[srt_idx], np.around(fp_g_scor[srt_idx],2),mk_gene_dict[gene_names[srt_idx]])
            #print 'positive rank: ',index, gene_names[srt_idx], np.around(fp_g_scor[srt_idx],2),mk_gene_dict[gene_names[srt_idx]]
        #print '\n'
        logging.debug('\n')
        for index, srt_idx in enumerate(neg_s_sort_idx[:n_top]):
            #print 'negative rank: ',index, gene_names[srt_idx], np.around(fp_g_scor[srt_idx],2),mk_gene_dict[gene_names[srt_idx]]
            logging.debug('negative rank: ',index, gene_names[srt_idx], np.around(fp_g_scor[srt_idx],2),mk_gene_dict[gene_names[srt_idx]])

        #head=["query", "significant", "p_value", "T", "Q", "Q&T", "precision", "recall", "term_id","domain", "group", "description", "depth", "intersection", "evcodes"]
        out_file_pos=out_folder+'/'+out_pa+'_pos.txt' 
        out_file_neg=out_folder+'/'+out_pa+'_neg.txt' 
        GO_ana_gene_list(pos_genes,out_file_pos)
        GO_ana_gene_list(neg_genes,out_file_neg)
        
    if model_ana_temp is None:
        model_ana_temp={}
    model_ana_temp['p_gt_dict']=p_gt_dict
    model_ana_temp['p_ge_dict']=p_ge_dict
    model_ana_temp['path_scor_dict']=path_scor_dict
    model_ana_temp['full_ts_es_dict']=full_ts_es_dict
    model_ana_temp['gene_names']=gene_names
    if 'mk_gene_list' not in model_ana_temp.keys():
        model_ana_temp['mk_gene_dict']=mk_gene_dict
        model_ana_temp['mk_gene_list']=mk_gene_list
    model_ana_temp['model']=model
    model_ana_temp['hid_var']=hid_var
    model_ana_temp['out_folder']=out_folder
    return model_ana_temp
def analyze_path_difference(model_ana_temp,path1,path2,append=""):
    #analyze AT1/AT2 diff gene here, path (5,4) and (7,4)
    global treutlein2014_mkgene
    global treutlein2014_mkgene_dict
    print path1
    print path2
    #p_gt_dict=model_ana_temp['p_gt_dict']
    #p_ge_dict=model_ana_temp['p_ge_dict']
    path_scor_dict=model_ana_temp['path_scor_dict']
    #full_ts_es_dict=model_ana_temp['full_ts_es_dict']
    gene_names=model_ana_temp['gene_names']
    out_folder = model_ana_temp['out_folder']

    #AT1_path = (7,4)
    #AT2_path = (5,4)
    abs_scor_diff=np.fabs([path_scor_dict[path1[1]][x]-path_scor_dict[path2[1]][x] for x in range(len(gene_names))])
    #print abs_scor_diff.shape
    srt_index = np.argsort(-abs_scor_diff)
    ntop = 300
    #cutoff = 0.75
    cut_ntop=20
    #print abs_scor_diff[abs_scor_diff>cutoff].shape
    #discovered_marker = []
    #discovered_marker_gene_dict=defaultdict(lambda:"")
    for i,(abs_scor,gene) in enumerate(zip(abs_scor_diff[srt_index],gene_names[srt_index])[:ntop]):
        if i < print_top_limit:
            print abs_scor, gene, treutlein2014_mkgene_dict[gene]
        else:
            logging.debug(str(abs_scor)+" "+ gene +" "+ str(treutlein2014_mkgene_dict[gene]))
        #print abs_scor, gene, marker_genes_dict[gene]
        #if abs_scor>cutoff:
        if i<20:
            model_ana_temp['mk_gene_list']+=[gene]
            val = model_ana_temp['mk_gene_dict'][gene]
            #if model_ana_temp['mk_gene_dict'][gene][0].endswith(path1[0]+'_'+path2[0]+'_scor_diff'):
            #    print 'already added'
            #    continue
            if val!=("",""):
                model_ana_temp['mk_gene_dict'][gene]=(val[0]+'|'+path1[0]+'_'+path2[0]+"_"+append+'_top'+str(i)+'_scor_diff',val[1]+", found by our model")
            else:
                model_ana_temp['mk_gene_dict'][gene]=(path1[0]+'_'+path2[0]+"_"+append+'_top'+str(i)+'_scor_diff', "found by our model")
            print model_ana_temp['mk_gene_dict'][gene][0]
    GO_ana_gene_list(gene_names[srt_index][:ntop],out_folder+'/'+path1[0]+"_"+path2[0]+"_scor_diff.txt")
    model_ana_temp['mk_gene_list']=list(set(model_ana_temp['mk_gene_list']))
    return model_ana_temp
def get_major_cell_types(model_ana_temp):
    #print model_ana_temp['hid_var'].keys()
    cps=model_ana_temp['hid_var']['cell_path']
    cls=model_ana_temp['hid_var']['cell_labels']
    cot=model_ana_temp['hid_var']['cell_ori_time']
    #print cot
    cls_2=np.array(map(lambda x:x[0]+'_16' if x[1]==16 and x[0] =='NA' else x[0], zip(cls,cot)))
    #print cls_2
    ret=[]
    for p in range(max(cps)+1):
        #print p
        #print cls_2[cps==p]
        if cls_2[cps==p].shape[0]==0:
            ret+=['X']
            continue
        lab,cnt=np.unique(cls_2[cps==p],return_counts=True)
        print p, lab
        if 'd5_earlyiN' in lab:
            ret+=['d5_earlyiN']
        elif 'd5_failedReprog' in lab:
            ret+=['d5_failedReprog']
        else:
            ret +=[ lab[cnt==max(cnt)][0]]
    return ret
#get_major_cell_types(model_ana_temp)

def plot_cont_marker_gexp(model_ana_temp,remove_path,rescale=True):
    gene_names=model_ana_temp['gene_names']
    p_gt_dict=model_ana_temp['p_gt_dict']
    p_ge_dict=model_ana_temp['p_ge_dict']
    full_ts_es_dict=model_ana_temp['full_ts_es_dict']
    mk_gene_dict=model_ana_temp['mk_gene_dict']
    mk_gene_list=model_ana_temp['mk_gene_list']
    #p_g_dict=model_ana_temp['p_g_dict']
    model = model_ana_temp['model']
    path_info = model['path_info']
    out_folder = model_ana_temp['out_folder']
    print mk_gene_list
    max_level= float(max([p['level'] for p in path_info]))+1
    path_cell_lab = get_major_cell_types(model_ana_temp)
    print 'mk_gene_list:',mk_gene_list
    for mk_gene in mk_gene_list[:]:
        if mk_gene not in gene_names.tolist():
            continue
        print mk_gene,mk_gene_dict[mk_gene]
        g_idx = gene_names.tolist().index(mk_gene)
        plt.clf()
        fig = plt.figure("all",figsize=(15,10))
        fig2= plt.figure("4th",figsize=(7,5),tight_layout=True)
        
        #plt.subtitle(mk_gene+str(mk_gene_dict[mk_gene]), fontsize="x-large")
        plt1 = fig.add_subplot(2,2,4)
        plt2 = fig.add_subplot(2,2,3)
        plt3 = fig.add_subplot(2,2,1)
        plt4 = fig.add_subplot(2,2,2)
        
        #plt.figure()
        #plt.subplot(131)
        for key,val in full_ts_es_dict.items():
            if key == remove_path:
                continue
            #print key
            #print len(val[0][g_idx])
            #print len(val[1][g_idx])
            x=val[0][g_idx]
            y=val[1][g_idx]
            px=[]
            py=[]
            for p in key:
                x1=path_info[p]['level']
                x2=path_info[p]['level']+1
                y1=model['g_param'][path_info[p]['Sp_idx']][g_idx]
                y2=model['g_param'][path_info[p]['Sc_idx']][g_idx]
                #px+=[x1,x2]
                #py+=[y1,y2]
                px += np.linspace(x1, x2, num=10, endpoint=True).tolist()
                py += np.linspace(y1, y2, num=10, endpoint=True).tolist()
            #px = [path['level']] 
            #py = [model['g_param'][path['Sp_idx']][g_idx]] 
            for tx,ty in zip(x,y):
                if ty==0:
                    continue
                px+=[tx]
                py+=[ty]
            #px+=[path['level']+1]
            #py+=[model['g_param'][path['Sc_idx']][g_idx]]
            #print px
            #print py
            xnew = np.linspace(min(px), max(px), num=100, endpoint=True)
            poly = np.poly1d(np.polyfit(px, py, 4))
            #t,c,k = splrep(px, py, k=3)
            #bspl = BSpline(px, py)
            ynew=poly(xnew)
            #ynew=bspl(xnew)
            #plt.plot(px, py, 'o', xnew, ynew, '-',label = str(key))
            #print key
            if rescale:
                xnew=[_/float(len(key)) for _ in xnew]
            plt.figure('all')
            plt1.plot(xnew, ynew, '-',label = str([_ for _ in reversed(list(key))]))
            plt.figure('4th')
            plt.plot(xnew, ynew, '-',label = str([_ for _ in reversed(list(key))])+' '+path_cell_lab[key[0]])
            #plt1.plot(xnew, ynew, '-',label = str(key))
            #plt3.plot(px, py, 'o',label = str(key))
        plt.figure('all')
        plt.xticks( fontsize = 12)
        plt.yticks( fontsize = 12)
        plt1.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 12)
        plt1.set_title("continuous poly d4 interpolation", fontsize = 20)
#         plt1.set_xticks(fontsize = 12)
#         plt1.set_yticks(fontsize = 12)
        plt.figure('4th')
        plt.xticks(fontsize = 16)
        plt.yticks(fontsize = 16)
#         #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 12)
#         plt.title("continuous poly d4 interpolation", fontsize = 20)
#         plt.xticks(fontsize = 12)
#         plt.yticks(fontsize = 12)
#         plt.xlabel("time",fontsize=15)
#         plt.ylabel("gene expression",fontsize=15)
        plt.figure('all')
        
        
        #plt3.legend()
        
        #plt2.figure()
        for index,path in enumerate(model['path_info']):
            if index not in p_gt_dict.keys():
                continue
            x1, y1 = [path['level'], path['level']+1], [model['g_param'][path['Sp_idx']][g_idx],model['g_param'][path['Sc_idx']][g_idx]] 
            
            #if rescale:
            #    x1=[_/float(path['level']+1) for _ in x1]
            
            plt2.plot(x1, y1, marker = 'o',label=path['ID'])
            #print len(p_gt_dict[index])
            x = np.array(p_gt_dict[index][g_idx])
            y = np.array(p_ge_dict[index][g_idx])
            
            #if rescale:
            #    x=[_/float(path['level']+1) for _ in x]
            
            plt3.plot(x, y, 'o',label=path['ID'])
            
            x_n0_add_line=[]
            y_n0_add_line=[]
            for tx,ty in zip(x,y):
                if ty==0:
                    continue
                x_n0_add_line+=[tx]
                y_n0_add_line+=[ty]
                
            x1=path['level']
            x2=path['level']+1
            y1=model['g_param'][path['Sp_idx']][g_idx]
            y2=model['g_param'][path['Sc_idx']][g_idx]    
            
            x_n0_add_line += np.linspace(x1, x2, num=10, endpoint=True).tolist()
            y_n0_add_line += np.linspace(y1, y2, num=10, endpoint=True).tolist()
            
            #if rescale:
            #    x_n0_add_line=[_/float(path['level']+1) for _ in x_n0_add_line]
            plt4.plot(x_n0_add_line, y_n0_add_line, 'o',label=path['ID'])
            
        plt2.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 12)
        plt3.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 12)
        plt4.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 12)
        plt2.set_title("discrete gene expression", fontsize = 20)
        plt3.set_title("original dataset", fontsize = 20)
        plt4.set_title("modified dataset for interpolation", fontsize = 20)
#         plt2.set_xticks(fontsize = 12)
#         plt2.set_yticks(fontsize = 12)
#         plt3.set_xticks(fontsize = 12)
#         plt3.set_yticks(fontsize = 12)
#         plt4.set_xticks(fontsize = 12)
#         plt4.set_yticks(fontsize = 12)
        
        
        
#         for x in range(0,max([p['level'] for p in path_info])+2):
#             plt1.axvline(x=x, color='k', linestyle=':')
#             plt2.axvline(x=x, color='k', linestyle=':')
#             plt3.axvline(x=x, color='k', linestyle=':')
#             plt4.axvline(x=x, color='k', linestyle=':')

        fig=plt.figure('all')
        fig.subplots_adjust(wspace=0.4)
        fig2=plt.figure('4th')
        fig.savefig(out_folder+'/'+mk_gene+"_"+mk_gene_dict[mk_gene][0]+'_continuous_all.png', bbox_inches='tight')
        fig2.savefig(out_folder+'/'+mk_gene+"_"+mk_gene_dict[mk_gene][0]+'_continuous_4th_nolegend.png', bbox_inches='tight')
        
        plt.figure('4th')
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 12)
        plt.title("continuous poly d4 interpolation", fontsize = 20)
        
        plt.xlabel("time",fontsize=15)
        plt.ylabel("gene expression",fontsize=15)
        fig3=plt.figure('4th')
        fig3.savefig(out_folder+'/'+mk_gene+"_"+mk_gene_dict[mk_gene][0]+'_continuous_4th.png', bbox_inches='tight')
        
        #plt.savefig(out_folder+'/'+mk_gene+"_"+mk_gene_dict[mk_gene][0]+'_continuous_all.png')
        plt.figure('all')
        
        plt.show()
        plt.figure('4th')
        #lt.savefig(out_folder+'/'+mk_gene+"_"+mk_gene_dict[mk_gene][0]+'_continuous_4th.png')
        plt.show()
        
def analyze_model(model_file,dataset='treutlein2014',append="",model_ana_temp=None):
    model_ana_temp=analyze_gene('model/'+model_file,model_ana_temp=model_ana_temp)
    plot_path_fig('model/'+model_file)
    if dataset=='treutlein2014':
        model_ana_temp = analyze_path_difference(model_ana_temp,('AT1',(7,4)),('AT2',(5,4)),append)
        model_ana_temp = analyze_path_difference(model_ana_temp,('ciliated',tuple([2,0])),('Clara',tuple([3,0])),append)
        plot_cont_marker_gexp(model_ana_temp,remove_path=tuple([6,1,0]))
    if dataset=='treutlein2016-2':
        model_ana_temp = analyze_path_difference(model_ana_temp,('Neuron',(8,6)),('d5_failed_intermediate',tuple([4,])),append)
        #model_ana_temp = analyze_path_difference(model_ana_temp,('ciliated',tuple([2,0])),('Clara',tuple([3,0])))
        plot_cont_marker_gexp(model_ana_temp,remove_path=tuple([3,1,0]))
    
    
    return model_ana_temp


def print_continuous_TF(model_file,pcut=0.1, print_no_show=False,print_no_gene=False):
    #out_file = open(model_file+'_tf_lists.txt','w')
    model,hid_var=ML.load_model(model_file)
    path_info=model["path_info"]
    out_df=pd.DataFrame()

    for i,path in enumerate(path_info):
        #print "path: ",path["ID"]," #TF: ", len(path["sorted_tf_tp"]) 
        print "#########################################"

        print "path: ",path["ID"]
        
        print "TF_no_gene: ",len(path["TF_no_gene"])
        print "TF_no_gene_pv: ",len(path["TF_no_gene_pv"])
        print "TF_no_show: ",len(path["TF_no_show"])
        print "TF_no_show_pv: ",len(path["TF_no_show_pv"])
        
        #tmp=path["sorted_tf_tp"]
        #print tmp[0]
        #idx=path["sorted_tf_tp"][:,2]<pcut
        #print idx
        print "top 10 expressed TF"
        count = 0
        app_vec=[]
        for tf,tp,pv,method in path["sorted_tf_tp"]:
            if pv < pcut:
                print "\t", tf,"\t", tp,"\t", pv,"\t", method
                count+=1
                app_vec += [tf.upper()+" "+str(tp)]
        
        while len(app_vec)<10:
            app_vec+=[""]
        out_df[i]=app_vec
        print " # top 10 expressed TF: ", count
        
        count=0
        if print_no_show:
            print "other expressed TF (no show)"
            for tf,pv in zip(path["TF_no_show"],path["TF_no_show_pv"]):
                if pv < pcut:
                    print "\t", tf,"\t", pv
                    count+=1
            print " # no show expressed TF: ", count
        count=0
        if print_no_gene:
            print "TF that is not in gene names"

            for tf,pv in zip(path["TF_no_gene"],path["TF_no_gene_pv"]):
                #print pv,pcut
                if pv < pcut:
                    print "\t", tf,"\t", pv
                    count+=1
            print " # no gene TF: ", count
        
        #print "\n".join(map(str,path["sorted_tf_tp"]))
    return out_df
def plot_cont_TF(model_file,data_file,select=None,searchDB=True,remove_zero=True,title=None):
    #select is the parameter for selecting path and genes of interest, format a list: [[path ids], gene1, gene2, ...]
    model,hid_var=ML.load_model(model_file)
    path_info=model['path_info']
    cell_names,cell_day,cell_labels,cell_exps,gene_names=ML.load_data(data_file,16000)
    cell_path = hid_var['cell_path']
    cell_labels = hid_var['cell_labels']
    cell_time=hid_var['cell_time']
    
    for i,p in enumerate(model['path_info']):
        if select is not None and p["ID"] not in select[0]:
            continue
        print 'path: ',i
        g_pa=model['g_param'][p['Sp_idx']]
        g_pb=model['g_param'][p['Sc_idx']]
        p_idx=cell_path==p["ID"]
        cell_exps_p=cell_exps[p_idx]
        cell_time_p=cell_time[p_idx]
        TF_list=p["TF_list"]
        sorted_tf_tp = p["sorted_tf_tp"]
        TF_tp_dict=defaultdict(lambda: " ")
        for tf,tp,pv,method in sorted_tf_tp:
            TF_tp_dict[tf] = tp
        if select is not None:
            TF_list = select[1:]
            print "selecting genes: ","\t".join(select[1:])
        elif searchDB:
            print "searh TcoF DB interactions"
            search_TcoF_DB(TF_list,out_folder+'/'+model_file.split("/")[1]+"_"+"path_"+str(i))
        for tf_name in TF_list:
            print "TF: ", tf_name
            tf_idx = gene_names.tolist().index(tf_name)
            x=[]
            y=[]
            x1=0
            x2=1
            y1=g_pa[tf_idx]
            y2=g_pb[tf_idx]
            x += np.linspace(x1, x2, num=10, endpoint=True).tolist()
            y += np.linspace(y1, y2, num=10, endpoint=True).tolist()
            for e_p, t_p in zip(cell_exps_p, cell_time_p ):
                if remove_zero and e_p[tf_idx]==0:
                    continue
                x+=[t_p]
                y+=[e_p[tf_idx]]
                
            xnew = np.linspace(min(x), max(x), num=100, endpoint=True)
            poly = np.poly1d(np.polyfit(x, y, 4))
            ynew=poly(xnew)
            plt.figure('all')
            plt.plot(xnew, ynew, '-',label = tf_name.upper()+" "+str(TF_tp_dict[tf_name]))
        
        plt.figure('all')    
        plt.xticks(fontsize = 16)
        plt.yticks(fontsize = 16)
        #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 12)
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.4),ncol=3, fancybox=True, shadow=True, prop={'size': 14})
        if title is None:
            plt.title("path "+str(i)+" all TF "+"poly d4 interpolation", fontsize = 20)
        else:    
            plt.title(title, fontsize = 20)
        plt.xlabel("time",fontsize=15)
        plt.ylabel("gene expression",fontsize=15)
        #ax = plt.figure().gca()
        plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
        
        
        
        fig=plt.figure('all')
        if select is None:
            fig.savefig(data_file+"_"+"path_"+str(i)+"_allTF"+'_continuous_plot.png', bbox_inches='tight')
        else:
            fig.savefig(data_file+"_"+"path_"+str(i)+"_"+"-".join(select[1:])+'_continuous_plot.png', bbox_inches='tight')
            
        plt.figure('all')
        plt.show()     
        plt.close()

def search_TcoF_DB(TF_list,out_file=None,verbose=True):
    df=pd.read_csv("TcoF_DBv2_human_mouse_tf_interactions.csv")
    #print df["Symbol1"]
    #print df["Symbol2"]
    ret=[]
    comb=[]
    for i in range(len(TF_list)-1):
        for j in range(i+1,len(TF_list)):
            comb+=[(i,j)]
            TF1=TF_list[i]
            TF2=TF_list[j]
            #The DB file is symmetric
            df2 = df[df["Symbol1"]==TF1.upper()  ]
            out_df = df2[df2["Symbol2"]==TF2.upper()]
            if out_df.shape[0]>0:
                ret+=[out_df]

                if verbose:
                    print TF1, TF2, out_df.shape[0] ,"interactions"
                if out_file is not None:
                    out_df.to_csv(out_file+"_TF_"+TF1+"_"+TF2+".csv",index=False)
    print "#combinations: ",len(comb)
    return ret
def calc_TcoF_DB_pval(model_file,data_file,out_folder='TF_analysis_result',calc_all_interaction=True):
    model,hid_var=ML.load_model(model_file)
    path_info=model['path_info']
    cell_names,cell_day,cell_labels,cell_exps,gene_names=ML.load_data(data_file,16000)
    #cell_path = hid_var['cell_path']
    #cell_labels = hid_var['cell_labels']
    #cell_time=hid_var['cell_time']
    dTD=ML.getdTD("tfDNA_predicted_100.txt.update")
    tmp = dTD.keys()
    all_tfs=[]
    for x in tmp:
        if x.lower() in gene_names:
            all_tfs+=[x.lower()]
    print "# of all tfs in gene names: ",len(all_tfs)
    #print all_tfs
    if calc_all_interaction:
        print "searh TcoF DB interactions for all tfs"
        all_tf_interactions=search_TcoF_DB(all_tfs,out_file=None, verbose=False)
        print "#interaction of all tfs: ", len(all_tf_interactions)
    sum_path_interaction=0
    sum_path_interaction_early=0
    sum_path_interaction_late=0
    sum_path_interaction_mix=0
    
    sum_combination_all=0
    sum_combination_early=0
    sum_combination_late=0
    sum_combination_mix=0
    for i,p in enumerate(model['path_info']):
        print 'path: ',i
        TF_list=p["TF_list"]
        sorted_tf_tp=p["sorted_tf_tp"]
        early_tfs=[]
        late_tfs=[]
        for (tf_name,tp,pv,method) in sorted_tf_tp:
            if tp == 0:
                early_tfs+=[tf_name]
            else:
                late_tfs+=[tf_name]
        print "early tfs: ",early_tfs
        print "late tfs: ",late_tfs
        print "searh TcoF DB all interactions"
        pi=search_TcoF_DB(TF_list,out_file=out_folder+"/path_"+str(i),verbose=None)
        print "all interaction: ",len(pi)
        sum_path_interaction+=len(pi)
        sum_combination_all+=len(TF_list)*(len(TF_list)-1)/2
        print "searh TcoF DB early interactions"
        pi_early=search_TcoF_DB(early_tfs,out_file=None,verbose=None)
        print "early interaction: ",len(pi_early)
        sum_path_interaction_early+=len(pi_early)
        sum_combination_early+=len(early_tfs)*(len(early_tfs)-1)/2
        print "searh TcoF DB late interactions"
        pi_late=search_TcoF_DB(late_tfs,out_file=None,verbose=None)
        print "late interaction: ",len(pi_late)
        sum_path_interaction_late+=len(pi_late)
        sum_combination_late+=len(late_tfs)*(len(late_tfs)-1)/2
        
        sum_path_interaction_mix+=len(pi)-len(pi_early)-len(pi_late)
        sum_combination_mix+=len(early_tfs)*len(late_tfs)
        
    print "sum of # of all interactions: ",sum_path_interaction,'/',sum_combination_all, '=',sum_path_interaction/float(sum_combination_all)
    print "sum of # of early interactions: ",sum_path_interaction_early,'/',sum_combination_early, '=',sum_path_interaction_early/float(sum_combination_early)
    print "sum of # of late interactions: ",sum_path_interaction_late,'/',sum_combination_late, '=',sum_path_interaction_late/float(sum_combination_late)
    print "sum of # of mix interactions: ",sum_path_interaction_mix,'/',sum_combination_mix, '=',sum_path_interaction_mix/float(sum_combination_mix)

            