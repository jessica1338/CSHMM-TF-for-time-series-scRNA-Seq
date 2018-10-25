import scdiff.scdiff as S
import argparse
def run_scdiff_init(data_file,tfdna=None,config="auto",large=None):
#     E=S.TabFile(data_file).read('\t')
#     print E[0][:3]
#     #global GL # Gene list global variable
#     S.GL=E[0][3:]
#     S.dTD={}
#     S.dMb={}
#     S.dTG={}
#     E=E[1:]
#     AllCells=[]
#     #S.SPECTRALIMIT=20
#     for i in E:
#             iid=i[0]
#             ti=float(i[1])     # time point for cell i
#             li=i[2]
#             ei=[float(item) for item in i[3:]] # expression for cell i
#             ci=S.Cell(iid,ti,ei,li) # cell i
#             AllCells.append(ci)
    #pdb.set_trace()
    #-----------------------------------------------------------------------
    # 1) : read in gene expression
    AllCells=[]
    print("reading cells...")
    with open(data_file,'r') as f:
        line_ct=0
        for line in f:
            #print line
            if line_ct==0:
                S.GL=line.strip().split("\t")[3:]
            else:
                line=line.strip().split("\t")
                iid=line[0]
                ti=float(line[1])
                li=line[2]
                ei=[round(float(item),2) for item in line[3:]]
                ci=S.Cell(iid,ti,ei,li,S.GL)
                AllCells.append(ci)
            line_ct+=1
            print('cell:'+str(line_ct))
            
    firstTime=min([float(item.T) for item in AllCells])
    
    
    G1=S.Graph(AllCells,tfdna,config,large,None,None)

#     if large:
#         G1=S.Graph(AllCells,'auto','True')  #Cells: List of Cell instances 
#     else:
#         G1=S.Graph(AllCells,'auto',None)  #Cells: List of Cell instances 
    out_file = open('init_cluster_'+data_file+'.txt','w')
    pairs=[]
    pairs2=[]
    print 'G1'
    for node in G1.Nodes:
        print 'ID: ', node.ID
        print 'index: ',G1.Nodes.index(node)
        ic=G1.Nodes.index(node)
        if node.P is not None:
            print 'P index: ',G1.Nodes.index(node.P)
            ip=G1.Nodes.index(node.P)
            pairs.append(str(ip)+' '+str(ic))
        else:
            print 'P index: ',None
        print 'ST: ',node.ST
        print 'T: ',node.T
        for cell in node.cells:
            print 'cell.ID: ', cell.ID
            pairs2.append(cell.ID+" "+str(ic))
    print pairs
    print pairs2
    out_file.write('\t'.join(pairs)+'\n')
    out_file.write('\t'.join(pairs2)+'\n')
    out_file.close()
if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-tf',"--tf_dna_file", help="specify the tf_target file")
    parser.add_argument('-d',"--data_file", help="specify the data file")
    parser.add_argument('-l',"--large_dataset", help="specify whether or not to speed up the initialization for large dataset (a different algorithm (PCA+k-means) will be applied)",type=int,default=0)
    args=parser.parse_args()

    large=None
    if args.large_dataset>1:
        large="True"
    if args.data_file is not None:
        data_file=args.data_file
        tfdna=args.tf_dna_file
        run_scdiff_init(data_file,tfdna,large=None) #if it takes too long, set large="True" for large dataset
    else:
        print "please use -d to specify the data file to run"
        print "please use -tf to specify the tf-target file to run"

