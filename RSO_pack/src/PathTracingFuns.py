# -*- coding: utf-8 -*-
import networkx as nx

def findRelaysInPath(path,ProDevices):
    sPathRs = [None] * (len(path)-1)
    rEdges = [(x['Bus1'],x['Bus2']) for x in ProDevices]
    for ii in range(len(path)-1):
        if((path[ii],path[ii+1]) in rEdges):
            sPathRs[ii] = rEdges.index((path[ii],path[ii+1]))
    return(sPathRs)    

def findDeviceInPath(path,ProDevices,DOC):
    sPathRs = [None] * (len(path)-1)
    rEdges = [(x['Bus1'],x['Bus2']) for x in ProDevices]
    for ii in range(len(path)-1):
        if(DOC==1):
            if((path[ii],path[ii+1]) in rEdges):
                sPathRs[ii] = rEdges.index((path[ii],path[ii+1]))
        else:
            if((path[ii],path[ii+1]) in rEdges):
                sPathRs[ii] = rEdges.index((path[ii],path[ii+1]))
            elif((path[ii+1],path[ii]) in rEdges):
                sPathRs[ii] = rEdges.index((path[ii+1],path[ii]))
    return sPathRs    


def faultPathDir(isProDevice,ProDevices):
    R = list(set([x for xx in isProDevice for x in xx if x is not None]))
    FD = [0]*len(ProDevices)
    for ii in range(len(R)):
        FD[R[ii]] = 1
    return FD 

                

def create_priamry_backup_from_paths(isProDevice,ProDevices):
    indx_path_relays = [x for x in isProDevice if x is not None]
    
    if(len(indx_path_relays)>1):
        kk=0
        pri_bac=[[]]*(len(indx_path_relays)-1)
        #bac=[[]]*(len(indx_path_relays)-1)
        #Pri_bac_type = [[]]*(len(indx_path_relays)-1)
        #bac_type = [[]]*(len(indx_path_relays)-1)
        for ii in range(len(indx_path_relays)-1,0,-1):
            pri_bac[kk] = [ProDevices[indx_path_relays[ii]]['Bus1'],ProDevices[indx_path_relays[ii]]['Bus2'],\
                           ProDevices[indx_path_relays[ii-1]]['Bus1'],ProDevices[indx_path_relays[ii-1]]['Bus2'],\
                           ProDevices[indx_path_relays[ii]]['Type'],ProDevices[indx_path_relays[ii-1]]['Type'],None,None]
            #bac[kk] = ()
            #Pri_bac_type[kk] = (ProDevices[indx_path_relays[ii]]['Type'],ProDevices[indx_path_relays[ii-1]]['Type'])
            #bac_type[kk] = ProDevices[indx_path_relays[ii-1]]['Type']
            kk=kk+1
        return pri_bac
    else:
        return []

        
def find_edgenode(G,ProDevices,Substation_bus):
    edgeNodeList = [n for (n,d) in G.degree if d==1]
    priRerlaysInPath = [None] * len(edgeNodeList)
    nodeDistance = [None] * len(edgeNodeList)
    for ii in range(len(edgeNodeList)):
        if(nx.has_path(G,Substation_bus,edgeNodeList[ii])):
            path = nx.shortest_path(G,Substation_bus,edgeNodeList[ii])
        else:
            path = []
        # find_relays_inpath(path,Pro_Devices);
        rInPath = findRelaysInPath(path,ProDevices)
        if(not rInPath == [None] * len(rInPath)):
            priRerlaysInPath[ii] = [x for x in rInPath if x is not None][-1]
            nodeDistance[ii] = len(path) - (rInPath.index([x for x in rInPath if x is not None][-1])+1)
    kk = 1
    fault_node = [None] * len(ProDevices)
    for ii in range(len(ProDevices)):
        fnode_list = [edgeNodeList[x] for x in range(len(priRerlaysInPath)) if priRerlaysInPath[x]==ii]
        fnode_dist = [nodeDistance[x] for x in range(len(priRerlaysInPath)) if priRerlaysInPath[x]==ii]
        if(len(fnode_list)!=0):
            fnode_ind = fnode_dist.index(max(fnode_dist))
            fault_node[kk] = fnode_list[fnode_ind]
            kk=kk+1
    
    return [x for x in fault_node if x is not None]

def find_faultpath_insys(fp,sb,G):
    sPaths =[[]]*len(sb)
    for ii in range(len(sb)):
        if(nx.has_path(G,sb[ii],fp)):
            sPaths[ii] = nx.shortest_path(G,sb[ii],fp)
        else:
            sPaths[ii] = []
            
    return sPaths        

def isPrimaryRelay(rPri,fPoint,RList,G):
    
    fault_path = nx.shortest_path(G,rPri,fPoint)            
    isPri = True
    BackupNum = 0
    
    for ii in range(len(RList)):
        pN1 = RList[ii]['Bus1']
        pN2 = RList[ii]['Bus2']
        
        pN1_match = pN1 in fault_path
        pN2_match = pN2 in fault_path
        
        if(pN1_match):
            pN1_ind = fault_path.index(pN1)
        else:
            pN1_ind = None
            
        if(pN2_match):
            pN2_ind = fault_path.index(pN2)
        else:
            pN2_ind = None
            
        if(pN1_match and pN2_match):
            if(pN1_ind+1 == pN2_ind):
                isPri = False
                BackupNum = BackupNum+1
                
    return (isPri,BackupNum)
    



