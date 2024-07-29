
import os
import pygad
import random
import networkx as nx
import matplotlib.pyplot as plt

from copy import deepcopy
from .on_gen import on_generation
from .OT_funs import OCTCC_Name
from .GenNxGraph import genGraph, index_dict
from .write_fit_fun import write_GA_fit_fun, write_Con 
from .PathTracingFuns import faultPathDir, isPrimaryRelay
from .Read_CSV_Functions import read_Relays_CSV_Data, read_Fault_CSV_Data, pdef
from .PathTracingFuns import find_faultpath_insys, findDeviceInPath, create_priamry_backup_from_paths
from .DirCalc import calcZ1Z0ANG

def runSettingsOptimizer(Main_dir,switchStates,switchLines,Device_Data_CSV,Fault_Data_CSV,Fault_Res,SysInfo_IN,Substation_bus,Min_Ip,enableIT,Force_NOIBR,CTI,OTmax,type_select=False,SetDir=False,Sho_Plots=True,GA_initial_seed=None):
    """
    
    Parameters
    ----------
    Main_dir : Str
        Main Dir path.
    Device_Data_CSV : List
        List of Devices and ther load values.
    Fault_Data_File : str
        Fault Data File location.
    SysInfo : Json 
        Json struct containg system data.
    Substation_bus : str
        Name of substation bus.
    enableIT : int
        1 - enable 50p/g, 0-dissable 50p/g.
    Force_NOIBR : int
        1 - No Dir settings, 0 - with Dir settings
    CTI : float
        Coordination time interval.
    OTmax : float
        maximum operating time.
    Sho_Plots : bool
        enable debuging plots.
    GA_initial_seed : List
        Previous solution to be used as the seed for the optimization

    Returns
    -------
    Relay_settings : list of dicts
        Optimized settings.
    info : Optimization Info for debug
        .

    """
    # %% Read Fault input Data
    fB = list(set([x['busNumber'] for x in Fault_Data_CSV]))
    fB.sort()
    faultBuses = [None]*len(fB)
    faultBusesLOC = [None]*len(fB)
    for ii in range(len(fB)):
        if('_closein' in fB[ii]):
            faultBuses[ii] = fB[ii].split('_closein')[0]
            faultBusesLOC[ii] = fB[ii].split('_closein')[1]
        else:
            faultBuses[ii] = fB[ii]
            faultBusesLOC[ii] = fB[ii]
    
    # %% Read system data
    SysInfo = deepcopy(SysInfo_IN)
    Relays = SysInfo['Relays']
    Recs = SysInfo['Recs']
    Fuses = SysInfo['Fuses']
    Lines = SysInfo['Lines']
    XFMRs = SysInfo['XFMRs']
    Buses = SysInfo['Buses']
    Pvs = SysInfo['Pvs']
    BESS = SysInfo['BESS']
    Gens = SysInfo['Gens']
    
    
    if(Force_NOIBR == 1):
        for ii in range(len(Device_Data_CSV)):
            Device_Data_CSV[ii]['PVs'] = '0'
    
    
    # %%Edit system based on settings
    for ii in range(len(Pvs)):
        if(Device_Data_CSV[ii]['PVs'] == '0'):
            Pvs[ii]['Enabled'] = False    
    
        
    if(sum([x['Enabled'] for x in Pvs] + [x['Enabled'] for x in BESS] + [x['Enabled'] for x in Gens])<=0):
        DOC = 0
    else:
        DOC = 1   
    
    # update switch states
    for ii in range(len(switchLines)):
        if(switchStates[ii]==1):
            #dssText.Command = 'edit line.'+switchLines[ii]+' enabled=true'
            Lines[index_dict(Lines,'Name',switchLines[ii])]['Enabled'] = True
        elif(switchStates[ii]==0):
            #dssText.Command = 'edit line.'+switchLines[ii]+' enabled=false'
            Lines[index_dict(Lines,'Name',switchLines[ii])]['Enabled'] = False
    
    # update Devices or Lines 
    for ii in range(len(Relays)):
        if(Relays[ii]['Name'].lower() in [x['RelayName'].lower() for x in Device_Data_CSV]):
            Relays[ii]['Enabled'] = Device_Data_CSV[index_dict(Device_Data_CSV,'RelayName',Relays[ii]['Name'])]['switchState']
    for ii in range(len(Recs)):
        if(Recs[ii]['Name'].lower() in [x['RelayName'].lower() for x in Device_Data_CSV]):
            Recs[ii]['Enabled'] = Device_Data_CSV[index_dict(Device_Data_CSV,'RelayName',Recs[ii]['Name'])]['switchState']
    # %% Generate Net Graph        
    Gr = genGraph(Lines,XFMRs,Buses,Relays,Recs,Fuses,Pvs,BESS,Gens,switchLines,switchStates)
    G = Gr[0]
    pos = Gr[1]
    Edges = Gr[2]
    
    #Node_list = list(G.nodes)
    #Edge_list = list(G.edges)
    if(Substation_bus in G.nodes):
        G.nodes[Substation_bus]['isSource'] = True
        G.nodes[Substation_bus]['sourceName'] = 'Sub'
        G.nodes[Substation_bus]['color'] = 'r'
        
    #Plot Graph
    Ecolors = [G[u][v]['color'] for u,v in G.edges]
    Ewidth = [G[u][v]['width'] for u,v in G.edges]
    Ncolors = [G.nodes[u]['color'] for u in G.nodes]
    Nsize = [G.nodes[u]['size'] for u in G.nodes]
    if(Sho_Plots):
        nx.draw(G,pos=pos,with_labels=1,edge_color=Ecolors,width = Ewidth,node_size=Nsize,node_color = Ncolors)
    
    # %% Find Source Buses
    kk=0
    sourceBuses = ['None'] * (len(Pvs)+len(BESS)+len(Gens)+1)
    for node in G.nodes:
        if(G.nodes[node]['isPV'] or G.nodes[node]['isBESS'] or G.nodes[node]['isSource'] ):
            sourceBuses[kk] = G.nodes[node]['Name']
            kk=kk+1
    del sourceBuses[kk:]
    
    # %% Generate Protective Device List
    
    # Relays
    nRelays = sum([nx.get_edge_attributes(G,'isRelay')[x]==True for x in nx.get_edge_attributes(G,'isRelay')])
    nRec = sum([nx.get_edge_attributes(G,'isRecloser')[x]==True for x in nx.get_edge_attributes(G,'isRecloser')]) 
    nFuse = sum([nx.get_edge_attributes(G,'isFuse')[x]==True for x in nx.get_edge_attributes(G,'isFuse')])
    
    if(DOC == 1):
        nRelays = nRelays*2
        nRec = nRec*2
    
    nProDevices = nRelays+nRec+nFuse
    print("Relays : %d, Recs: %d, Fuses: %d,Total Devices: %d"  % (nRelays,nRec,nFuse,nProDevices))
    
    if(nProDevices == 0):
        print('No deives to set')
        return [],[]
    ProDevices = [dict.fromkeys(['Bus1','Bus2','phases','Line','Direction','Name','Type','Vpu','IL','In','Ip','Inp','Imax3ph','Imin3ph','Vmax3ph','Vmin3ph','Igmax3ph','Igmin3ph',\
                  'ImaxLL','IminLL','VmaxLL','VminLL','IgmaxLL','IgminLL','ImaxSLG','IminSLG','VmaxSLG','VminSLG','IgmaxSLG','IgminSLG','IT','ITg','Oind']) for number in range(nProDevices)]
    
    kk = 0;
    for ii in range(len(Relays)):
        # Relays
        if(Relays[ii]['MonitoredObj'].split('.')[1].lower() in [Name[2] for Name in list(G.edges.data('Name'))] and Relays[ii]['Enabled']):
            ProDevices[kk]['Bus1'] = Relays[ii]['Bus1']
            ProDevices[kk]['Bus2'] = Relays[ii]['Bus2']
            ProDevices[kk]['Line'] = [Edges[x]['Name'] for x in range(len(Edges))].index(Relays[ii]['MonitoredObj'].split('.')[1].lower())
            ProDevices[kk]['phases'] = Edges[ProDevices[kk]['Line']]['numPhases']
            if DOC==1:
                ProDevices[kk]['Direction'] = 'F' 
                ProDevices[kk]['Name'] = Relays[ii]['Name']
                ProDevices[kk]['Type'] = 'Relay_TOC' #add Logic later
                kk=kk+1
                
                ProDevices[kk]['Bus1'] = ProDevices[kk-1]['Bus2']
                ProDevices[kk]['Bus2'] = ProDevices[kk-1]['Bus1']
                ProDevices[kk]['Line'] = ProDevices[kk-1]['Line']
                ProDevices[kk]['phases'] = ProDevices[kk-1]['phases']
                ProDevices[kk]['Direction'] = 'R' 
                ProDevices[kk]['Name'] = Relays[ii]['Name']
                ProDevices[kk]['Type'] = 'Relay_TOC' #add Logic later
                kk=kk+1
            else:
                ProDevices[kk]['Direction'] = 'E' #add logic later
                ProDevices[kk]['Name'] = Relays[ii]['Name']
                ProDevices[kk]['Type'] = 'Relay_TOC' #add Logic later
                kk=kk+1
        
        
    # Rec
    for ii in range(len(Recs)):
        if(Recs[ii]['MonitoredObj'].split('.')[1].lower() in [Name[2] for Name in list(G.edges.data('Name'))] and Recs[ii]['Enabled']):
            ProDevices[kk]['Bus1'] = Recs[ii]['Bus1']
            ProDevices[kk]['Bus2'] = Recs[ii]['Bus2']
            ProDevices[kk]['Line'] = [Edges[x]['Name'] for x in range(len(Edges))].index(Recs[ii]['MonitoredObj'].split('.')[1].lower())
            ProDevices[kk]['phases'] = Edges[ProDevices[kk]['Line']]['numPhases']
            if DOC==1:
                ProDevices[kk]['Direction'] = 'F' 
                ProDevices[kk]['Name'] = Recs[ii]['Name']
                ProDevices[kk]['Type'] = 'Rec_TOC' #add Logic later
                kk=kk+1
                
                ProDevices[kk]['Bus1'] = ProDevices[kk-1]['Bus2']
                ProDevices[kk]['Bus2'] = ProDevices[kk-1]['Bus1']
                ProDevices[kk]['Line'] = ProDevices[kk-1]['Line']
                ProDevices[kk]['phases'] = ProDevices[kk-1]['phases']
                ProDevices[kk]['Direction'] = 'R' 
                ProDevices[kk]['Name'] = Recs[ii]['Name']
                ProDevices[kk]['Type'] = 'Rec_TOC' #add Logic later
                kk=kk+1
            else:
                ProDevices[kk]['Direction'] = 'E' #add logic later
                ProDevices[kk]['Name'] = Recs[ii]['Name']
                ProDevices[kk]['Type'] = 'Rec_TOC' #add Logic later
                kk=kk+1
    # Fuse
    for ii in range(len(Fuses)):
        if(Fuses[ii]['MonitoredObj'].lower().split('.')[1] in [Name[2] for Name in list(G.edges.data('Name'))]):
            ProDevices[kk]['Bus1'] = Fuses[ii]['Bus1']
            ProDevices[kk]['Bus2'] = Fuses[ii]['Bus2']
            ProDevices[kk]['Line'] = [Edges[x]['Name'] for x in range(len(Edges))].index(Fuses[ii]['MonitoredObj'].lower().split('.')[1])
            ProDevices[kk]['phases'] = Edges[ProDevices[kk]['Line']]['numPhases']
            ProDevices[kk]['Direction'] = 'E' #add logic later
            ProDevices[kk]['Name'] = Fuses[ii]['Name']
            ProDevices[kk]['Type'] = Fuses[ii]['Type'] #add Logic later
            kk=kk+1
            
    # %% Generate or Get Fault bus Locations
    # rfaultBuses = list(set([x['Bus1'] for x in ProDevices]+[x['Bus2'] for x in ProDevices]))        
    # rfaultNodes = find_edgenode(G,ProDevices,Substation_bus) 
    # faultBuses = list(set(rfaultBuses+rfaultNodes))
    
    faultBusPhases = [None]*len(faultBuses)
    for ii in range(len(faultBuses)):
        faultBusPhases[ii] = Buses[index_dict(Buses,'Name',faultBuses[ii])]['nodes']
        
    # %% Update Load voltage and currents 
    for ii in range(len(ProDevices)):
        #Linfo = getLineVI(dssCircuit,Edges[ProDevices[ii]['Line']]['Name'])
        Dev_ind = index_dict(Device_Data_CSV,'RelayName',ProDevices[ii]['Name'])
        LineRealPower = Device_Data_CSV[Dev_ind]['P']
        if(ProDevices[ii]['Direction'] == 'E'):
            Idir = 1;
        else:
            if(LineRealPower>0 and ProDevices[ii]['Direction'] == 'F'):
                Idir = 1
            elif(LineRealPower>0 and ProDevices[ii]['Direction'] == 'R'):
                Idir = 0
            elif(LineRealPower<0 and ProDevices[ii]['Direction'] == 'F'):    
                Idir = 0
            elif(LineRealPower<0 and ProDevices[ii]['Direction'] == 'R'):     
                Idir = 1
        ProDevices[ii]['Vpu'] = Device_Data_CSV[Dev_ind]['Va_mag']/(4.16e3/(3**0.5)) #min([x for x in Linfo[0] if x!=0])
        ProDevices[ii]['IL'] = Device_Data_CSV[Dev_ind]['Ia_mag'] * Idir # max(Linfo[1])*Idir
        ProDevices[ii]['In'] = Device_Data_CSV[Dev_ind]['In'] * Idir #Linfo[2]*Idir
        ProDevices[ii]['Vabc'] = Device_Data_CSV[Dev_ind]['Vabc']
        ProDevices[ii]['Iabc'] = Device_Data_CSV[Dev_ind]['Iabc']
        ProDevices[ii]['V012'] = Device_Data_CSV[Dev_ind]['V012']
        ProDevices[ii]['I012'] = Device_Data_CSV[Dev_ind]['I012']
        ProDevices[ii]['Z012'] = Device_Data_CSV[Dev_ind]['Z012']
        if(DOC == 1 and ProDevices[ii]['Direction'] == 'F' and SetDir):
            Z1,Z0,Z2,Z12T = calcZ1Z0ANG(ProDevices[ii]['Name'],Device_Data_CSV,Fault_Data_CSV,G,ProDevices[ii]['Bus1'],ProDevices[ii]['Bus2'],0)
            ProDevices[ii]['Z1MAG'] = Z1[0]
            ProDevices[ii]['Z1ANG'] = Z1[1]
            ProDevices[ii]['Z0MAG'] = Z0[0]
            ProDevices[ii]['Z0ANG'] = Z0[1]
            ProDevices[ii]['Z2F'] = Z2[0]
            ProDevices[ii]['Z2R'] = Z2[1]
            #print (ProDevices[ii]['Name']+' --- Z1: '+ str(ProDevices[ii]['Z1ANG'])+' Z0: '+str(ProDevices[ii]['Z0ANG'])+' Z2F: '+str(ProDevices[ii]['Z2F'])+' Z2R: '+str(ProDevices[ii]['Z2R']))
        elif(DOC == 1 and ProDevices[ii]['Direction'] == 'R' and SetDir):
            # find F and copy
            ProDevices[ii]['Z1MAG'] = ProDevices[ii-1]['Z1MAG']
            ProDevices[ii]['Z1ANG'] = ProDevices[ii-1]['Z1ANG']
            ProDevices[ii]['Z0MAG'] = ProDevices[ii-1]['Z0MAG']
            ProDevices[ii]['Z0ANG'] = ProDevices[ii-1]['Z0ANG']
            ProDevices[ii]['Z2F'] = ProDevices[ii-1]['Z2F']
            ProDevices[ii]['Z2R'] = ProDevices[ii-1]['Z2R']
            #print (ProDevices[ii]['Name']+' --- Z1: '+ str(ProDevices[ii]['Z1ANG'])+' Z0: '+str(ProDevices[ii]['Z0ANG'])+' Z2F: '+str(ProDevices[ii]['Z2F'])+' Z2R: '+str(ProDevices[ii]['Z2R']))
        else:
            pass # no need for Dir 

    # %% Path Tracing for Coordiantion pairs   
    faultPaths = [[]]*len(faultBuses) 
    faultDir = [[]]*len(faultBuses)
    
    isProDevice = [[[]]*len(sourceBuses) for i in range(len(faultBuses))]
    
    #isProDevice = [[[]]*len(sourceBuses)]*len(faultBuses) 
    
    pri_bac = [[]]*(2**len(ProDevices)) 
    #bac = [[]]*(2**len(ProDevices))
    
    # ii = faultbus
    # jj = sourceBuses
    # kk = Node in path
    Pairs_len = 0
    for ii in range(len(faultBuses)):
        faultPaths[ii] = find_faultpath_insys(faultBuses[ii],sourceBuses,G)
        for jj in range(len(faultPaths[ii])):
            if(len(faultPaths[ii][jj])!=0):
                #print((ii,jj))
                isProDevice[ii][jj] = findDeviceInPath(faultPaths[ii][jj],ProDevices,DOC)
                Temp = create_priamry_backup_from_paths(isProDevice[ii][jj],ProDevices)
                if(len(Temp)>0):
                    for pp in range(len(Temp)):
                        if(Temp[pp] not in pri_bac):
                            #print(Temp[pp])
                            pri_bac[Pairs_len] = Temp[pp]
                            Pairs_len = Pairs_len+1
                            
        faultDir[ii] = faultPathDir(isProDevice[ii],ProDevices)  
    
    del pri_bac[Pairs_len:]
    Pairs = sorted(pri_bac,key=lambda y: y[0])
    
    if(len(Pairs)==0):
        print('No Coordination Pairs Found in the System, ensure at least 1 fault in each relay protection zone')
    
    # %% Build M2
    # M2: [1    , 2     , 3   , 4   , 5   , 6   , 7    , 8    , 9     , 10    , 11    , 12    , 13     , 14     ,15     ,16      ,17      , 18     ,19       , 20      ,21       ,22     ]
    # M2: [Fault, Device, Imax, Imin, Vmax, Vmin, I0max, I0min, IMaxLL, IminLL, VMaxLL, VminLL, I0MaxLL, I0MinLL,IMaxSLG, IminSLG, VMaxSLG, VminSLG, I0MaxSLG, I0MinSLG,isPrimary,backunm]
    M2 = [[None]*22 for i in range( (len(ProDevices)*len(faultBuses)*6) )]
    
    kk=0;
    for ii in range(len(faultBuses)):
        DevIds = list(set([x for xs in isProDevice[ii] for x in xs if x != None]))
        if(len(DevIds)==0):
            continue
        
        for rr in DevIds:
            M2[kk][0] = ii # faultBus number
            M2[kk][1] = rr # device Number
            
            # determin fault bus name
            if(faultBuses[ii] == faultBusesLOC[ii]):
                fault_bus_name = faultBuses[ii]
            else:
                fault_bus_name = faultBuses[ii]+'_closein'+faultBusesLOC[ii]     
            
            str_r1 = Fault_Res[0]
            
            if(len(Fault_Res)==1):
                str_r2 = Fault_Res[0] # if only 1 resistance avaliable use it for min and max
            elif(len(Fault_Res)>1):
                str_r2 = Fault_Res[-1] # if multiple resistacne avaliable use first and last ( assumed  inscresing order)
            
            if(len(faultBusPhases[ii])>2):
                str_t = 'TPH_'
                FDI1 = [Fault_Data_CSV.index(Fault) for Fault in Fault_Data_CSV if 
                                  (Fault['busNumber'] == fault_bus_name and 
                                   Fault['FaultType'] == (str_t+str_r1) and 
                                   Fault['Relay'].lower() == ProDevices[rr]['Name'])][0]
                
                FDI2 = [Fault_Data_CSV.index(Fault) for Fault in Fault_Data_CSV if 
                                  (Fault['busNumber'] == fault_bus_name and 
                                   Fault['FaultType'] == (str_t+str_r2) and 
                                   Fault['Relay'].lower() == ProDevices[rr]['Name'])][0]

                M2[kk][2] = max(Fault_Data_CSV[FDI1]['Ia_mag'],Fault_Data_CSV[FDI1]['Ib_mag'],Fault_Data_CSV[FDI1]['Ic_mag'])
                M2[kk][3] = max(Fault_Data_CSV[FDI2]['Ia_mag'],Fault_Data_CSV[FDI2]['Ib_mag'],Fault_Data_CSV[FDI2]['Ic_mag'])
                
                M2[kk][4] = min(Fault_Data_CSV[FDI1]['Va_mag'],Fault_Data_CSV[FDI1]['Vb_mag'],Fault_Data_CSV[FDI1]['Vc_mag'])
                M2[kk][5] = min(Fault_Data_CSV[FDI2]['Va_mag'],Fault_Data_CSV[FDI2]['Vb_mag'],Fault_Data_CSV[FDI2]['Vc_mag'])

                M2[kk][6] = Fault_Data_CSV[FDI1]['In_mag']
                M2[kk][7] = Fault_Data_CSV[FDI2]['In_mag']

            if(len(faultBusPhases[ii])>1):
               str_t = 'BC_'
               FDI1 = [Fault_Data_CSV.index(Fault) for Fault in Fault_Data_CSV if 
                                   (Fault['busNumber'] == fault_bus_name and 
                                    Fault['FaultType'] == (str_t+str_r1) and 
                                    Fault['Relay'].lower() == ProDevices[rr]['Name'])][0]
                 
               FDI2 = [Fault_Data_CSV.index(Fault) for Fault in Fault_Data_CSV if 
                                   (Fault['busNumber'] == fault_bus_name and 
                                    Fault['FaultType'] == (str_t+str_r2) and 
                                    Fault['Relay'].lower() == ProDevices[rr]['Name'])][0]
               
               
               M2[kk][8] = max(Fault_Data_CSV[FDI1]['Ib_mag'],Fault_Data_CSV[FDI1]['Ic_mag'])
               M2[kk][9] = max(Fault_Data_CSV[FDI2]['Ib_mag'],Fault_Data_CSV[FDI2]['Ic_mag'])
               
               M2[kk][10] = min(Fault_Data_CSV[FDI1]['Vb_mag'],Fault_Data_CSV[FDI1]['Vc_mag'])
               M2[kk][11] = min(Fault_Data_CSV[FDI2]['Vb_mag'],Fault_Data_CSV[FDI2]['Vc_mag'])
               
               M2[kk][12] = Fault_Data_CSV[FDI1]['In_mag']
               M2[kk][13] = Fault_Data_CSV[FDI2]['In_mag']

            if(len(faultBusPhases[ii])>0):
                 str_t = 'SLG_A_'
                 FDI1 = [Fault_Data_CSV.index(Fault) for Fault in Fault_Data_CSV if 
                                   (Fault['busNumber'] == fault_bus_name and 
                                    Fault['FaultType'] == (str_t+str_r1) and 
                                    Fault['Relay'].lower() == ProDevices[rr]['Name'])][0]
                 
                 FDI2 = [Fault_Data_CSV.index(Fault) for Fault in Fault_Data_CSV if 
                                   (Fault['busNumber'] == fault_bus_name and 
                                    Fault['FaultType'] == (str_t+str_r2) and 
                                    Fault['Relay'].lower() == ProDevices[rr]['Name'])][0]
                 
                 F1Iabc = [Fault_Data_CSV[FDI1]['Ia_mag'],Fault_Data_CSV[FDI1]['Ib_mag'],Fault_Data_CSV[FDI1]['Ic_mag']]
                 F2Iabc = [Fault_Data_CSV[FDI2]['Ia_mag'],Fault_Data_CSV[FDI2]['Ib_mag'],Fault_Data_CSV[FDI2]['Ic_mag']]
                 
                 M2[kk][14] = max(x for x in F1Iabc if x!=None)
                 M2[kk][15] = max(x for x in F2Iabc if x!=None)
                 
                 F1Vabc = [Fault_Data_CSV[FDI1]['Va_mag'],Fault_Data_CSV[FDI1]['Vb_mag'],Fault_Data_CSV[FDI1]['Vc_mag']]
                 F2Vabc = [Fault_Data_CSV[FDI2]['Va_mag'],Fault_Data_CSV[FDI2]['Vb_mag'],Fault_Data_CSV[FDI2]['Va_mag']]
                 M2[kk][16] = min(x for x in F1Vabc if x!=None) #Fault_Data_CSV[FDI1]['Va_mag']
                 M2[kk][17] = min(x for x in F2Vabc if x!=None) #Fault_Data_CSV[FDI2]['Va_mag']
             
                 M2[kk][18] = Fault_Data_CSV[FDI1]['In_mag']
                 M2[kk][19] = Fault_Data_CSV[FDI2]['In_mag']
                 
            kk = kk + 1
    
    del M2[kk:]
    
    for kk in range(len(M2)):
        fPoint = faultBuses[M2[kk][0]]
        rPri = ProDevices[M2[kk][1]]['Bus2']
        (isPri,BacNum) = isPrimaryRelay(rPri,fPoint,ProDevices,G)
        if(isPri):
            M2[kk][20] = 1
        elif(BacNum>=1):
            M2[kk][21] = BacNum
    
    for ii in range(len(Pairs)):
        Pairs[ii][6] = [ProDevices.index(Dev) for Dev in ProDevices if Dev['Bus1'] == Pairs[ii][0] and Dev['Bus2'] == Pairs[ii][1]][0]
        Pairs[ii][7] = [ProDevices.index(Dev) for Dev in ProDevices if Dev['Bus1'] == Pairs[ii][2] and Dev['Bus2'] == Pairs[ii][3]][0]
    
    # %% Setup Pickups
    PDN = ['Imax3ph','Imin3ph','Vmax3ph','Vmin3ph','Igmax3ph','Igmin3ph',
           'ImaxLL','IminLL','VmaxLL','VminLL','IgmaxLL','IgminLL',
           'ImaxSLG','IminSLG','VmaxSLG','VminSLG','IgmaxSLG','IgminSLG']
    for ii in range(len(ProDevices)):
        Farray = [x for x in M2 if (x[1]==ii and (x[20] == 1 or x[21]==1))]
        # if(len(Farray)==0):
        #     continue
            
        # else:
        for TT in range(0,13,6): 
            Iarr1 = [abs(F[TT+2]) for F in Farray if F[TT+2]!=None]
            Iarr2 = [abs(F[TT+3]) for F in Farray if F[TT+3]!=None]
            
            Varr1 = [abs(F[TT+4]) for F in Farray if F[TT+2]!=None]
            Varr2 = [abs(F[TT+5]) for F in Farray if F[TT+3]!=None]
            
            Narr1 = [abs(F[TT+6]) for F in Farray if F[TT+2]!=None]
            Narr2 = [abs(F[TT+7]) for F in Farray if F[TT+3]!=None]
            
            if(len(Iarr1)==0):
                MaxFault = 0
                Indmax = 0
                Iarr1 = [0]
                Varr1 = [0]
                Narr1 = [0]
                print('Error: Device'+str(ii)+' Iarr1 empty for Type '+str(TT)+'\n')
            else:            
                MaxFault = max(Iarr1)
                Indmax = Iarr1.index(MaxFault)
                
            ProDevices[ii][PDN[TT+0]] = MaxFault
            ProDevices[ii][PDN[TT+2]] = Varr1[Indmax]
            ProDevices[ii][PDN[TT+4]] = Narr1[Indmax]
    
            if(len(Iarr2)==0):
                MinFault = 0
                Indmin = 0
                Iarr2 = [0]
                Varr2 = [0]
                Narr2 = [0]
                print('Error: Device'+str(ii)+' Iarr2 empty for Type '+str(TT)+'\n')
            else:
                MinFault = min(Iarr2)
                Indmin = Iarr2.index(MinFault)
                
            ProDevices[ii][PDN[TT+1]] = MinFault
            ProDevices[ii][PDN[TT+3]] = Varr2[Indmin]
            ProDevices[ii][PDN[TT+5]] = Narr2[Indmin]
        
        # need to check wahts ok and not ok
        #MinFault_PA = [x for x in [ProDevices[ii]['Imin3ph'],ProDevices[ii]['IminLL'],ProDevices[ii]['IminSLG']] if (x!=None and x>0.1)]
        #MinFault_NA = [x for x in [ProDevices[ii]['Igmin3ph'],ProDevices[ii]['IgminLL'],ProDevices[ii]['IgminSLG']] if (x!=None and x>0.1)]
        
        MinFault_PA = [x for x in [ProDevices[ii]['Imin3ph'],ProDevices[ii]['IminLL']] if (x!=None and x>0.1)]
        MinFault_NA = [x for x in [ProDevices[ii]['IgminSLG']] if (x!=None and x>0.1)]
        #[x for x in [ProDevices[ii]['Igmin3ph'],ProDevices[ii]['IgminLL'],ProDevices[ii]['IgminSLG']] if (x!=None and x>0.1)]
        
        if(len(MinFault_PA)==0):
            MinFault_SLG = [x for x in [ProDevices[ii]['IminSLG']] if (x!=None and x>0.1)]
            if(len(MinFault_SLG)==0):
                MinFault_P = ProDevices[ii]['IL']*2
            else:
                MinFault_P = min(MinFault_SLG)
        else:
            MinFault_P = min(MinFault_PA)
            
        if(len(MinFault_NA)==0):
            MinFault_N = ProDevices[ii]['In']*2
        else:
            MinFault_N = min(MinFault_NA)
              
        if(ProDevices[ii]['Type'] == 'Relay_TOC' or ProDevices[ii]['Type'] == 'Relay_DT' or ProDevices[ii]['Type'] == 'Rec'):
            if(ProDevices[ii]['IL'] > 0):
                # Phase Pickup
                if(MinFault_P>ProDevices[ii]['IL']*1.5):
                    ProDevices[ii]['Ip'] = ProDevices[ii]['IL']*1.25
                else:
                    ProDevices[ii]['Ip'] = MinFault_P*0.5
                # Ground Pickup
                if(MinFault_N>ProDevices[ii]['In']*1.5):
                    ProDevices[ii]['Inp'] = ProDevices[ii]['In']*1.25
                else:
                    ProDevices[ii]['Inp'] = MinFault_N*0.5
            else:
                # Phase Pickup
                ProDevices[ii]['Ip'] = MinFault_P*0.1
                ProDevices[ii]['Inp'] = MinFault_N*0.1
                
            ProDevices[ii]['Ip'] = round(ProDevices[ii]['Ip'],1)
            ProDevices[ii]['Inp'] = round(ProDevices[ii]['Inp'],1)    
                
            # check Pickup Needs To be updated to make variable
            if(ProDevices[ii]['Ip']<=Min_Ip[0]):
                ProDevices[ii]['Ip'] = Min_Ip[0]
            if(ProDevices[ii]['Inp']<=Min_Ip[1]):
                ProDevices[ii]['Inp'] = Min_Ip[1]
            
            if(enableIT==1):
                MaxFault_Pri = [x for xs in [[F[2],F[3],F[8],F[9],F[14],F[15]] for F in Farray if F[20]==1] for x in xs if x != None]
                MaxFault_Bac = [x for xs in [[F[2],F[3],F[8],F[9],F[14],F[15]] for F in Farray if F[21]==1] for x in xs if x != None]
                
                MaxFault0_Pri = [x for xs in [[F[6],F[7],F[12],F[13],F[18],F[19]] for F in Farray if F[20]==1] for x in xs if x != None]
                MaxFault0_Bac = [x for xs in [[F[6],F[7],F[12],F[13],F[18],F[19]] for F in Farray if F[20]==1] for x in xs if x != None]
                
                MaxPri = max(MaxFault_Pri)
                MaxBac = max(MaxFault_Bac)
                MaxPri0= max(MaxFault0_Pri)
                MaxBac0= max(MaxFault0_Bac)
                if(len(MaxBac)!=0):
                    if(MaxPri/MaxBac >= 1.3):
                        ProDevices[ii]['IT'] = MaxBac * 1.25
                    else:
                        ProDevices[ii]['IT'] = 0
                else:
                    ProDevices[ii]['IT'] = 0
                    
                if(len(MaxBac0)!=0):
                    if(MaxPri0/MaxBac0 >= 1.3):
                        ProDevices[ii]['ITg'] = MaxBac * 1.25
                    else:
                        ProDevices[ii]['ITg'] = 0
                else:
                    ProDevices[ii]['ITg'] = 0
            else:
                ProDevices[ii]['IT'] = 0
                ProDevices[ii]['ITg'] = 0  
        elif(ProDevices[ii]['Type'] == 'Fuse'):
            ProDevices[ii]['Ip'] = 'T10'

    # %% Find Relevent Devices 
    Relevent_Devs = list(set([P[6] for P in Pairs]+[P[7] for P in Pairs]))
    nOCDevs = 0
    nCtrlVars = 0
    for dev in Relevent_Devs:
        if('Relay' in ProDevices[dev]['Type']):
            nOCDevs+=1
            nCtrlVars += 4        
        elif('Rec' in ProDevices[dev]['Type']):
            nOCDevs+=2 
            nCtrlVars += 8  
        else:
            pass
        
    ProDevRef = [[] for i in range(0,nOCDevs)]
    kk=0
    for ii in Relevent_Devs:
        if('Relay' in ProDevices[ii]['Type']):
            ProDevRef[kk] = deepcopy(ProDevices[ii])
            ProDevRef[kk]['Oind'] = ii
            kk = kk + 1
        elif('Rec' in ProDevices[ii]['Type']):
            ProDevRef[kk] = deepcopy(ProDevices[ii])
            ProDevRef[kk]['Type'] = ProDevRef[kk]['Type'] +'_'+'Fast'
            ProDevRef[kk]['Oind'] = ii
            kk = kk + 1
            
            ProDevRef[kk] = deepcopy(ProDevices[ii])
            ProDevRef[kk]['Type'] = ProDevRef[kk]['Type'] +'_'+'Slow'
            ProDevRef[kk]['Oind'] = ii
            kk = kk + 1
        else:
            ProDevRef[kk] = deepcopy(ProDevices[ii])
            ProDevRef[kk]['Oind'] = ii
            kk = kk + 1
 
    # %% Pick device type
    if(type_select):
        dt_th = 5    
        for ii in range(len(ProDevRef)):
            
            if('Fuse' in [x[4] for x in Pairs if x[7] == ProDevRef[ii]['Oind']]):
                ProDevRef[ii]['Type'] = 'Relay_TOC'
            else:
                Mp_min = min(ProDevRef[ii]['Imin3ph'],ProDevRef[ii]['IminLL'])/(ProDevRef[ii]['Ip']) 
                Mp_max = max(ProDevRef[ii]['Imax3ph'],ProDevRef[ii]['ImaxLL'],ProDevRef[ii]['ImaxSLG'])/(ProDevRef[ii]['Ip']) 
                Mn_min = ProDevRef[ii]['IgminSLG']/(ProDevRef[ii]['In']*1.25) 
                Mn_max = ProDevRef[ii]['IgmaxSLG']/(ProDevRef[ii]['In']*1.25) 
                ProDevRef[ii]['MX'] = [Mp_min,Mp_max,Mn_min,Mn_max]
                
                if(Mp_min != 0 and min([Mp_min,Mn_min])<=dt_th and Mp_max/Mp_min <=100):
                    ProDevRef[ii]['Type'] = 'Relay_DT' 
                # elif(min([Mp_min,Mn_min])>dt_th and Mp_max/Mp_min >100):
                #     ProDevRef[ii]['Type'] = 'Relay_Dist'  
                else:
                    ProDevRef[ii]['Type'] = 'Relay_TOC'
    
    # %% Check if Optimization is needed/ossible
    if(len(Pairs)==0 or len(ProDevRef)==0):
        print('No coordiantion pairs found')
        Relay_settings= [dict.fromkeys(['Name','From','To','PickupI','TDS','TOC','PickupI0','TDSg','TOCg','VR','IT','IT0']) for number in range(len(ProDevices))]
        for ii in range(len(ProDevices)):
            Relay_settings[ii]['Name'] = ProDevices[ii]['Name']+'_'+ProDevices[ii]['Direction']
            Relay_settings[ii]['From'] = ProDevices[ii]['Bus1']
            Relay_settings[ii]['To'] = ProDevices[ii]['Bus2']
            Relay_settings[ii]['PickupI'] = ProDevices[ii]['Ip']
            Relay_settings[ii]['TDS'] = 1
            Relay_settings[ii]['TOC'] = OCTCC_Name(1)
            Relay_settings[ii]['PickupI0'] = ProDevices[ii]['Inp']
            Relay_settings[ii]['TDSg'] = 1
            Relay_settings[ii]['TOCg'] = OCTCC_Name(1)
            Relay_settings[ii]['VR'] = False
            Relay_settings[ii]['IT'] = ProDevices[ii]['IT']
            Relay_settings[ii]['IT0'] = ProDevices[ii]['ITg']
        info = {'SWs':switchStates,
                'DOC':DOC,'Force_NOIBR':Force_NOIBR,
                'rerun':-1,'CTImin':0,
                'Pro':ProDevices,
                'CTI':CTI,
                'fobj':0,
                'Load':Device_Data_CSV,
                'G':G,
                'sol': [],
                'M2':M2,
                'FBuses':faultBuses}
        return Relay_settings,info
    else:
        print('Calculating Settings')
    # %% Create control variables and limits
    CVMin = [None]*nCtrlVars
    CVMax = [None]*nCtrlVars
    CtrlDev  = [0]*len(ProDevRef)
    kk=0
    for ii in range(len(ProDevRef)):
        if('Relay_TOC' == ProDevRef[ii]['Type']):
            
            CVMin[kk] = 5
            CVMax[kk] = 150
            CVMin[kk+1] = 5
            CVMax[kk+1] = 150
            if(ProDevRef[ii]['Ip'] > ProDevRef[ii]['IL']):
                CVMin[kk+2] = 1
                CVMax[kk+2] = 5
            else:
                CVMin[kk+2] = 7
                CVMax[kk+2] = 11
            
            if(ProDevRef[ii]['Inp'] > ProDevRef[ii]['In']):
                CVMin[kk+3] = 1
                CVMax[kk+3] = 5
            else:
                CVMin[kk+3] = 7
                CVMax[kk+3] = 11   
            kk=kk+4
            CtrlDev[ii] = kk
        
        elif('Relay_DT' == ProDevRef[ii]['Type']):
            CVMin[kk] = 0.02
            CVMax[kk] = 60        
            CVMin[kk+1] = 0.02
            CVMax[kk+1] = 60
            if(ProDevRef[ii]['Ip'] >ProDevRef[ii]['IL']):
                CVMin[kk+2] = 6
                CVMax[kk+2] = 6
            else:
                CVMin[kk+2] = 12
                CVMax[kk+2] = 12
            
            if(ProDevRef[ii]['Inp'] > ProDevRef[ii]['In']):
                CVMin[kk+3] = 6
                CVMax[kk+3] = 6
            else:
                CVMin[kk+3] = 12
                CVMax[kk+3] = 12   
            kk=kk+4
            CtrlDev[ii] = kk
        elif('Rec'in ProDevRef[ii]['Type']):
            CVMin[kk] = 5
            CVMax[kk] = 150        
            CVMin[kk+1] = 5
            CVMax[kk+1] = 150
            if(ProDevRef[ii]['Ip'] > ProDevRef[ii]['IL']):
                CVMin[kk+2] = 1
                CVMax[kk+2] = 5
            else:
                CVMin[kk+2] = 7
                CVMax[kk+2] = 11
            
            if(ProDevRef[ii]['Inp'] > ProDevRef[ii]['In']):
                CVMin[kk+3] = 1
                CVMax[kk+3] = 5
            else:
                CVMin[kk+3] = 7
                CVMax[kk+3] = 11   
            kk=kk+4
            CtrlDev[ii] = kk
        elif('Dist' in ProDevRef[ii]['Type']):
            # pretend TOC for testing 
            CVMin[kk] = 5
            CVMax[kk] = 150        
            CVMin[kk+1] = 5
            CVMax[kk+1] = 150
            if(ProDevRef[ii]['Ip'] > ProDevRef[ii]['IL']):
                CVMin[kk+2] = 1
                CVMax[kk+2] = 5
            else:
                CVMin[kk+2] = 7
                CVMax[kk+2] = 11
            
            if(ProDevRef[ii]['Inp'] > ProDevRef[ii]['In']):
                CVMin[kk+3] = 1
                CVMax[kk+3] = 5
            else:
                CVMin[kk+3] = 7
                CVMax[kk+3] = 11   
            kk=kk+4
            CtrlDev[ii] = kk
            
        else:
            print('Error unknown Relay type')
            
            
    # %% Write Objective  
    fundir = Main_dir
    print('Fun Dir = '+fundir+'\n')
    print('Writing Obj File')
    
    Objfile = 'F1'
    Obj_fileID = write_GA_fit_fun(CtrlDev,ProDevRef,Pairs,M2,fundir,Objfile,CTI,OTmax)
    con_fileID = write_Con(CtrlDev,ProDevRef,Pairs,M2,fundir,Objfile,CTI)
    
    # %% Run GA
    pop_size = 50
    cti_min = 0
    rerun = 0
    fobj = -100000

    if(GA_initial_seed!=None and len(GA_initial_seed)==len(CVMin)):
        initpop = [[0]*len(CVMax) for x in range(pop_size)]
        for ii in range(pop_size-10):
            for jj in range(len(CVMax)):
                if(type(CVMin[jj]) == int):
                    initpop[ii][jj] = random.randint(CVMin[jj], CVMax[jj])
                else:
                    initpop[ii][jj] = random.random()*CVMax[jj]
        for ii in range(50-10,50):
            initpop[ii] = GA_initial_seed
    else:
        initpop = None
    
    # Fake int pop
    
    
    
    while(rerun<2 and (cti_min<0.24 or fobj<-1000)):
        if(rerun>=2):
            initpop = None
            
        last_fitness = 0 
        import importlib
        import sys
        # Reload Fit and con functions 
        spec = importlib.util.spec_from_file_location(Objfile, os.path.join(Main_dir,f'{Objfile}.py'))
        module = importlib.util.module_from_spec(spec)
        sys.modules[Objfile] = module
        spec.loader.exec_module(module)
        fitness_func = module.fitness_func
        
        spec1 = importlib.util.spec_from_file_location(Objfile+'_con', os.path.join(Main_dir,f'{Objfile}_con.py'))
        module1 = importlib.util.module_from_spec(spec1)
        sys.modules[Objfile+'_con'] = module1
        spec1.loader.exec_module(module1)
        Con_func = module1.Con_func
        print('Running optimizer')
        gs = [{"low": CVMin[x], "high": CVMax[x]} for x in range(len(CVMin))]
        ga_instance = pygad.GA(num_generations=10000,
                                num_parents_mating=5,
                                sol_per_pop=50,
                                num_genes=len(CVMax),
                                gene_space=gs,
                                mutation_by_replacement=True,
                                fitness_func= lambda ga_instance, solution, solution_idx : fitness_func(solution,solution_idx),
                                on_generation=on_generation,
                                gene_type=int,
                                save_solutions = False,
                                #initial_population = initpop,
                                stop_criteria = ["saturate_50"])
        
        ga_instance.run()
        
        if(Sho_Plots):
            ga_instance.plot_fitness()
        
        solution, solution_fitness, solution_idx = ga_instance.best_solution(ga_instance.last_generation_fitness)
        #print("Solution", solution)
        print("Fitness value of the best solution = {solution_fitness}".format(solution_fitness=solution_fitness))
        
        Pen,gac,Ttot = Con_func(solution,1)
        fobj = solution_fitness;
        print("OT : Avg = {OTavg},\t min = {OTmin},\t Max = {OTmax}".format(OTavg=sum(Ttot)/len(Ttot),
                                                                            OTmin=min(Ttot),
                                                                            OTmax=max(Ttot)))
        print("CTI: Avg = {OTavg},\t min = {OTmin},\t Max = {OTmax}".format(OTavg=(sum(gac)/len(gac))+CTI,
                                                                            OTmin=min(gac)+CTI,
                                                                            OTmax=max(gac)+CTI))   
        cti_min = min(gac)+CTI
        rerun = rerun+1
    # %% Process Results
    
    Relay_settings= [dict.fromkeys(['Name','From','To','PickupI','TDS','TOC','PickupI0','TDSg','TOCg','VR','IT','IT0']) for number in range(len(ProDevRef))]
    
    for ii in range(len(ProDevRef)):
        if('Relay_TOC' == ProDevRef[ii]['Type'] or 'Relay_DT' == ProDevRef[ii]['Type']):
            Relay_settings[ii]['Name'] = ProDevRef[ii]['Name']+'_'+ProDevRef[ii]['Direction']
            Relay_settings[ii]['From'] = ProDevRef[ii]['Bus1']
            Relay_settings[ii]['To'] = ProDevRef[ii]['Bus2']
            Relay_settings[ii]['PickupI'] = ProDevRef[ii]['Ip']
            Relay_settings[ii]['TDS'] = solution[CtrlDev[ii]-4]/10
            Relay_settings[ii]['TOC'] = OCTCC_Name(solution[CtrlDev[ii]-2])
            Relay_settings[ii]['PickupI0'] = ProDevRef[ii]['Inp']
            Relay_settings[ii]['TDSg'] = solution[CtrlDev[ii]-3]/10
            Relay_settings[ii]['TOCg'] = OCTCC_Name(solution[CtrlDev[ii]-1])
            Relay_settings[ii]['VR'] = [solution[CtrlDev[ii]-2]>6 or solution[CtrlDev[ii]-1]>6][0]
            Relay_settings[ii]['IT'] = ProDevRef[ii]['IT']
            Relay_settings[ii]['IT0'] = ProDevRef[ii]['ITg']
            if(DOC == 1 and SetDir):
                Relay_settings[ii]['Z1MAG'] = ProDevRef[ii]['Z1MAG']
                Relay_settings[ii]['Z1ANG'] =ProDevRef[ii]['Z1ANG']
                Relay_settings[ii]['Z0MAG'] = ProDevRef[ii]['Z0MAG']
                Relay_settings[ii]['Z0ANG'] =ProDevRef[ii]['Z0ANG']
                Relay_settings[ii]['Z2F'] =ProDevRef[ii]['Z2F']
                Relay_settings[ii]['Z2R'] =ProDevRef[ii]['Z2R']
                
    info = {'SWs':switchStates,'DOC':DOC,'Force_NOIBR':Force_NOIBR,'rerun':rerun,'CTImin':cti_min,'Pro':ProDevRef,'CTI':CTI,'fobj':fobj,'Load':Device_Data_CSV,'G':G,'sol': solution,'M2':M2,'FBuses':faultBuses}
    return Relay_settings,info


