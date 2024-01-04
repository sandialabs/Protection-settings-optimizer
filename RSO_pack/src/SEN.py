# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 11:08:44 2022

@author: maste
"""
import os
import pygad
import networkx as nx
from copy import deepcopy
from .on_gen import on_generation
from .OT_funs import OCTCC_Name,OCTCC_Num
from .GenNxGraph import genGraph, index_dict
from .write_fit_fun import write_GA_fit_fun, write_Con 
from .PathTracingFuns import faultPathDir, isPrimaryRelay
from .Read_CSV_Functions import read_Relays_CSV_Data, read_Fault_CSV_Data
from .PathTracingFuns import find_faultpath_insys, findDeviceInPath, create_priamry_backup_from_paths

import importlib
import sys

def checkSensivity(Main_dir,IEEE123_SysInfo,old_info,settings,Force_NOIBR,Device_Data_File,Fault_Data_File,Substation_bus,enableIT,OT_max_allow):
    # %% get curretn config
    os.chdir(Main_dir)
    Device_Data_CSV = read_Relays_CSV_Data(Device_Data_File)
    switchStates = [1,1,1,1]
    switchLines = ['sw2','sw3','sw7','sw4']
    #               RTL1 RTL2 RTL3 RTL4
    
    Run_Flag = False
    Error_code = 0
    Error_Value = 0
    
    for ii in range(len(Device_Data_CSV)):
        if(Device_Data_CSV[ii]['switchState'] == False):
            if(Device_Data_CSV[ii]['RelayName'] == 'RTL1'):
                switchStates[0] = 0
            elif(Device_Data_CSV[ii]['RelayName'] == 'RTL2'):
                switchStates[1] = 0
            elif(Device_Data_CSV[ii]['RelayName'] == 'RTL3'):
                switchStates[2] = 0
            elif(Device_Data_CSV[ii]['RelayName'] == 'RTL4'):
                switchStates[3] = 0
            else:
                pass
    
    # %% check if switch state has change
    if(switchStates != old_info['SWs']):
        Run_Flag = True;
        Error_code = 'SEN: Config changed from' + str(old_info['SWs']) + ' to ' + str(switchStates)
        Error_Value = 0
        return Run_Flag,Error_code,Error_Value
    
    

    # %% get old PV data and check if similar
    if(Force_NOIBR == 1):
        for ii in range(len(Device_Data_CSV)):
            Device_Data_CSV[ii]['PVs'] = '0'

    if(old_info['Force_NOIBR'] == Force_NOIBR):
        pass
    else:
        old_pv = [x['PVs'] for x in old_info['Load']]
        new_pv = [y['PVs'] for y in Device_Data_CSV]
        if(len(old_pv) != len(new_pv)):
            Run_Flag = True;
            Error_code = 'SEN: change in PV status'
            Error_Value = -1
            return Run_Flag,Error_code,Error_Value
        else:
            for ii in range(len(old_pv)):
                if(old_pv[ii] != new_pv[ii]):
                    Run_Flag = True;
                    Error_code = 'SEN: change in PV status'
                    Error_Value = ii
                    return Run_Flag,Error_code,Error_Value

    # %% get old load data nd check if similar
    old_IL = [x['Ia_mag'] for x in old_info['Load']]
    new_IL = [y['Ia_mag'] for y in Device_Data_CSV]
    
    old_LIn = [x['Ia_mag'] for x in old_info['Load']]
    new_LIn= [y['Ia_mag'] for y in Device_Data_CSV]
    
    
    for ii in range(len(old_IL)):
        if( (abs(old_IL[ii]-new_IL[ii])/old_IL[ii]) > 0.15):
            Run_Flag = True;
            Error_code = 'SEN: Load current change Imax'
            Error_Value = (abs(old_IL[ii]-new_IL[ii])/old_LIn[ii])
            return Run_Flag,Error_code,Error_Value
        
        elif( (abs(old_LIn[ii]-new_LIn[ii])/old_LIn[ii]) > 0.15 ):
            Run_Flag = True;
            Error_code = 'SEN: Load curretn Change In'
            Error_Value = (abs(old_LIn[ii]-new_LIn[ii])/old_LIn[ii])
            return Run_Flag,Error_code,Error_Value

    # %% Read Fault input file
    Fault_Data_CSV = read_Fault_CSV_Data(Fault_Data_File)
    
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
            
    # %% get previous settigns and test
    old_settings = settings

    # Read system data
    Relays = IEEE123_SysInfo['Relays']
    Recs = IEEE123_SysInfo['Recs']
    Fuses = IEEE123_SysInfo['Fuses']
    Lines = IEEE123_SysInfo['Lines']
    XFMRs = IEEE123_SysInfo['XFMRs']
    Buses = IEEE123_SysInfo['Buses']
    Pvs = IEEE123_SysInfo['Pvs']
    BESS = IEEE123_SysInfo['BESS']
    Gens = IEEE123_SysInfo['Gens']

    # %%Edit system based on settings
    for ii in range(len(Pvs)):
        if(Device_Data_CSV[ii]['PVs'] == '0'):
            Pvs[ii]['Enabled'] = False    
    
        
    if(sum([x['Enabled'] for x in Pvs] + [x['Enabled'] for x in BESS] + [x['Enabled'] for x in Gens])<=0):
        DOC = 0
    else:
        DOC = 1   
        
    for ii in range(len(switchLines)):
        if(switchStates[ii]==1):
            #dssText.Command = 'edit line.'+switchLines[ii]+' enabled=true'
            Lines[index_dict(Lines,'Name',switchLines[ii])]['Enabled'] = True
        elif(switchStates[ii]==0):
            #dssText.Command = 'edit line.'+switchLines[ii]+' enabled=false'
            Lines[index_dict(Lines,'Name',switchLines[ii])]['Enabled'] = False    



    # %% Generate Net Graph        
    Gr = genGraph(Lines,XFMRs,Buses,Relays,Recs,Fuses,Pvs,BESS,Gens,switchLines,switchStates)
    G = Gr[0]
    pos = Gr[1]
    Edges = Gr[2]
    
    Node_list = list(G.nodes)
    Edge_list = list(G.edges)
    if(Substation_bus in G.nodes):
        G.nodes[Substation_bus]['isSource'] = True
        G.nodes[Substation_bus]['sourceName'] = 'Sub'

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
    ProDevices = [dict.fromkeys(['Bus1','Bus2','phases','Line','Direction','Name','Type','Vpu','IL','In','Ip','Inp','Imax3ph','Imin3ph','Vmax3ph','Vmin3ph','Igmax3ph','Igmin3ph',\
                  'ImaxLL','IminLL','VmaxLL','VminLL','IgmaxLL','IgminLL','ImaxSLG','IminSLG','VmaxSLG','VminSLG','IgmaxSLG','IgminSLG','IT','ITg','Oind']) for number in range(nProDevices)]
    
    kk = 0;
    for ii in range(len(Relays)):
        # Relays
        if(Relays[ii]['MonitoredObj'].split('.')[1].lower() in [Name[2] for Name in list(G.edges.data('Name'))]):
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
        if(Recs[ii]['MonitoredObj'].split('line.')[1].lower() in [Name[2] for Name in list(G.edges.data('Name'))]):
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
        if(Fuses[ii]['MonitoredObj'].lower().split('line.')[1] in [Name[2] for Name in list(G.edges.data('Name'))]):
            ProDevices[kk]['Bus1'] = Fuses[ii]['Bus1']
            ProDevices[kk]['Bus2'] = Fuses[ii]['Bus2']
            ProDevices[kk]['Line'] = [Edges[x]['Name'] for x in range(len(Edges))].index(Fuses[ii]['MonitoredObj'].lower().split('.')[1])
            ProDevices[kk]['phases'] = Edges[ProDevices[kk]['Line']]['numPhases']
            ProDevices[kk]['Direction'] = 'E' #add logic later
            ProDevices[kk]['Name'] = Fuses[ii]['Name']
            ProDevices[kk]['Type'] = Fuses[ii]['Type'] #add Logic later
            kk=kk+1

    # %% Generate or Get Fault bus Locations
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

    # %% Build M2
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
            
            str_r1 = 'R0_05'
            str_r2 = 'R1'
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
                 
                 M2[kk][14] = Fault_Data_CSV[FDI1]['Ia_mag']
                 M2[kk][15] = Fault_Data_CSV[FDI2]['Ia_mag']
                 
                 M2[kk][16] = Fault_Data_CSV[FDI1]['Va_mag']
                 M2[kk][17] = Fault_Data_CSV[FDI2]['Va_mag']
             
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
                print('Error: Device'+str(ii)+' Iarr1 empty \n')
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
                print('Error: Device'+str(ii)+' Iarr2 empty \n')
            else:
                MinFault = min(Iarr2)
                Indmin = Iarr2.index(MinFault)
                
            ProDevices[ii][PDN[TT+1]] = MinFault
            ProDevices[ii][PDN[TT+3]] = Varr2[Indmin]
            ProDevices[ii][PDN[TT+5]] = Narr2[Indmin]
        
        MinFault_PA = [x for x in [ProDevices[ii]['Imin3ph'],ProDevices[ii]['IminLL'],ProDevices[ii]['IminSLG']] if (x!=None and x>0.1)]
        MinFault_NA = [x for x in [ProDevices[ii]['Igmin3ph'],ProDevices[ii]['IgminLL'],ProDevices[ii]['IgminSLG']] if (x!=None and x>0.1)]
        if(len(MinFault_PA)==0):
            MinFault_P = 0
        else:
            MinFault_P = min(MinFault_PA)
            
        if(len(MinFault_NA)==0):
            MinFault_N = 0
        else:
            MinFault_N = min(MinFault_NA)
              
        if(ProDevices[ii]['Type'] == 'Relay_TOC' or ProDevices[ii]['Type'] == 'Relay_DT' or ProDevices[ii]['Type'] == 'Rec'):
            if(ProDevices[ii]['IL'] > 0):
                # Phase Pickup
                if(MinFault_P>ProDevices[ii]['IL']*1.5):
                    ProDevices[ii]['Ip'] = ProDevices[ii]['IL']*1.25
                else:
                    ProDevices[ii]['Ip'] = MinFault_P*0.9
                # Ground Pickup
                if(MinFault_N>ProDevices[ii]['In']*1.5):
                    ProDevices[ii]['Inp'] = ProDevices[ii]['In']*1.25
                else:
                    ProDevices[ii]['Inp'] = MinFault_N*0.9
            else:
                # Phase Pickup
                ProDevices[ii]['Ip'] = MinFault_P*0.1
                ProDevices[ii]['Inp'] = MinFault_N*0.1
                
            ProDevices[ii]['Ip'] = round(ProDevices[ii]['Ip'],1)
            ProDevices[ii]['Inp'] = round(ProDevices[ii]['Inp'],1)    
                
            # check Pickup
            if(ProDevices[ii]['Ip']<=60):
                ProDevices[ii]['Ip'] = 60
            if(ProDevices[ii]['Inp']<=5):
                ProDevices[ii]['Inp'] = 5
            
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

    # %% Create control variables and limits
    CtrlDev  = [0]*len(ProDevRef)
    kk=0
    
    
    for ii in range(len(ProDevRef)):
        
        if('Relay_TOC' == ProDevRef[ii]['Type']):
            kk=kk+4
            CtrlDev[ii] = kk
        elif('Relay_DT' == ProDevRef[ii]['Type']):
            kk=kk+4
            CtrlDev[ii] = kk
        elif('Rec'in ProDevRef[ii]['Type']):
            kk=kk+4
            CtrlDev[ii] = kk
        elif('Dist' in ProDevRef[ii]['Type']):
            kk=kk+4
            CtrlDev[ii] = kk
            
        else:
            print('Error unknown Relay type')

    XX = [None] * CtrlDev[-1]

    for ii in range(len(ProDevRef)):
        Set_ind = [settings.index(x) for x in settings if x['From'] == ProDevRef[ii]['Bus1'] and x['To'] == ProDevRef[ii]['Bus2'] ][0]
        #print(Set_ind)
        ProDevRef[ii]['Ip'] = settings[Set_ind]['PickupI']
        ProDevRef[ii]['Inp'] = settings[Set_ind]['PickupI0']
        XX[CtrlDev[ii]-4] = settings[Set_ind]['TDS'] * 10
        XX[CtrlDev[ii]-3] = settings[Set_ind]['TDSg'] * 10
        XX[CtrlDev[ii]-2] = OCTCC_Num(settings[Set_ind]['TOC'],settings[Set_ind]['VR'])
        XX[CtrlDev[ii]-1] = OCTCC_Num(settings[Set_ind]['TOCg'],settings[Set_ind]['VR'])
    #print(XX)
        
    Objfile = 'SEN1'
    CTI = old_info['CTI']
    fundir = Main_dir
    write_Con(CtrlDev,ProDevRef,Pairs,M2,fundir,Objfile,CTI)
    
    spec1 = importlib.util.spec_from_file_location(Objfile+'_con', os.path.join(Main_dir,f'{Objfile}_con.py'))
    module1 = importlib.util.module_from_spec(spec1)
    sys.modules[Objfile+'_con'] = module1
    spec1.loader.exec_module(module1)
    Con_func = module1.Con_func
    
    Pen,gac,Ttot = Con_func(XX,1)

    print("OT : Avg = {OTavg},\t min = {OTmin},\t Max = {OTmax}".format(OTavg=sum(Ttot)/len(Ttot),
                                                                        OTmin=min(Ttot),
                                                                        OTmax=max(Ttot)))
    print("CTI: Avg = {OTavg},\t min = {OTmin},\t Max = {OTmax}".format(OTavg=(sum(gac)/len(gac))+CTI,
                                                                        OTmin=min(gac)+CTI,
                                                                        OTmax=max(gac)+CTI))   
    cti_min = min(gac)+CTI
    if(cti_min<0.98*CTI):
        Run_Flag = True
        Error_code = 'SEN: minimum CTI Volated'
        Error_Value = cti_min 
        return Run_Flag,Error_code,Error_Value
    elif(max(Ttot)>OT_max_allow):
        Run_Flag = True
        Error_code = 'SEN: Maximum operating time > OT_Max'
        Error_Value = max(Ttot) 
        return Run_Flag,Error_code,Error_Value
    
    # none match return
    return Run_Flag,Error_code,Error_Value

