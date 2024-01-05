# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 01:58:04 2023

@author: maste
"""

import os
import json
import RSO_pack

# %% 
pwd = os.getcwd()
Main_dir = pwd
# import data for RONM JSON file
jsonFile = os.path.join(pwd,'Examples','RONM_Data','output.ieee13.faults.json') 
Main_dir = pwd
f = open(jsonFile) 
jsonDict = json.load(f)
f.close()

# %% Get SysInfo from Json 
proSettings = jsonDict['Protection settings']
powerFlow = jsonDict['Powerflow output']

SysInfo = RSO_pack.getRONMSysInfo(proSettings,powerFlow,ignore_fuse=True)
# Get Timeline Data form Json
devTimeLine = jsonDict['Device action timeline']
(sW_Status,sW_Names) = RSO_pack.get_Sw_Status(jsonDict['Powerflow output'],devTimeLine)

# %% Calculate devie Data 
Buses = SysInfo['Buses']
Relay_list = [x['MonitoredObj'].split('line.')[1] for x in SysInfo['Relays']] 
Reclo_list = [x['MonitoredObj'].split('line.')[1] for x in SysInfo['Recs']]
Fuses_list = [x['MonitoredObj'].split('line.')[1] for x in SysInfo['Fuses']]
devTypes = ['relay']*len(Relay_list)+['recloser']*len(Reclo_list)+['fuse']*len(Fuses_list)
devLines = Relay_list+Reclo_list+Fuses_list
devNames = [x['Name'] for x in SysInfo['Relays']]+[x['Name'] for x in SysInfo['Recs']]+[x['Name'] for x in SysInfo['Fuses']]
dev_BusV = [Buses[RSO_pack.index_dict(Buses,'Name',x['Bus1'])]['kV']*1e3 for x in SysInfo['Relays']]+[Buses[RSO_pack.index_dict(Buses,'Name',x['Bus1'])]['kV']*1e3 for x in SysInfo['Recs']]+[Buses[RSO_pack.index_dict(Buses,'Name',x['Bus1'])]['kV']*1e3 for x in SysInfo['Fuses']]

# create settings list 
settings = []*len(devTimeLine)

# %% itterate per step and calcualte settings  
for Ts in range(len(devTimeLine)):
    print('------======'+str(Ts)+'======------\n')
    # get switch states
    switchLines = sW_Names
    switchStates = sW_Status[Ts]
    # Collect Load Flow Deta 
    Device_Data_CSV = RSO_pack.getRONMDeviceData(Ts,powerFlow,devTypes,devNames,devLines,dev_BusV,SysInfo,sW_Status,sW_Names)

    # collect Fault Data
    faultBuses = list(jsonDict['Fault currents'][1].keys())
    faultBusesLOC = faultBuses
    faultBusPhases = [None]*len(faultBuses)
    for ii in range(len(faultBuses)):
        faultBusPhases[ii] = Buses[RSO_pack.index_dict(Buses,'Name',faultBuses[ii])]['nodes']
    
    Fault_Data_CSV = RSO_pack.getRONMFaultData(Ts, faultBuses, devTypes, devNames, devLines, dev_BusV, jsonDict['Fault currents'][Ts])
    
    # %% Run optimizer
    
    # settigns 
    Force_NOIBR = 0
    DOC = 1
    enableIT = 0
    CTI = 0.25
    OTmax = 10
    Sho_Plots=False
    type_select = False
    Min_Ip = [1,1]
    Fault_Res = list(set(['R'+x['FaultType'].split('_R')[1] for x in Fault_Data_CSV]))
    Substation_bus = 'sourcebus'
    initpop = None
    SetDir=False
    
    settings,old_info = RSO_pack.runSettingsOptimizer(Main_dir,switchStates,switchLines,Device_Data_CSV,Fault_Data_CSV,Fault_Res,SysInfo,Substation_bus,Min_Ip,enableIT,Force_NOIBR,CTI,OTmax,type_select,SetDir,Sho_Plots,GA_initial_seed=initpop)



