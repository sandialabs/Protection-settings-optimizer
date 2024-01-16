# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 22:10:50 2023

@author: maste
"""
# Example OpenDSS Input 
import os
import RSO_pack
import opendssdirect as dss

pwd = os.getcwd()

# Start Run IEEE 34 bus system
dss.run_command('clear')
dss.run_command('compile '+ os.path.join(pwd,'Examples','IEEE34_OpenDSS','IEEE34Test.dss'))
dss.run_command('set maxcontroliter = 500')
dss.run_command('solve')

dssText = dss.Text
dssCircuit = dss.Circuit

# collect system info from OpenDSS
SysInfo = RSO_pack.getSysInfo(dssCircuit)

# %% collect steady state Data
Buses = SysInfo['Buses']
devLines = [x['MonitoredObj'].split('Line.')[1] for x in SysInfo['Relays']]+[x['MonitoredObj'].split('Line.')[1] for x in SysInfo['Recs']]
devNames = [x['Name'] for x in SysInfo['Relays']]+[x['Name'] for x in SysInfo['Recs']]
dev_BusV = [Buses[RSO_pack.index_dict(Buses,'Name',x['Bus1'])]['kV']*1e3 for x in SysInfo['Relays'] ]+[Buses[RSO_pack.index_dict(Buses,'Name',x['Bus1'])]['kV']*1e3 for x in SysInfo['Recs'] ]

Device_Data_CSV = RSO_pack.getDeviceData(dssCircuit,devNames,devLines,dev_BusV)

# %% collect fault Data 
faultBuses = [x['Name'] for x in Buses]
faultBusPhases = [None]*len(faultBuses)
for ii in range(len(faultBuses)):
    faultBusPhases[ii] = Buses[RSO_pack.index_dict(Buses,'Name',faultBuses[ii])]['nodes']

Fres = ['0.001','1']
Fts = ['3ph','SLG','LL']

FData = RSO_pack.getFaultInfo(dssCircuit,dssText,faultBuses,faultBusPhases,Fres,Fts,devLines,devNames,dev_BusV)
Fault_File_loc = os.path.join(pwd,'FData.csv')
FData.to_csv(Fault_File_loc,index=False,header=False)
Fault_Data_CSV = RSO_pack.read_Fault_CSV_Data(Fault_File_loc)


# %% Program settigns 
Force_NOIBR = 1
DOC = 0
enableIT = 0
CTI = 0.25
OTmax = 10
Sho_Plots=1
type_select = False
Fault_Res = ['R0_001','R1']
Min_Ip = [0.1,0.1]
Substation_bus = 'sourcebus'
initpop = None
SetDir=False

Main_dir = pwd

switchLines = ['Sw1','Sw2a','Sw2b','Sw3','Sw4','Sw5','Sw6']
switchStates = [1  ,1   ,1    ,1   ,1    ,1    ,1    ]
settings,old_info = RSO_pack.runSettingsOptimizer(Main_dir,
                                                  switchStates,
                                                  switchLines,
                                                  Device_Data_CSV,
                                                  Fault_Data_CSV,
                                                  Fault_Res,
                                                  SysInfo,
                                                  Substation_bus,
                                                  Min_Ip,
                                                  enableIT,
                                                  Force_NOIBR,
                                                  CTI,
                                                  OTmax,
                                                  type_select,
                                                  SetDir,
                                                  Sho_Plots,
                                                  GA_initial_seed=initpop)

    
