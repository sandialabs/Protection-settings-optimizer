# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 12:40:32 2022

@author: maste
"""
# CAPE inpit method
import os
import json
import time
import pandas
import RSO_pack

startTime = time.time()
pwd = os.getcwd()
Main_dir = pwd

# %% Program settigns 
Force_NOIBR = 0
DOC = 1
enableIT = 0
CTI = 0.25
OTmax = 10
Sho_Plots=False

# %% Input Data files
Substation_bus = 'sourcebus'
initpop = None
Device_Data_File = os.path.join(pwd,'Examples','Adapt_Data','outputV1_624.csv')

Fault_Data_File = os.path.join(pwd,'Examples','Adapt_Data','fault_report_624.csv')

# load system data
jsonFile = os.path.join(Main_dir,'Examples','IEEE123_SysInfo.json') 
f = open(jsonFile) 
IEEE123_SysInfo = json.load(f)

# %% Process input files

switchStates = [1,1,1,1]
switchLines = ['sw2','sw3','sw7','sw4']

Device_Data_CSV = RSO_pack.read_Relays_CSV_Data(Device_Data_File)
#               RTL1 RTL2 RTL3 RTL4
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
# get fault data from file 
Fault_Data_CSV = RSO_pack.read_Fault_CSV_Data(Fault_Data_File)



# %% calculate settings
type_select = True
SetDir=True
Min_Ip = [60,12]
Fault_Res = list(set(['R'+x['FaultType'].split('_R')[1] for x in Fault_Data_CSV]))

settings,old_info = RSO_pack.runSettingsOptimizer(Main_dir,
                                                  switchStates,
                                                  switchLines,
                                                  Device_Data_CSV,
                                                  Fault_Data_CSV,
                                                  Fault_Res,
                                                  IEEE123_SysInfo,
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

DevSet = pandas.DataFrame(settings)
print(DevSet)
# Write setting to file
outfile_name = 'Settings_DT_V2.csv'
# %%
executionTime = (time.time() - startTime)
print(executionTime)
RSO_pack.writeSettingsToExcel(old_info,IEEE123_SysInfo,settings,outfile_name)
        


