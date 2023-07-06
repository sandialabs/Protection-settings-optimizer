# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 14:36:15 2022

@author: maste
"""
import pandas as pa
import numpy as np
from .GenNxGraph import index_dict

def writeSettingsToExcel(old_info,IEEE123_SysInfo,settings,outfile_name):
    RName = ['R1','RTL1','R2','R3','RTL2','R4','RTL4','R5','R6','RTL3']
    RIp_max =[1000,1536.0,1000,1000,1024.0,1000,1000.0,1000,1000,1000.0]
    if old_info['DOC'] == 1:
        Name = [x for xs in [[x+'_F',x+'_R'] for x in RName] for x in xs]
        Ip_max = [x for xs in [[x,x] for x in RIp_max] for x in xs]
    else:
        Name = [x+'_E' for x in RName]
        Ip_max = RIp_max
    
    settings_table = [None]*len(Name)
    for ii in range(len(Name)):
        ind = index_dict(settings,'Name',Name[ii])
        if ind == None:
            Name_ = Name[ii].split('_')[0]
            Dir_ = Name[ii].split('_')[1]
            ind_FT = index_dict(IEEE123_SysInfo['Relays'],'Name',Name_)
            if(Dir_ == 'R'):
                To = IEEE123_SysInfo['Relays'][ind_FT]['Bus1']
                From = IEEE123_SysInfo['Relays'][ind_FT]['Bus2']
            else:
                From = IEEE123_SysInfo['Relays'][ind_FT]['Bus1']
                To = IEEE123_SysInfo['Relays'][ind_FT]['Bus2']
            
            if old_info['DOC'] == 0:
                settings_table[ii] ={'Name':Name[ii],'From':From,'To':To,'PickupI':Ip_max[ii],
                                    'TDS':1,'TOC':'U1: moderately inverse (OC)','PickupI0':Ip_max[ii],'TDSg':1,'TOCg':'U1: moderately inverse (OC)',
                                     'VR': str(np.nan),'IT': str(np.nan),'IT0' : str(np.nan)}
            else:
                settings_table[ii] ={'Name':Name[ii],'From':From,'To':To,'PickupI':Ip_max[ii],
                                    'TDS':1,'TOC':'U1: moderately inverse (OC)','PickupI0':Ip_max[ii],'TDSg':1,'TOCg':'U1: moderately inverse (OC)',
                                     'VR': str(np.nan),'IT': str(np.nan),'IT0' : str(np.nan),'Z1MAG':1,'Z1ANG': 45,'Z0MAG':1,'Z0ANG':45,'Z2F':0.1,'Z2R':0.1}
        else:
            settings_table[ii] = settings[ind]
            settings_table[ii]['Name'] = Name[ind]
            # edit VR
            if(settings_table[ii]['VR'] == False):
                settings_table[ii]['VR'] = str(np.nan)
            else:
                settings_table[ii]['VT'] = 0.95
            # edit IT 
            if(settings_table[ii]['IT'] == 0):
                settings_table[ii]['IT'] = str(np.nan)
            # edit IT0 
            if(settings_table[ii]['IT0'] == 0):
                settings_table[ii]['IT0'] = str(np.nan)

    output = pa.DataFrame(settings_table)
    output.to_csv(outfile_name,index=False)