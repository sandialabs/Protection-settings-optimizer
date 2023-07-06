# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 15:26:55 2022

@author: maste
"""
import numpy as np
from ..src.GenNxGraph import index_dict

def pdef(mag,ang):
    A = mag*np.exp(1j*ang*(np.pi/180))
    return A

def abc2012(Vabc):
    a = 1*np.exp(120*(np.pi/180)*1j)
    #A = np.array([[1,1,1],[1,a**2,a],[1,a,a**2]])
    #Ainv = np.linalg.inv(A)
    
    V0 = (1/3) * (Vabc[0]+Vabc[1]+Vabc[2])
    V1 = (1/3) * (Vabc[0] + a*Vabc[1] + a**2*Vabc[2]) ;
    V2 = (1/3) * (Vabc[0] + a**2*Vabc[1] + a*Vabc[2]) ;
    return [V0,V1,V2]

# %% convert RONM system Data
def getRONMSysInfo(proSettings,powerFlow,ignore_fuse=True):
    """
    

    Parameters
    ----------
    proSettings : Dict
        JSON Dictionary contaning system data from RONM output
    ignore_fuse : Int
        Flag to ignore fuses.

    Returns
    -------
    SysInfo : Dict
        Dicttionary containg Bus, Line, transformer, Relay, Recloser, fuse, PV, Battery and generator data.

    """
    # Buses
    lenBuses = len(proSettings['network_model']['bus']) # Number of Buses
    Buses = [dict.fromkeys(['Name','nodes','numPhases','X','Y']) for number in range(lenBuses)]
    for ii in range(lenBuses):    # get Name, Nodes, numbr of phases and coordinates if avaliable for each bus
        Buses[ii]['Name'] = proSettings['network_model']['bus'][ii]['name']
        Buses[ii]['nodes'] = proSettings['network_model']['bus'][ii]['phases']
        Buses[ii]['numPhases'] = proSettings['network_model']['bus'][ii]['nphases']
        Buses[ii]['kV'] = np.mean([x['bus'][Buses[ii]['Name']]['voltage (V)'] for x in powerFlow if np.mean(x['bus'][Buses[ii]['Name']]['voltage (V)'])>0])
        if(np.isnan(Buses[ii]['kV'])):
            print('error in Bus data for '+ Buses[ii]['Name'])
            Buses[ii]['kV'] = 1.0
        if('X' in proSettings['network_model']['bus'][ii].keys()):
            Buses[ii]['X'] = proSettings['network_model']['bus'][ii]['X']
            Buses[ii]['Y'] = proSettings['network_model']['bus'][ii]['Y']
        else:
            Buses[ii]['X'] = np.nan
            Buses[ii]['Y'] = np.nan
    
    # Lines
    lenLines = len(proSettings['network_model']['line']) # Number of Lines
    Lines = [dict.fromkeys(['Name','Bus1','Bus2','Enabled','isSwitch','numPhases']) for number in range(lenLines)]
    for ii in range(lenLines):   # get Line names, connections, status , can be switched, number of phases, Length, impedance
        Lines[ii]['Name'] = proSettings['network_model']['line'][ii]['name']
        Lines[ii]['Bus1'] = proSettings['network_model']['line'][ii]['f_bus']
        Lines[ii]['Bus2'] = proSettings['network_model']['line'][ii]['t_bus']
        Lines[ii]['Enabled'] = proSettings['network_model']['line'][ii]['status']
        Lines[ii]['isSwitch'] = proSettings['network_model']['line'][ii]['switch']
        Lines[ii]['numPhases'] = proSettings['network_model']['line'][ii]['nphases']
        Lines[ii]['Length'] = 1
        #Lines[ii]['Rpu'] = proSettings['network_model']['line'][ii]['rs'][0][0]
        #Lines[ii]['Xpu'] = proSettings['network_model']['line'][ii]['xs'][0][0]

    # XFMRS
    lenXFMR = len(proSettings['network_model']['transformer']) # Number of transformers
    XFMRs = [dict.fromkeys(['Name','Bus1','Bus2','Enabled','numPhases']) for number in range(lenXFMR)]
    for ii in range(lenXFMR):   # get Name, Bus connections, status, number of phases and impedance 
        XFMRs[ii]['Name'] = proSettings['network_model']['transformer'][ii]['name']
        XFMRs[ii]['Bus1'] = proSettings['network_model']['transformer'][ii]['buses'][0]
        XFMRs[ii]['Bus2'] = proSettings['network_model']['transformer'][ii]['buses'][1]
        XFMRs[ii]['Enabled'] = proSettings['network_model']['transformer'][ii]['status']
        XFMRs[ii]['numPhases'] = proSettings['network_model']['transformer'][ii]['nphases']
        XFMRs[ii]['Rpu'] = 0;
        XFMRs[ii]['Xpu'] = 0;
    
    # Pvs, BESSs, Gens
    nPV = 0 # number of PV systems
    nBESS = 0 # number of bat systems
    nGENS = 0 # number of generators 
    for ii in range(len(proSettings['network_model']['source'])): # count the number of each source type
        if proSettings['network_model']['source'][ii]['type'].casefold()=='pvsystem':
            nPV+=1
        elif proSettings['network_model']['source'][ii]['type'].casefold()=='storage':
            nBESS+=1
        elif proSettings['network_model']['source'][ii]['type'].casefold()=='generator':
            nGENS+=1

    # Read PV Data
    if nPV == 0:    # if no PVs in systems return empty dict
        PVs = {}
    else:
        kk=0
        # get PV systme Name, Location, status and number of phases
        PVs = [dict.fromkeys(['Name','Bus','Enabled','numPhases']) for number in range(nPV)]
        for ii in range(len(proSettings['network_model']['source'])):
            if(proSettings['network_model']['source'][ii]['type'].casefold()=='pvsystem'):
                PVs[kk]['Name'] = proSettings['network_model']['source'][ii]['name']
                PVs[kk]['Bus'] = proSettings['network_model']['source'][ii]['bus']
                PVs[kk]['Enabled'] = proSettings['network_model']['source'][ii]['status']
                PVs[kk]['numPhases'] = proSettings['network_model']['source'][ii]['nphases']
                kk = kk + 1
            else:
                pass
    
    # Read BESS Data
    if nBESS == 0:    # if no storage in systems return empty dict
        BESSs = {}
    else:
        kk=0
        # get storage systme Name, Location, status and number of phases
        BESSs = [dict.fromkeys(['Name','Bus','Enabled','numPhases']) for number in range(nBESS)]
        for ii in range(len(proSettings['network_model']['source'])):
            if(proSettings['network_model']['source'][ii]['type'].casefold()=='storage'):
                BESSs[kk]['Name'] = proSettings['network_model']['source'][ii]['name']
                BESSs[kk]['Bus'] = proSettings['network_model']['source'][ii]['bus']
                BESSs[kk]['Enabled'] = proSettings['network_model']['source'][ii]['status']
                BESSs[kk]['numPhases'] = proSettings['network_model']['source'][ii]['nphases']
                kk = kk + 1
            else:
                pass
    
    # Read Generator Data
    if nGENS == 0:    # if no storage in systems return empty dict
        GENSs = {}
    else:
        kk=0
        # get generator Name, Location, status and number of phases
        GENSs = [dict.fromkeys(['Name','Bus','Enabled','numPhases']) for number in range(nGENS)]
        for ii in range(len(proSettings['network_model']['source'])):
            if(proSettings['network_model']['source'][ii]['type'].casefold()=='generator'):
                GENSs[kk]['Name'] = proSettings['network_model']['source'][ii]['name']
                GENSs[kk]['Bus'] = proSettings['network_model']['source'][ii]['bus']
                GENSs[kk]['Enabled'] = proSettings['network_model']['source'][ii]['status']
                GENSs[kk]['numPhases'] = proSettings['network_model']['source'][ii]['nphases']
                kk = kk + 1
            else:
                pass
    
    # Read Relays, Reclosers and fuses
    networkmodel = proSettings['network_model'] 
    nSRelays = 0
    nSRecs = 0
    nSFuses = 0
    for ii in range(len(networkmodel['protection'])):
        if(networkmodel['protection'][ii]['type'] == 'relay' and ('fuse_' not in networkmodel['protection'][ii]['name'])):
            nSRelays += 1
        elif(networkmodel['protection'][ii]['type'] == 'recloser'):
            nSRecs += 1
        elif(networkmodel['protection'][ii]['type'] == 'fuse'):
            nSFuses += 1
        else:
            print('Error:Device Type not found')
    
    # Read Relay Data
    if(nSRelays>0):
        kk = 0
        Relays = [dict.fromkeys(['Name','Bus1','Bus2','MonitoredObj','Enabled']) for number in range(nSRelays)]
        for ii in range(len(networkmodel['protection'])):
            if(networkmodel['protection'][ii]['type'] == 'relay'):
                Relays[kk]['Name'] = networkmodel['protection'][ii]['name']
                Relays[kk]['Enabled'] = 1
                Relays[kk]['MonitoredObj'] = networkmodel['protection'][ii]['location']
                Relays[kk]['Bus1'] = Lines[index_dict(Lines,'Name',networkmodel['protection'][ii]['location'].split('.')[1])]['Bus1']
                Relays[kk]['Bus2'] = Lines[index_dict(Lines,'Name',networkmodel['protection'][ii]['location'].split('.')[1])]['Bus2']
                kk+=1
    else:
        Relays = {}
    
    # Read Recloser Data
    if(nSRecs>0):
        kk = 0
        Recs = [dict.fromkeys(['Name','Bus1','Bus2','MonitoredObj','Enabled']) for number in range(nSRecs)]
        for ii in range(len(networkmodel['protection'])):
            if(networkmodel['protection'][ii]['type'] == 'recloser'):
                Recs[kk]['Name'] = networkmodel['protection'][ii]['name']
                Recs[kk]['Enabled'] = 1
                Recs[kk]['MonitoredObj'] = networkmodel['protection'][ii]['location']
                Recs[kk]['Bus1'] = Lines[index_dict(Lines,'Name',networkmodel['protection'][ii]['location'].split('.')[1])]['Bus1']
                Recs[kk]['Bus2'] = Lines[index_dict(Lines,'Name',networkmodel['protection'][ii]['location'].split('.')[1])]['Bus2']
                kk+=1
    else:
        Recs={}
    
    # Read Fuse Data Fuses
    if(nSFuses>0 and ignore_fuse==0):
        kk = 0
        Fuses = [dict.fromkeys(['Name','Bus1','Bus2','MonitoredObj','Enabled']) for number in range(nSFuses)]
        for ii in range(len(networkmodel['protection'])):
            if(networkmodel['protection'][ii]['type'] == 'fuse'):
                Fuses[kk]['Name'] = networkmodel['protection'][ii]['name']
                Fuses[kk]['Enabled'] = 1
                Fuses[kk]['MonitoredObj'] = networkmodel['protection'][ii]['location']
                Fuses[kk]['Bus1'] = Lines[index_dict(Lines,'Name',networkmodel['protection'][ii]['location'].split('.')[1])]['Bus1']
                Fuses[kk]['Bus2'] = Lines[index_dict(Lines,'Name',networkmodel['protection'][ii]['location'].split('.')[1])]['Bus2']
                kk+=1
    else:
        Fuses={}

    SysInfo = {"Relays": Relays,"Recs": Recs ,"Fuses": Fuses,"Lines": Lines, "XFMRs": XFMRs, "Buses": Buses, "Pvs": PVs, "BESS": BESSs,"Gens": GENSs}

    return SysInfo


# %% Convert RONM Fault Data to ADAPT format
def getRONMFaultData(Ts,faultBuses,devTypes,devNames,devLines,dev_BusV,Fault_Info):
    
    # list of recoreded lines 
    Line_list = list(Fault_Info[faultBuses[0]]['3p']['1']['line'].keys())
    switch_list = list(Fault_Info[faultBuses[0]]['3p']['1']['switch'].keys())
    
    Data_len = len(devNames)*len(faultBuses)*6
    
    kk = 0
    Fault_Data = [dict.fromkeys(['busNumber','Relay','FaultType','Ia_mag','Ia_ang','Ib_mag','Ib_ang','Ic_mag','Ic_ang','Va_mag','Va_ang','Vb_mag','Vb_ang','Vc_mag','Vc_ang','In_mag','In_ang','P','Q','Z012']) for number in range(Data_len)]
    for Fbus in faultBuses:
        # Fault location 
        for Ftype in list(Fault_Info[Fbus].keys()):
            # Fault type
            if(Ftype == '3p'):
                ft_str = 'TPH_'
            elif(Ftype == 'll'):
                ft_str = 'BC_'
            elif(Ftype == 'lg'):
                ft_str = 'SLG_A_'
            elif(Ftype == 'llg'):
                ft_str = 'DLG_'
            else:
                ft_str = 'Err';
                print('Error fault type unknown for '+ Fbus+Ftype)
            # Resistacne 
            for Fres in list(Fault_Info[Fbus][Ftype].keys()):
                Fstr = ft_str+'R'+Fres # fault type and resistance
                # Recored Data for all devices
                for dev_i in range(len(devLines)):
                    if(devLines[dev_i] in Line_list):
                        list_type = 'line' # if in line list look in line data
                    elif(devLines[dev_i] in switch_list):
                        list_type = 'switch' # if in switch list look in switch data
                    else:
                        print('Error: Fault Data not recorded, fault Bus '+ Fbus +' ,Type '+ Fbus+Ftype + 'Device on Line ' + devLines[dev_i] +' Device Name' + devNames[dev_i])
                        return ValueError
                    # collect Data
                    Fault_Data[kk]['busNumber'] = Fbus
                    Fault_Data[kk]['Relay'] = devNames[dev_i]
                    Fault_Data[kk]['FaultType'] = Fstr
                    # Curretns (A / deg)
                    if(len(Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['|V| (V)'])==3):
                        Fault_Data[kk]['Ia_mag'] = Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['|I| (A)'][0]
                        Fault_Data[kk]['Ia_ang'] = Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['theta'][0]
                        
                        Fault_Data[kk]['Ib_mag'] = Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['|I| (A)'][1]
                        Fault_Data[kk]['Ib_ang'] = Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['theta'][1]
                        
                        Fault_Data[kk]['Ic_mag'] = Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['|I| (A)'][2]
                        Fault_Data[kk]['Ic_ang'] = Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['theta'][2]
                        # voltages (pu)
                        Fault_Data[kk]['Va_mag'] = (Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['|V| (V)'][0]*1e3) / dev_BusV[dev_i]
                        Fault_Data[kk]['Va_ang'] = Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['phi (deg)'][0]
                        
                        Fault_Data[kk]['Vb_mag'] = (Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['|V| (V)'][1]*1e3) / dev_BusV[dev_i]
                        Fault_Data[kk]['Vb_ang'] = Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['phi (deg)'][1]
                        
                        Fault_Data[kk]['Vc_mag'] = (Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['|V| (V)'][2]*1e3) / dev_BusV[dev_i]
                        Fault_Data[kk]['Vc_ang'] = Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['phi (deg)'][2]

                    elif(len(Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['|V| (V)'])==2):
                        Fault_Data[kk]['Ia_mag'] = Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['|I| (A)'][0]
                        Fault_Data[kk]['Ia_ang'] = Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['theta'][0]
                        
                        Fault_Data[kk]['Ib_mag'] = Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['|I| (A)'][1]
                        Fault_Data[kk]['Ib_ang'] = Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['theta'][1]
                        
                        Fault_Data[kk]['Ic_mag'] = 0
                        Fault_Data[kk]['Ic_ang'] = 0
                        # voltages (pu)
                        Fault_Data[kk]['Va_mag'] = (Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['|V| (V)'][0]*1e3) / dev_BusV[dev_i]
                        Fault_Data[kk]['Va_ang'] = Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['phi (deg)'][0]
                        
                        Fault_Data[kk]['Vb_mag'] = (Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['|V| (V)'][1]*1e3) / dev_BusV[dev_i]
                        Fault_Data[kk]['Vb_ang'] = Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['phi (deg)'][1]
                        
                        Fault_Data[kk]['Vc_mag'] = 0
                        Fault_Data[kk]['Vc_ang'] = 0

                    else:
                        Fault_Data[kk]['Ia_mag'] = Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['|I| (A)'][0]
                        Fault_Data[kk]['Ia_ang'] = Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['theta'][0]
                        
                        Fault_Data[kk]['Ib_mag'] = 0
                        Fault_Data[kk]['Ib_ang'] = 0
                        
                        Fault_Data[kk]['Ic_mag'] = 0
                        Fault_Data[kk]['Ic_ang'] = 0
                        # voltages (pu)
                        Fault_Data[kk]['Va_mag'] = (Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['|V| (V)'][0]*1e3) / dev_BusV[dev_i]
                        Fault_Data[kk]['Va_ang'] = Fault_Info[Fbus][Ftype][Fres][list_type][devLines[dev_i]]['phi (deg)'][0]
                        
                        Fault_Data[kk]['Vb_mag'] = 0
                        Fault_Data[kk]['Vb_ang'] = 0
                        
                        Fault_Data[kk]['Vc_mag'] = 0
                        Fault_Data[kk]['Vc_ang'] = 0
                        
                    Va = pdef(Fault_Data[kk]['Va_mag'],Fault_Data[kk]['Va_ang'])*dev_BusV[dev_i]
                    Vb = pdef(Fault_Data[kk]['Vb_mag'],Fault_Data[kk]['Vb_ang'])*dev_BusV[dev_i]
                    Vc = pdef(Fault_Data[kk]['Vc_mag'],Fault_Data[kk]['Vc_ang'])*dev_BusV[dev_i]
                    
                    Ia = pdef(Fault_Data[kk]['Ia_mag'],Fault_Data[kk]['Ia_ang'])
                    Ib = pdef(Fault_Data[kk]['Ib_mag'],Fault_Data[kk]['Ib_ang'])
                    Ic = pdef(Fault_Data[kk]['Ic_mag'],Fault_Data[kk]['Ic_ang'])
                    
                    V012 = abc2012([Va,Vb,Vc])
                    I012 = abc2012([Ia,Ib,Ic])
                    
                    Fault_Data[kk]['Vabc'] = [Va,Vb,Vc]
                    Fault_Data[kk]['Iabc'] = [Ia,Ib,Ic]
                    Fault_Data[kk]['V012'] = V012
                    Fault_Data[kk]['I012'] = I012
                    
                    Fault_Data[kk]['In_mag'] = abs(Fault_Data[kk]['I012'][0])
                    Fault_Data[kk]['In_ang'] = np.angle(Fault_Data[kk]['I012'][0],deg=True)
                    
                    Sa = Va * np.conj(Ia)
                    Sb = Vb * np.conj(Ib)
                    Sc = Vc * np.conj(Ic)
                    
                    Fault_Data[kk]['P'] = np.real(Sa+Sb+Sc)
                    Fault_Data[kk]['Q'] = np.imag(Sa+Sb+Sc)
                    Fault_Data[kk]['Z012'] = [np.divide(V012[0],I012[0],where=I012[0]!=0),
                                                  np.divide(V012[1],I012[1],where=I012[1]!=0),
                                                  np.divide(V012[2],I012[2],where=I012[2]!=0)]
                    kk = kk+1
    return Fault_Data[:kk]


# %% Convert RONM Device Data to ADAPT format
def getRONMDeviceData(Ts,powerFlow,devTypes,devNames,devLines,dev_BusV,SysInfo,sW_Status,sW_Names):
    Device_Data_CSV = [dict.fromkeys(['RelayName','Ia_mag','Ia_ang','Va_mag','Va_ang','P','Q','switchState','PVs','In']) for number in range(len(devNames))]
    for dev_i in range(len(devNames)):
        Device_Data_CSV[dev_i]['RelayName'] = devNames[dev_i] # RelayName
        
        # Find Switch State
        Line_Num = index_dict(SysInfo['Lines'],'Name',devLines[dev_i])
        Switch_Num = sW_Names.index(devNames[dev_i]) if devNames[dev_i] in sW_Names else None
        
        if(Line_Num == None and Switch_Num == None):
            Device_Data_CSV[dev_i]['switchState'] =  0
        elif(Line_Num != None):
            Device_Data_CSV[dev_i]['switchState'] = SysInfo['Lines'][Line_Num]['Enabled']
            if(Switch_Num != None):
                Device_Data_CSV[dev_i]['switchState'] = sW_Status[Ts][Switch_Num]
            else:
                Device_Data_CSV[dev_i]['switchState'] = 0
        
        Vmag = [x*1e3 for x in powerFlow[Ts]['protection'][devTypes[dev_i]+'.'+devNames[dev_i]]['voltage (kV)']] 
        Vang = powerFlow[Ts]['protection'][devTypes[dev_i]+'.'+devNames[dev_i]]['phi (deg)']
        V = [None]*len(Vmag)*2
        I = [None]*len(Vmag)*2
        
        P = [x*1e3 for x in powerFlow[Ts]['protection'][devTypes[2]+'.'+devNames[2]]['real power flow (kW)'] ]
        Q = [x*1e3 for x in powerFlow[Ts]['protection'][devTypes[2]+'.'+devNames[2]]['reactive power flow (kVar)'] ]
        S = [x+y*1j for x,y in zip(P,Q)]
        
        for ii in range(len(Vmag)):
            V[ii*2] = Vmag[ii]
            V[ii*2+1] = Vang[ii]
            I_calc = np.conj(S[ii]/pdef(Vmag[ii],Vang[ii]))
            I[ii*2] = abs(I_calc)
            I[ii*2+1] = np.angle(I_calc)*(180/np.pi)
        
        Imag = I[0::2]
        
        Device_Data_CSV[dev_i]['Va_mag'] = min(Vmag) # Va_mag
        Device_Data_CSV[dev_i]['Va_ang'] = Vang[0] # Va_ang
        # max Phase current 
        Device_Data_CSV[dev_i]['Ia_mag'] = max(Imag)
        Device_Data_CSV[dev_i]['Ia_ang'] = I[1]
        
        if(len(V) == 6): # 3ph device
            pha = [V[0],V[1],
                   V[2],V[3],
                   V[4],V[5],
                   I[0],I[1],
                   I[2],I[3],
                   I[4],I[5] ]
        elif(len(V) == 4): #2 phase device
            pha = [V[0],V[1],
                   V[2],V[3],
                   0,0,
                   I[0],I[1],
                   I[2],I[3],
                   0,0]
        elif(len(V) == 2): # single phase device
            pha = [V[0]/dev_BusV[dev_i],V[1],
                   0,0,
                   0,0,
                   I[0],I[1],
                   0,0,
                   0,0 ]

        Device_Data_CSV[dev_i]['Vabc'] = [pdef(pha[0],pha[1]),pdef(pha[2],pha[3]),pdef(pha[4],pha[5])]
        Device_Data_CSV[dev_i]['Iabc'] = [pdef(pha[6],pha[7]),pdef(pha[8],pha[9]),pdef(pha[10],pha[11])]
        Device_Data_CSV[dev_i]['V012'] = abc2012(Device_Data_CSV[dev_i]['Vabc'])
        Device_Data_CSV[dev_i]['I012'] = abc2012(Device_Data_CSV[dev_i]['Iabc'])
        Device_Data_CSV[dev_i]['Z012'] = np.divide(Device_Data_CSV[dev_i]['V012'],Device_Data_CSV[dev_i]['I012'])
        
        Device_Data_CSV[dev_i]['In'] = abs(pdef(pha[6],pha[7]) + pdef(pha[8],pha[9]) + pdef(pha[10],pha[11]))
        
        Sa = pdef(pha[0],pha[1]) * np.conj( pdef(pha[6],pha[7]) )
        Sb = pdef(pha[2],pha[3]) * np.conj( pdef(pha[8],pha[9]) )
        Sc = pdef(pha[4],pha[5]) * np.conj( pdef(pha[10],pha[11]) )
        
        Device_Data_CSV[dev_i]['P'] = np.real(Sa+Sb+Sc)
        Device_Data_CSV[dev_i]['Q'] = np.imag(Sa+Sb+Sc)
        Device_Data_CSV[dev_i]['PVs'] = 1
        
    return Device_Data_CSV


# %% Get Switch states
def get_Sw_Status(powerFlow,devTimeLine):
    # Get Swtich Names
    sW_Names=list(powerFlow[0]['switch'].keys())
    sw_cols = len(sW_Names)
    sw_rows = len(devTimeLine)
    # Get Switch Status for each switch for each timestep 
    sW_Status = np.empty((sw_rows,sw_cols,))
    for ii in range(sw_rows):
        for jj in range(sw_cols):
            sWstate = devTimeLine[ii]['Switch configurations'][sW_Names[jj]]
            if sWstate.casefold()=='open':
                sW_Status[ii][jj] = 0
            elif sWstate.casefold()=='closed':
                sW_Status[ii][jj] = 1
            else:
                sW_Status[ii][jj] = 2
    return (sW_Status,sW_Names)