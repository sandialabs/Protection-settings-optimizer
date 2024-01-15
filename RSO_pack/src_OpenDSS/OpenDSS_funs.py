# -*- coding: utf-8 -*-
"""
Functions to interface wiht opendss
"""

import numpy as np
import pandas as pd
import opendssdirect as dss

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

def getDeviceData(dssCircuit,devNames,devLines,dev_BusV):
    # dict of devices
    Device_Data_CSV = [dict.fromkeys(['RelayName','Ia_mag','Ia_ang','Va_mag','Va_ang','P','Q','switchState','PVs','In']) for number in range(len(devNames))]
    # Vmag,Imag,In,P,V,I
    for dev_i in range(len(devNames)):
        SS_Data = getLineVI(dssCircuit, devLines[dev_i])
        Device_Data_CSV[dev_i]['RelayName'] = devNames[dev_i]
        # voltage
        V = SS_Data[5]
        I = SS_Data[6]
        
        Vmin_ind = SS_Data[0].index(min(SS_Data[0]))
        Imax_ind = SS_Data[1].index(max(SS_Data[1]))
        
        Device_Data_CSV[dev_i]['Va_mag'] = SS_Data[0][Vmin_ind]
        Device_Data_CSV[dev_i]['Va_ang'] = V[(Vmin_ind*2)+1]
        # current 
        Device_Data_CSV[dev_i]['Ia_mag'] = SS_Data[1][Imax_ind]
        Device_Data_CSV[dev_i]['Ia_ang'] = I[(Imax_ind*2)+1]
        
        if(len(V)==12):
            pha = [V[0],V[1],
                   V[2],V[3],
                   V[4],V[5],
                   I[0],I[1],
                   I[2],I[3],
                   I[4],I[5] ]
        elif(len(V)==8):
            pha = [V[0],V[1],
                   V[2],V[3],
                   0,0,
                   I[0],I[1],
                   I[2],I[3],
                   0,0]
        elif(len(V) == 4):
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
        Device_Data_CSV[dev_i]['switchState'] = True
    
    return Device_Data_CSV
        

def getFaultInfo(dssCircuit,dssText,faultBuses,faultBusPhases,Fres,Fts,devLines,devNames,dev_BusV):
    """
    

    Parameters
    ----------
    dssCircuit : OpenDSS circuit Object
        COM interface to OpenDSS circuit
    dssText : OpenDSS circuit Object
        Text interface to OpenDSS
    faultBuses : List
        list of fault bus names.
    faultBusPhases : List
        number of Phases for each fault bus.
    Fres : List
        List of fault resistacnes to test.
    Fts : List
        Name of Fault types, supported types ('3ph','LL','SLG').
    devLines : List
        List of lines with protection devices.
    devNames : List
        Name of protection devices corresponding to the device list.
    dev_BusV : List 
        Rated voltage (L-N) of the protection device.

    Returns
    -------
    FData : Pandas DataFrame
        Datafram contining the fault current observed by devices for each fault.

    """
    dss.Circuit = dssCircuit

    dssText.command = 'set maxcontroliter = 500'
    dssText.Command = 'solve'

    #create fault to be moved around
    dssText.command = 'New Fault.F1 enabled=false'

    FData = pd.DataFrame()
    for res in Fres:
        for Ft in Fts:
            for ii in range(len(faultBuses)):
                NewData = simFaults(faultBuses[ii],Ft,faultBusPhases[ii],res,devLines,devNames,dev_BusV,dssCircuit,dssText)
                if(type(NewData)!=int):
                    FData = pd.concat([FData,NewData],)
                else:
                    print('Cannot Get '+Ft+' fault data for bus '+faultBuses[ii])
    return FData

# %% greate systeminfo json file
def getSysInfo(dssCircuit):
    """
    Create system Json from Opendss

    Parameters
    ----------
    dssCircuit : OpenDSS circuit object
        OpenDSS circuit object wiht the circuit compiled

    Returns
    -------
    SysInfo : dict
        dict contaning system information (Relays,Reclosers,Fuses,Lines,Transformers,Buses, Pvs .

    """
    dss.Circuit = dssCircuit
    # collect network data
    Relays = getRelayInfo(dssCircuit)
    Recs = getRecloserInfo(dssCircuit)
    Fuses = getFuseInfo(dssCircuit)
    Lines = getLineInfo(dssCircuit)
    XFMRs = getTransformerInfo(dssCircuit)
    Buses = getBusInfo(dssCircuit)
    Pvs = getPvInfo(dssCircuit)
    BESS = getBESSInfo(dssCircuit)
    Gens = getGeneratorInfo(dssCircuit)
    
    
    SysInfo = {"Relays": Relays,"Recs": Recs ,"Fuses": Fuses,"Lines": Lines, "XFMRs": XFMRs, "Buses": Buses, "Pvs": Pvs, "BESS": BESS,"Gens": Gens}
    return SysInfo


# %% get battery info 
def getBESSInfo(dssCircuit):
    """
    Get Battery Energy systems information from OpenDSS

    Parameters
    ----------
    dssCircuit : OpenDSS circuit Object
        A compiled circuit in OpenDSS

    Returns
    -------
    BESSs : List
        List of Battery energy systems.

    """
    dssCircuit.SetActiveClass('Storage')    # set class in opendss
    BESS_names = dss.ActiveClass.AllNames() # gent list of bat systems
    nBESS = len(BESS_names)
    if nBESS == 0:      # if no Battery systems return empty dict
        BESSs = {}; 
        return BESSs
    else:               # get Battery system Name, Location, status and number of phases 
        BESSs = [dict.fromkeys(['Name','Bus','Enabled','numPhases']) for number in range(nBESS)]
        for ii in range(nBESS):
            dssCircuit.SetActiveElement('Storage.'+BESS_names[ii])
            BESSs[ii]['Name'] = dss.CktElement.Name().split('.')[1]
            BESSs[ii]['Bus'] = dss.CktElement.BusNames()[0]
            BESSs[ii]['Enabled'] = dss.CktElement.Enabled()
            BESSs[ii]['numPhases'] = dss.Properties.Value('phases')
        return(BESSs)

# %% get Generator info 
def getGeneratorInfo(dssCircuit):
    """
    Get generator information from OpenDSS

    Parameters
    ----------
    dssCircuit : OpenDSS circuit Object
        A compiled circuit in OpenDSS

    Returns
    -------
    Gens : List
        List of generators in the systems.

    """
    Gen_names = dss.Generators.AllNames() # get list of generators in the system
    if len(Gen_names) == 0 or Gen_names[0] == 'NONE': # if no generators in systems return empty dict
        Gens = {}
        return(Gens)
    else:
        # get generator Name, Location, status and number of phases 
        Gens = [dict.fromkeys(['Name','Bus','Enabled','numPhases']) for number in range(len(Gen_names))]
        for ii in range(len(Gen_names)):
            dssCircuit.SetActiveElement('generator.'+Gen_names[ii])
            Gens[ii]['Name'] = dss.CktElement.Name.split()('.')[1]
            Gens[ii]['Bus'] = dss.CktElement.BusNames()[0]
            Gens[ii]['Enabled'] = dss.CktElement.Enabled()
            Gens[ii]['numPhases'] = dss.Properties.Value('phases')
        return(Gens)

# %% get PV info 
def getPvInfo(dssCircuit):
    """
    Get PV systems information from OpenDSS

    Parameters
    ----------
    dssCircuit : OpenDSS circuit Object
        A compiled circuit in OpenDSS

    Returns
    -------
    PVs : List
        List of PVs in the systems.

    """
    dssCircuit.SetActiveClass('PVSystem')
    Pv_names = dss.PVsystems.AllNames()    # get list of PVs in the system names
    nPV = len(Pv_names)
    if nPV == 0:    # if no PVs in systems return empty dict
        PVs = {}
        return(PVs)
    else:          # get PV systme Name, Location, status and number of phases
        PVs = [dict.fromkeys(['Name','Bus','Enabled','numPhases']) for number in range(len(Pv_names))]
        for ii in range(len(Pv_names)):
            dssCircuit.SetActiveElement('PVSystem.'+Pv_names[ii])
            PVs[ii]['Name'] = dss.CktElement.Name().split('.')[1]
            PVs[ii]['Bus'] = dss.CktElement.BusNames()[0]
            PVs[ii]['Enabled'] = dss.CktElement.Enabled()
            PVs[ii]['numPhases'] = dss.Properties.Value('phases')
            return(PVs)

# %% get Bus info 
def getBusInfo(dssCircuit):
    """
    Get Bus information from OpenDSS

    Parameters
    ----------
    dssCircuit : OpenDSS circuit Object
        A compiled circuit in OpenDSS

    Returns
    -------
    Buses : List
        List of Buses in the systems.

    """
    Bus_names = dssCircuit.AllBusNames() # get all bus names
    Buses = [dict.fromkeys(['Name','nodes','numPhases','X','Y']) for number in range(len(Bus_names))]
    for ii in range(len(Bus_names)):    # get Name, Nodes, numbr of phases and coordinates if avaliable for each bus
        dssCircuit.SetActiveBus(Bus_names[ii])
        Buses[ii]['Name'] = dss.Bus.Name().split('.')[0]
        Buses[ii]['nodes'] = dss.Bus.Nodes()
        Buses[ii]['numPhases'] = dss.Bus.NumNodes()
        Buses[ii]['kV'] = dss.Bus.kVBase()
        if(dss.Bus.Coorddefined()):
            Buses[ii]['X'] = dss.Bus.X()
            Buses[ii]['Y'] = dss.Bus.Y()
        else:
            Buses[ii]['X'] = np.nan
            Buses[ii]['Y'] = np.nan
    return(Buses)

# %% Get Transformer information 
def getTransformerInfo(dssCircuit):
    """
    Get Transformer Information  from OpenDSS

    Parameters
    ----------
    dssCircuit : OpenDSS circuit Object
        A compiled circuit in OpenDSS

    Returns
    -------
    XFMRs : List
        List of Transformers in the systems.

    """
    XFMR_names = dss.Transformers.AllNames() # get all transformer names from opendss
    XFMRs = [dict.fromkeys(['Name','Bus1','Bus2','Enabled','numPhases']) for number in range(len(XFMR_names))]
    
    for ii in range(len(XFMR_names)):   # get Name, Bus connections, status, number of phases and impedance 
        dssCircuit.SetActiveElement('transformer.'+XFMR_names[ii])
        XFMRs[ii]['Name'] = dss.CktElement.Name().split('.')[1]
        XFMRs[ii]['Bus1'] = dss.CktElement.BusNames()[0]
        XFMRs[ii]['Bus2'] = dss.CktElement.BusNames()[1]
        XFMRs[ii]['Enabled'] = dss.CktElement.Enabled()
        XFMRs[ii]['numPhases'] = dss.Properties.Value('phases')
        XFMRs[ii]['Rpu'] = float(dss.Properties.Value('%R'))
        XFMRs[ii]['Xpu'] = float(dss.Properties.Value('X12'))
    return(XFMRs)

# %% get Line info
def getLineInfo(dssCircuit):
    """
    Get Line Information  from OpenDSS

    Parameters
    ----------
    dssCircuit : OpenDSS circuit Object
        A compiled circuit in OpenDSS

    Returns
    -------
    Lines : List
        List of Lines in the systems.

    """
    Line_names = dss.Lines.AllNames() # get all line Names
    Lines = [dict.fromkeys(['Name','Bus1','Bus2','Enabled','isSwitch','numPhases']) for number in range(len(Line_names))]
    
    for ii in range(len(Line_names)):   # get Line names, connections, status , can be switched, number of phases, Length, impedance
        dssCircuit.SetActiveElement('line.'+Line_names[ii])
        Lines[ii]['Name'] = dss.CktElement.Name().split('.')[1]
        Lines[ii]['Bus1'] = dss.CktElement.BusNames()[0]
        Lines[ii]['Bus2'] = dss.CktElement.BusNames()[1]
        Lines[ii]['Enabled'] = dss.CktElement.Enabled()
        Lines[ii]['isSwitch'] = dss.Properties.Value('Switch')
        Lines[ii]['numPhases'] = dss.Properties.Value('phases')
        Lines[ii]['Length'] = float(dss.Properties.Value('length'))
        Lines[ii]['Rpu'] = dss.Lines.RMatrix()[0] * Lines[ii]['Length']
        Lines[ii]['Xpu'] = dss.Lines.XMatrix()[0] * Lines[ii]['Length']
    return(Lines)

# %% Get Fuse info 
def getFuseInfo(dssCircuit):
    """
    

    Parameters
    ----------
    dssCircuit : OpenDSS circuit Object
        A compiled circuit in OpenDSS

    Returns
    -------
    Fuses : List
        List of fuses in the system 

    """
    Fuse_names = dss.Fuses.AllNames() # get all fuse names from OpenDSS
    if(len(Fuse_names) == 0 or Fuse_names[0] == 'NONE'): # return empty dict if no fuses in the system 
        Fuses = {}
        return Fuses
    else: # get fuse Name, status, Type, Rating, location
        Fuses = [dict.fromkeys(['Name','Bus1','Bus2','MonitoredObj','Enabled','Type','Rating']) for number in range(len(Fuse_names))]
        for ii in range(len(Fuse_names)):
            dssCircuit.SetActiveElement('Fuse.'+Fuse_names[ii])
            Fuses[ii]['Name'] = dss.CktElement.Name().split('.')[1]
            Fuses[ii]['Enabled'] = dss.CktElement.Enabled()
            Fuses[ii]['Type'] = dss.Properties.Value('FuseCurve')
            Fuses[ii]['Rating'] = dss.Properties.Value('RatedCurrent')
            Fuses[ii]['MonitoredObj'] = dss.Properties.Value('MonitoredObj')
            #change Active Circuit to Fuse and get bus info
            dssCircuit.SetActiveElement(Fuses[ii]['MonitoredObj']) 
            Fuses[ii]['Bus1'] = dss.CktElement.BusNames()[0].split('.')[0]
            Fuses[ii]['Bus2'] = dss.CktElement.BusNames()[1].split('.')[0]
            
        return Fuses
    
# %% Get Rec info
def getRecloserInfo(dssCircuit):
    """
    

    Parameters
    ----------
    dssCircuit : OpenDSS circuit Object
        A compiled circuit in OpenDSS

    Returns
    -------
    Rec : List
        List of reclosers in the system 

    """
    Rec_names = dss.Reclosers.AllNames() # get all recloser names
    if(len(Rec_names) == 0 or Rec_names[0] == 'NONE'): # if no reclosers in system retutn empty dict
        Rec = {}
        return Rec
    else:
        # get recloser name, location, status from opendss
        Rec = [dict.fromkeys(['Name','Bus1','Bus2','MonitoredObj','Enabled']) for number in range(len(Rec_names))]
        for ii in range(len(Rec_names)):
            dssCircuit.SetActiveElement('Recloser.'+Rec_names[ii])
            Rec[ii]['Name'] = dss.CktElement.Name().split('.')[1]
            Rec[ii]['Enabled'] = dss.CktElement.Enabled()
            Rec[ii]['MonitoredObj'] = dss.Properties.Value('MonitoredObj')
            # change Active Circuit to Reclosers 
            dssCircuit.SetActiveElement(Rec[ii]['MonitoredObj'])
            Rec[ii]['Bus1'] = dss.CktElement.BusNames()[0].split('.')[0]
            Rec[ii]['Bus2'] = dss.CktElement.BusNames()[1].split('.')[0]
        return Rec # return list of reclosers

# %% Get Relay info
def getRelayInfo(dssCircuit):
    """
    Get relay info from opendss

    Parameters
    ----------
    dssCircuit : OpenDSS circuit Object
        A compiled circuit in OpenDSS

    Returns
    -------
    Relays : TYPE
        List of relays in the system.

    """
    Relay_names = dss.Relays.AllNames()
    if(Relay_names[0] == 'NONE'): # if no relays in the system return empty dict
        Relays = {}
        return Relays
    else:
        # get relay name, location, status for each relay in the system 
        Relays = [dict.fromkeys(['Name','Bus1','Bus2','MonitoredObj','Enabled']) for number in range(len(Relay_names))]
        for ii in range(len(Relay_names)):
            dssCircuit.SetActiveElement('relay.'+Relay_names[ii])
            Relays[ii]['Name'] = dss.CktElement.Name().split('.')[1]
            Relays[ii]['Enabled'] = dss.CktElement.Enabled()
            Relays[ii]['MonitoredObj'] = dss.Properties.Value('MonitoredObj')
            # change Active Circuit to Relay 
            dssCircuit.SetActiveElement(Relays[ii]['MonitoredObj'])
            Relays[ii]['Bus1'] = dss.CktElement.BusNames()[0].split('.')[0]
            Relays[ii]['Bus2'] = dss.CktElement.BusNames()[1].split('.')[0]
    return(Relays) # return list of relays

# %% Get Load info 
def getLoadInfo(dssCircuit):
    """
    Get Load info

    Parameters
    ----------
    dssCircuit : OpenDSS circuit Object
        A compiled circuit in OpenDSS

    Returns
    -------
    Loads : list of dict
        List of loads in the system.

    """
    Load_names = dssCircuit.Loads.AllNames # get all load names 
    if(Load_names[0] == 'NONE'): # if no loads in the system return empty dict
        Loads = {}
        return Loads
    else:
        # get load name, status, location, connection (nodes and phases), voltage, powerfactor, powers, model type, connection type (wye or delta
        Loads = [dict.fromkeys(['Name','Bus','Nodes','Enabled','phases','kV','pf','kW','kvar','kVA','model','conn']) for number in range(len(Load_names))]
        for ii in range(len(Load_names)):
            dssCircuit.setActiveElement('load.'+Load_names[ii])
            Loads[ii]['Name'] = dss.CktElement.Name().split('.')[1]
            Loads[ii]['Enabled'] = dss.CktElement.Enabled()
            Loads[ii]['Bus'] = dss.CktElement.BusNames()[0].split('.')[0]
            Loads[ii]['Nodes'] = dss.CktElement.NodeOrder()
            Loads[ii]['phases'] = dss.CktElement.Properties.Value('phases')
            Loads[ii]['kV'] = dss.CktElement.Properties.Value('kV')
            Loads[ii]['pf'] = dss.CktElement.Properties.Value('pf')
            Loads[ii]['kW'] = dss.CktElement.Properties.Value('kW')
            Loads[ii]['kvar'] = dss.CktElement.Properties.Value('kvar')
            Loads[ii]['kVA'] = dss.CktElement.Properties.Value('kVA')
            Loads[ii]['model'] = dss.CktElement.Properties.Value('model')
            Loads[ii]['conn'] = dss.CktElement.Properties.Value('model')
            
        return  Loads

# %% Get Load Data

def getLoadVI(dssCircuit,Load):
    """
    Get voltage, current and power info form a specified load

    Parameters
    ----------
    dssCircuit : OpenDSS circuit Object
        A compiled and solved circuit in OpenDSS
    Load : str
        Load Name

    Returns
    -------
    Vmag : List
        Magnitude of Vabc.
    Imag : List
        Magnitude of Iabc
    P : List
        Real power Pabc
    Q : List
        Reactive power Qabc
    V : List
        Vabc phasors
    I : List
        Iabc phasors

    """
    dssCircuit.SetActiveElement('load.'+Load)
    nNodes = 2*int(dss.CktElement.Properties('phases').val)
    V = dss.CktElement.VoltagesMagAng()
    I = dss.CktElement.CurrentsMagAng()
    S = dss.CktElement.Powers()
    dssCircuit.SetActiveBus(dssCircuit.ActiveCktElement.BusNames[0])
    Vbase = dss.Bus.kVBase()*1000
    P = [0]*int(nNodes/2)
    Q = [0]*int(nNodes/2)
    Vmag = [0]*int(nNodes/2)
    Imag = [0]*int(nNodes/2)
    for ii in range(0,nNodes,2):
        Vmag[int(ii/2)] = V[ii]/Vbase
        Imag[int(ii/2)] = I[ii]
        P[int(ii/2)] = S[ii]
        Q[int(ii/2)] = S[ii+1]
    return (Vmag,Imag,P,Q,V,I)

# %% Get Line Data
def getLineVI(dssCircuit,Line):
    """
    

    Parameters
    ----------
    dssCircuit : OpenDSS circuit Object
        A compiled and solved circuit in OpenDSS
    Line : TYPE
        Line Name

    Returns
    -------
    Vmag : List
        Magnitude of Vabc.
    Imag : List
        Magnitude of Iabc
    In : complex Float
        neutral current .
    P : List
        Real power for each phase.
    V : List
        Vabc phasors.
    I : List
        Iabc phasors.

    """
    dssCircuit.SetActiveElement('line.'+Line)
    nNodes = len(dss.CktElement.NodeOrder())
    V = dss.CktElement.VoltagesMagAng()
    I = dss.CktElement.CurrentsMagAng()
    S = dss.CktElement.Powers()
    Iseq = dss.CktElement.SeqCurrents()
    dssCircuit.SetActiveBus(dss.CktElement.BusNames()[0])
    Vbase = dss.Bus.kVBase()*1000
    P = [0]*int(nNodes/2)
    Q = [0]*int(nNodes/2)
    Vmag = [0]*int(nNodes/2)
    Imag = [0]*int(nNodes/2)
    for ii in range(0,nNodes,2):
        Vmag[int(ii/2)] = V[ii]/Vbase
        Imag[int(ii/2)] = I[ii]
        P[int(ii/2)] = S[ii]
        Q[int(ii/2)] = S[ii+1]
    In = 3*Iseq[0] 
    
    return (Vmag,Imag,In,P,Q,V,I)       

# %% simFaults 
def simFaults(faultBuses,Type,faultBusPhases,Fres,devLines,devNames,dev_BusV,dssCircuit,dssText):
    data = [[None]*17 for i in devLines]
    #for ii in range(len(faultBuses)):
    if(Type == '3ph' and len(faultBusPhases)>2):
        Ftype = 'TPH'+'_R'+Fres
        Ftype1 = '.1.2.3 phases=3 '
        Ftype2 = '.0.0.0 '

    elif(Type == 'LL' and len(faultBusPhases)>1):
        Ftype = 'BC'+'_R'+Fres
        #Ftype1 = '.2 phases=1 '
        #Ftype2 = '.3 '
        Ftype1 = '.'+str(faultBusPhases[len(faultBusPhases)-2])+' phases=1 '
        Ftype2 = '.'+str(faultBusPhases[len(faultBusPhases)-1])+' phases=1 '

    elif(Type == 'SLG' and len(faultBusPhases)>0):
        Ftype = 'SLG_A'+'_R'+Fres
        Ftype1 = '.'+str(faultBusPhases[0])+' phases=1 '
        Ftype2 = '.0 '
    else:
        return -1
    
    # Simulate and collect
    dssText.command = 'edit Fault.F1 bus2='+faultBuses+Ftype2+' bus1='+faultBuses+Ftype1+' r='+Fres+' enabled=true ontime=0.05'
    dssText.command = 'set mode=dynamic controlmode=time'
    dssText.command = 'set stepsize=0.0002s number=500'
    dssText.Command = 'solve'
    if(dss.Solution.Converged()):
        for ii in range(len(devLines)):
            data[ii][0] = faultBuses
            data[ii][1] = devNames[ii]
            Vmag,Imag,In,P,Q,V,I = getLineVI(dssCircuit,devLines[ii])
            data[ii][2] = Ftype.replace('.','_')
            # 3ph line
            if(len(I)==12):
                data[ii][3] = I[0]
                data[ii][4] = I[1]
                data[ii][5] = I[2]
                data[ii][6] = I[3]
                data[ii][7] = I[4]
                data[ii][8] = I[5]
                data[ii][9] = V[0]/(dev_BusV[ii])
                data[ii][10] = V[1]
                data[ii][11] = V[2]/(dev_BusV[ii])
                data[ii][12] = V[3]
                data[ii][13] = V[4]/(dev_BusV[ii])
                data[ii][14] = V[5]
                I0 = pdef(I[0],I[1])+pdef(I[2],I[3])+pdef(I[4],I[5])
                data[ii][15] = np.abs(I0)
                data[ii][16] = np.angle(I0)*(180/np.pi)
            # 2phase line
            elif(len(I)==8):
                data[ii][3] = I[0]
                data[ii][4] = I[1]
                data[ii][5] = I[2]
                data[ii][6] = I[3]
                data[ii][7] = 0#None
                data[ii][8] = 0#None
                data[ii][9] = V[0]/(dev_BusV[ii])
                data[ii][10] = V[1]
                data[ii][11] = V[2]/(dev_BusV[ii])
                data[ii][12] = V[3]
                data[ii][13] = 0#None
                data[ii][14] = 0#None
                data[ii][15] = max(I[0],I[2]) #np.abs(I0)
                data[ii][16] = 0 #np.angle(I0)*(180/np.pi)
            # single phase line
            elif(len(I)==4):
                data[ii][3] = I[0]
                data[ii][4] = I[1]
                data[ii][5] = 0#None
                data[ii][6] = 0#None
                data[ii][7] = 0#None
                data[ii][8] = 0#None
                data[ii][9] = V[0]/(dev_BusV[ii])
                data[ii][10] = V[1]
                data[ii][11] = 0#None
                data[ii][12] = 0#None
                data[ii][13] = 0#None
                data[ii][14] = 0#None
                data[ii][15] = I[0] #np.abs(I0)
                data[ii][16] = I[1] #np.angle(I0)*(180/np.pi)
    else:
        print('Failed to Converge for a '+ Ftype +' fault on bus'+faultBuses[ii]+'with res='+Fres)
        
    return pd.DataFrame(data)
