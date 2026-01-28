# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 14:56:07 2026

@author: tpatel
"""

import re
import numpy as np
import networkx as nx
from ..src.GenNxGraph import index_dict

def get(raven,target,default=0):
    return default if len(raven_find(raven,target)) == 0 else raven_find(raven,target)[0][1]

def has(raven,target):
    return len(raven_find(raven,target)) != 0

def raven_find(raven,target,path=""):
    """
    Input: raven JSON dictionary and a target member
    Output: list of all instances containing that target member with path to that instance 
    """
    found = []
    if isinstance(raven,dict):
        for key, val in list(raven.items()):
            #target found
            if key == target:
                found.append((path +"/"+str(key),val))
            #recursive step
            elif isinstance(val,dict):
                found += raven_find(val,target,path+"/"+str(key))
            elif isinstance(val,list):
                for item in val:
                    if isinstance(item,dict):
                        found += raven_find(item,target,path+"/"+str(key))
    return found

def get_switches_RAVENS(data,Buses,HCE_HC=0):
    Switches = [dict.fromkeys({'Name','Enabled','From','To','Type','Info','phasecode','kV','RAVENS_ID'}) for Sw in data['PowerSystemResource']['Equipment']['ConductingEquipment']['Switch']]
    ii = 0
    for Switch_name in data['PowerSystemResource']['Equipment']['ConductingEquipment']['Switch']:
        Switch = data['PowerSystemResource']['Equipment']['ConductingEquipment']['Switch'][Switch_name] 
        
        Switches[ii]['Name'] = Switch['IdentifiedObject.name'] 
        # check if Switch state is in the data if not assumen closed
        if('Switch.open' in Switch.keys()):
            SwitchOpen = Switch['Switch.open']
        else:
            if('Switch.normalOpen' in Switch.keys()):
                SwitchOpen = Switch['Switch.normalOpen']
            else:
                SwitchOpen = False
        Switches[ii]['Enabled'] =  Switch['Equipment.inService'] 
        Switches[ii]['State'] = not SwitchOpen
        From_str = [x['Terminal.ConnectivityNode'] for x in Switch['ConductingEquipment.Terminals'] if x['ACDCTerminal.sequenceNumber'] == 1][0]
        To_str = [x['Terminal.ConnectivityNode'] for x in Switch['ConductingEquipment.Terminals'] if x['ACDCTerminal.sequenceNumber'] == 2][0]
        Switches[ii]['From'] = re.search("ConnectivityNode::\'(.+)\'",From_str).group(1)
        Switches[ii]['To'] = re.search("ConnectivityNode::\'(.+)\'",To_str).group(1)
        Switches[ii]['Type'] = Switch['Ravens.cimObjectType']
        Switches[ii]['Info'] = re.search("SwitchInfo::\'(.+)\'",Switch['PowerSystemResource.AssetDatasheet']).group(1)
        Switches[ii]['phasecode'] = [x['Terminal.phases'] for x in Switch['ConductingEquipment.Terminals'] if x['ACDCTerminal.sequenceNumber'] == 1][0]
        Bus_idx = index_dict(Buses,'Name', Switches[ii]['From'])
        Switches[ii]['kV'] = Buses[Bus_idx]['kV']
        #Switches[ii]['RAVENS_ID'] = 
        ii = ii+1
        
        
    if(HCE_HC==1):
        Switches.append({'From': 'BattT1',
                          'To': 'bess1600',
                          'Type': 'Breaker',
                          'Info': 'DEFAULT',
                          'phasecode': 'PhaseCode.ABC',
                          'RAVENS_ID': None,
                          'kV': 0.332,
                          'Enabled': True,
                          'Name': 'bess1600',
                          'State': True})
    return Switches

def update_states_RAVENS(data,Switches,Lines,XFMRs,Pvs,BESS,Step_num):
    Statuses = data['AnalysisResult']['OptimalPowerFlow']['OperationsResult.Statuses']
    for status in Statuses:
        Equipment_Type = status['ArStatus.ConductingEquipment'].split('::')[0]
        # update Lines States Based on powerflow results                 
        if(Equipment_Type == 'ACLineSegment'):
            # Updatae line status
            step_idx = index_dict(status['AnalysisResultData.Curve']['AnalysisResultCurve.CurveDatas'],'ArCurveData.xvalue',Step_num)
            line_InService =  status['AnalysisResultData.Curve']['AnalysisResultCurve.CurveDatas'][step_idx]['ArCurveData.DataValues']['AvStatus.inService']
            line_Name  = status['ArStatus.ConductingEquipment'].split('::')[1]
            line_idx = index_dict(Lines,'Name',line_Name[1:-1])
            if(line_InService):
                Lines[line_idx]['Enabled']  = True
            else:
                Lines[line_idx]['Enabled']  = False
        # update BESS States Based on powerflow results                 
        elif(Equipment_Type == 'BatteryUnit'):
            step_idx = index_dict(status['AnalysisResultData.Curve']['AnalysisResultCurve.CurveDatas'],'ArCurveData.xvalue',Step_num)
            bess_InService =  status['AnalysisResultData.Curve']['AnalysisResultCurve.CurveDatas'][step_idx]['ArCurveData.DataValues']['AvStatus.inService']
            bess_Name  = status['ArStatus.ConductingEquipment'].split('::')[1]
            bess_idx = index_dict(BESS,'Name',bess_Name[1:-1])
            if(bess_InService):
                BESS[bess_idx]['Enabled']  = True
            else:
                BESS[bess_idx]['Enabled']  = False
                
        elif(Equipment_Type == 'EnergyConsumer'):
            pass
        # update Pvs States Based on powerflow results                 
        elif(Equipment_Type == 'PhotoVoltaicUnit'):
            step_idx = index_dict(status['AnalysisResultData.Curve']['AnalysisResultCurve.CurveDatas'],'ArCurveData.xvalue',Step_num)
            pv_InService =  status['AnalysisResultData.Curve']['AnalysisResultCurve.CurveDatas'][step_idx]['ArCurveData.DataValues']['AvStatus.inService']
            pv_Name  = status['ArStatus.ConductingEquipment'].split('::')[1]
            pv_idx = index_dict(Pvs,'Name',pv_Name[1:-1])
            if(pv_InService):
                Pvs[pv_idx]['Enabled']  = True
            else:
                Pvs[pv_idx]['Enabled']  = False
        # update Tx States Based on powerflow results                         
        elif(Equipment_Type == 'PowerTransformer'):
            step_idx = index_dict(status['AnalysisResultData.Curve']['AnalysisResultCurve.CurveDatas'],'ArCurveData.xvalue',Step_num)
            tx_InService =  status['AnalysisResultData.Curve']['AnalysisResultCurve.CurveDatas'][step_idx]['ArCurveData.DataValues']['AvStatus.inService']
            tx_Name  = status['ArStatus.ConductingEquipment'].split('::')[1]
            tx_idx = index_dict(XFMRs,'Name',tx_Name[1:-1])
            if(tx_InService):
                XFMRs[tx_idx]['Enabled']  = True
            else:
                XFMRs[tx_idx]['Enabled']  = False
        # update Switch States Based on powerflow results                 
        elif(Equipment_Type == 'Switch'):
            step_idx = index_dict(status['AnalysisResultData.Curve']['AnalysisResultCurve.CurveDatas'],'ArCurveData.xvalue',Step_num)
            switch_InService =  status['AnalysisResultData.Curve']['AnalysisResultCurve.CurveDatas'][step_idx]['ArCurveData.DataValues']['AvStatus.inService']
            switch_Name  = status['ArStatus.ConductingEquipment'].split('::')[1]
            Switch_stat = [x for x in data['AnalysisResult']['OptimalPowerFlow']['OperationsResult.Switches'] if x['ArSwitch.Switch'].split('::')[1] == switch_Name][0]
            switch_Open = Switch_stat['AnalysisResultData.Curve']['AnalysisResultCurve.CurveDatas'][step_idx]['ArCurveData.DataValues']['AvSwitch.open']
            Switch_idx = index_dict(Switches,'Name',switch_Name[1:-1])
            
            if(switch_InService and not switch_Open):
                Switches[Switch_idx]['Enabled'] = True
                Switches[Switch_idx]['State'] = True
            elif(switch_InService and switch_Open):
                Switches[Switch_idx]['Enabled'] = True
                Switches[Switch_idx]['State'] = False
            elif(not switch_InService):
                Switches[Switch_idx]['Enabled'] = False
                Switches[Switch_idx]['State'] = False
            else:
                print('switch state for {} could not be determined'.format(switch_Name))
            
        else:
            print('Error: could not update Equipment status')
    return Switches,Lines,XFMRs,Pvs,BESS

def get_Lines_RAVENS(data,Buses):
    Lines = [dict.fromkeys({'Name', 'Bus1', 'Bus2', 'Enabled', 'isSwitch', 'numPhases', 'kV', 'Length', 'Rpu', 'Xpu'}) # 'R1pu', 'R0pu', 'X1pu', 'X0pu'
             for Line in data['PowerSystemResource']['Equipment']['ConductingEquipment']['Conductor']['ACLineSegment'] ]
    ii = 0
    for Line_name in data['PowerSystemResource']['Equipment']['ConductingEquipment']['Conductor']['ACLineSegment'] :
        Line = data['PowerSystemResource']['Equipment']['ConductingEquipment']['Conductor']['ACLineSegment'][Line_name]
        Lines[ii]['Name'] = Line['IdentifiedObject.name']
        From_str = [x['Terminal.ConnectivityNode'] for x in Line['ConductingEquipment.Terminals'] if x['ACDCTerminal.sequenceNumber'] == 1][0]
        To_str = [x['Terminal.ConnectivityNode'] for x in Line['ConductingEquipment.Terminals'] if x['ACDCTerminal.sequenceNumber'] == 2][0]
        Lines[ii]['Bus1'] = re.search("ConnectivityNode::\'(.+)\'",From_str).group(1)
        Lines[ii]['Bus2'] = re.search("ConnectivityNode::\'(.+)\'",To_str).group(1)
        Lines[ii]['Enabled'] = Line['Equipment.inService']
        Lines[ii]['isSwitch'] = False
        phasecode = [x['Terminal.phases'] for x in Line['ConductingEquipment.Terminals'] if x['ACDCTerminal.sequenceNumber'] == 1][0]
        Lines[ii]['numPhases'] = len(re.search("PhaseCode.(.+)",phasecode).group(1))
        if('ConductingEquipment.BaseVoltage' in Line.keys()):
            Lines[ii]['kV'] = Line['ConductingEquipment.BaseVoltage']
        else:
            Bus_idx = index_dict(Buses,'Name', Lines[ii]['Bus1'])
            Lines[ii]['kV'] =  Buses[Bus_idx]['kV'];
        Lines[ii]['Length'] = Line['Conductor.length']
        Lines[ii]['Rpu'] = 1
        Lines[ii]['Xpu'] = 1
        Lines[ii]['phasecode'] = re.search("PhaseCode.(.+)",phasecode).group(1)
        ii = ii+1
    return Lines

def get_XFMRs_RAVENS(data,HCE_HC=0):
    XFMRs = [dict.fromkeys({'Name', 'Bus1', 'Bus2', 'Enabled', 'numPhases', 'MVA', 'Rpu', 'Xpu'}) 
             for Tx in data['PowerSystemResource']['Equipment']['ConductingEquipment']['PowerTransformer']]
    ii = 0
    for Tx_name in data['PowerSystemResource']['Equipment']['ConductingEquipment']['PowerTransformer']:
        Tx = data['PowerSystemResource']['Equipment']['ConductingEquipment']['PowerTransformer'][Tx_name]
        XFMRs[ii]['Name'] = Tx['IdentifiedObject.name']
        From_str = [x['Terminal.ConnectivityNode'] for x in Tx['ConductingEquipment.Terminals'] if x['ACDCTerminal.sequenceNumber'] == 1][0]
        To_str = [x['Terminal.ConnectivityNode'] for x in Tx['ConductingEquipment.Terminals'] if x['ACDCTerminal.sequenceNumber'] == 2][0]
        XFMRs[ii]['Bus1'] = re.search("ConnectivityNode::\'(.+)\'",From_str).group(1)
        XFMRs[ii]['Bus2'] = re.search("ConnectivityNode::\'(.+)\'",To_str).group(1)
        XFMRs[ii]['Enabled'] = Tx['Equipment.inService']
        phasecode = [x['Terminal.phases'] for x in Tx['ConductingEquipment.Terminals'] if x['ACDCTerminal.sequenceNumber'] == 1][0]
        XFMRs[ii]['numPhases'] = len(re.search("PhaseCode.(.+)",phasecode).group(1))
        XFMRs[ii]['MVA'] = 100
        XFMRs[ii]['Rpu'] = 1
        XFMRs[ii]['Xpu'] = 1
        XFMRs[ii]['phasecode'] = re.search("PhaseCode.(.+)",phasecode).group(1)
        ii = ii + 1
        
    if(HCE_HC==1):
        XFMRs[index_dict(XFMRs,'Name','BattT1')]['Bus2'] = 'bess1600'
    return XFMRs

def get_Buses_RAVENS(data,BaseVoltages=[],HCE_HC=0): 
    nBuses = len(set([x['ArVoltage.ConnectivityNode'] for x in  data['AnalysisResult']['OptimalPowerFlow']['OperationsResult.Voltages']]))
    Buses = [dict.fromkeys({'Name', 'nodes', 'numPhases', 'X', 'Y', 'kV'}) for x in range(nBuses) ]
    ii = 0
    for Bus_phase in data['AnalysisResult']['OptimalPowerFlow']['OperationsResult.Voltages']:
        # check if bus allreay exists 
        Bus_Name = re.search("ConnectivityNode::\'(.+)\'",Bus_phase['ArVoltage.ConnectivityNode']).group(1)
        Bus_in_list = Bus_Name in [Bus['Name'] for Bus in Buses]
        if(Bus_in_list):
            # Bus allready in list
            Bus_idx = None
            Bus_idx = index_dict(Buses,'Name', Bus_Name)
            New_node = ord((re.search("SinglePhaseKind.(.+)",Bus_phase['AnalysisResultData.phase']).group(1)).upper())-64
            if(New_node in Buses[Bus_idx]['nodes']):
                pass # Not sure why there are reapeated values 
                #print("repeted entry for Bus {Bus} phase {ph}".format(Bus = Bus_Name, ph = str(New_node)) )
            else:
                Buses[Bus_idx]['numPhases'] += 1
                Buses[Bus_idx]['nodes'].append(New_node)
                Buses[Bus_idx]['nodes'].sort()
        else:
            # bus not in list
            Buses[ii]['Name'] = Bus_Name
            Buses[ii]['numPhases'] = 1 
            Buses[ii]['nodes'] = [ord((re.search("SinglePhaseKind.(.+)",Bus_phase['AnalysisResultData.phase']).group(1)).upper())-64]
            # estimate Bus nominal voltage 
            Buses[ii]['kV'] = 0 
            if(has(Bus_phase,'ArCurveData.DataValues')):
                Varr = np.array([x['ArCurveData.DataValues']['AvVoltage.v'] for x in Bus_phase['AnalysisResultData.Curve']['AnalysisResultCurve.CurveDatas']])
                Vmean = Varr[np.nonzero(Varr)].mean()
                if(np.isnan(Vmean) ):
                    Buses[ii]['kV'] = 0; 
                else:
                    Buses[ii]['kV'] = float(Vmean) / 1e3
                    #Buses[ii]['kV'] = min(BaseVoltages, key=lambda x:abs(x-Vmean)) / 1e3
            Buses[ii]['X'] = np.nan
            Buses[ii]['Y'] = np.nan
            ii += 1 
    if(HCE_HC == 1):
        Buses.append({'numPhases': 3,'nodes': [1, 2, 3],'kV': 0.332,'Name': 'bess1600','X':np.nan,'Y':np.nan})
    return Buses

def get_BESS_RAVENS(data):
    # Battries (teated as sources)
    batt_data = [x for x in data['AnalysisResult']['OptimalPowerFlow']['OperationsResult.PowerFlows'] if 'BatteryUnit::' in x['ArPowerFlow.ConductingEquipment']]
    nBatts = len(set([x['ArPowerFlow.ConductingEquipment'] for x in data['AnalysisResult']['OptimalPowerFlow']['OperationsResult.PowerFlows'] if 'BatteryUnit::' in x['ArPowerFlow.ConductingEquipment']]))
    BESS = [dict.fromkeys({'Name', 'Bus', 'Enabled', 'numPhases'}) for x in range(nBatts)]
    ii = 0
    for Bat_phase in batt_data:
        Batt_Name = re.search("BatteryUnit::\'(.+)\'",Bat_phase['ArPowerFlow.ConductingEquipment']).group(1)
        Batt_in_list = Batt_Name in [BES['Name'] for BES in BESS]
        if(Batt_in_list):
            Bus_idx = None
            Bus_idx = index_dict(BESS,'Name', Batt_Name)
            BESS[Bus_idx]['numPhases'] += 1
        else:
            Batt_sys_data = get(data['PowerSystemResource']['Equipment']['ConductingEquipment']['EnergyConnection'],Batt_Name,default=-1)
            if(Batt_sys_data == -1):
                raise KeyError("Unique entry for BESS {BESS_Name} Not found in PowerSystemResource".format(BESS_Name = Batt_Name))
            else:
                BESS[ii]['Name'] = Batt_Name
                Batt_terminals =get(Batt_sys_data,'Terminal.ConnectivityNode')
                BESS[ii]['Bus'] = re.search("ConnectivityNode::\'(.+)\'",Batt_terminals).group(1)
                BESS[ii]['Enabled'] = Batt_sys_data['Equipment.inService']
                BESS[ii]['numPhases'] = 1
                ii += 1
    return BESS

def get_Pvs_RAVENS(data):
    # PVS (teated as sources)
    PV_data = [x for x in data['AnalysisResult']['OptimalPowerFlow']['OperationsResult.PowerFlows'] if 'PhotoVoltaicUnit::' in x['ArPowerFlow.ConductingEquipment']]
    nPVs = len(set([x['ArPowerFlow.ConductingEquipment'] for x in data['AnalysisResult']['OptimalPowerFlow']['OperationsResult.PowerFlows'] if 'PhotoVoltaicUnit::' in x['ArPowerFlow.ConductingEquipment']]))
    Pvs = [dict.fromkeys({'Name', 'Bus', 'Enabled', 'numPhases'}) for x in range(nPVs)]
    ii = 0
    for PV_phase in PV_data:
        PV_Name = re.search("PhotoVoltaicUnit::\'(.+)\'",PV_phase['ArPowerFlow.ConductingEquipment']).group(1)
        PV_in_list = PV_Name in [Pv['Name'] for Pv in Pvs]
        if(PV_in_list):
            Bus_idx = None
            Bus_idx = index_dict(Pvs,'Name', PV_Name)
            Pvs[Bus_idx]['numPhases'] += 1
        else:
            PV_sys_data = get(data['PowerSystemResource']['Equipment']['ConductingEquipment']['EnergyConnection'],PV_Name,default=-1)
            if(PV_sys_data == -1):
                raise KeyError("Unique entry for BESS {PV_Name} Not found in PowerSystemResource".format(PV_Name = PV_Name))
            else:
                Pvs[ii]['Name'] = PV_Name
                PV_terminals =get(PV_sys_data,'Terminal.ConnectivityNode')
                Pvs[ii]['Bus'] = re.search("ConnectivityNode::\'(.+)\'",PV_terminals).group(1)
                Pvs[ii]['Enabled'] = PV_sys_data['Equipment.inService']
                Pvs[ii]['numPhases'] = 1
                ii += 1    
    return Pvs

def get_Relays_RAVENS(data,Switches):
    nRelays = len([x['Type'] for x in Switches if (x['Type'] == 'Breaker' or x['Type'] == 'ProtectedSwitch')])
    if(nRelays == 0):
        Relays = []
        return Relays
    
    Relays = [dict.fromkeys( {'Name','Bus1','Bus2','MonitoredObj','Enabled'}) for x in range(nRelays)]
    ii=0
    for switch in Switches:
        if(switch['Type'] == 'Breaker' or switch['Type'] == 'ProtectedSwitch'):
            Relays[ii]['Name'] = switch['Name']
            Relays[ii]['Bus1'] = switch['From']
            Relays[ii]['Bus2'] = switch['To']
            Relays[ii]['Enabled'] = switch['State']
            Relays[ii]['MonitoredObj'] = "line." + switch['Name']
            ii += 1
    return Relays

def get_Recs_RAVENS(data,Switches):
    nRecs = len([x['Type'] for x in Switches if (x['Type'] == 'Recloser')])
    if(nRecs == 0):
        Recs = []
        return Recs
    
    Recs = [dict.fromkeys( {'Name','Bus1','Bus2','MonitoredObj','Enabled'}) for x in range(nRecs)]
    ii=0
    for switch in Switches:
        if(switch['Type'] == 'Recloser'):
            Recs[ii]['Name'] = switch['Name']
            Recs[ii]['Bus1'] = switch['From']
            Recs[ii]['Bus2'] = switch['To']
            Recs[ii]['Enabled'] = switch['State']
            Recs[ii]['MonitoredObj'] = "MonitoredObj." + switch['Name']
            ii += 1
    return Recs

def get_Fuses_RAVENS(data,Switches):
    nFuses = len([x['Type'] for x in Switches if (x['Type'] == 'Fuse')])
    if(nFuses == 0):
        Fuses = []
        return Fuses
    Fuses = [dict.fromkeys( {'Name', 'Bus1', 'Bus2', 'MonitoredObj', 'Enabled', 'Type', 'Rating', 'OT'}) for x in range(nFuses)]
    ii=0
    for switch in Switches:
        if(switch['Type'] == 'Fuse'):
            Fuses[ii]['Name'] = switch['Name']
            Fuses[ii]['Bus1'] = switch['From']
            Fuses[ii]['Bus2'] = switch['To']
            Fuses[ii]['Enabled'] = switch['State']
            Fuses[ii]['MonitoredObj'] = "MonitoredObj." + switch['Name']
            Fuses[ii]['Type'] = switch['Info'].split('-')[1]
            Fuses[ii]['Rating'] = switch['Info'].split('-')[0]
            # Need to implinet Fuse curves in some way later, Ignore for now
            #  Fuses[ii]['OT'] = {'xTC': None, 'yTC':None,'xMM':None,'yMM':None}
            ii += 1
    return Fuses
    
def add_Switches_as_lines(Lines,Switches,Buses):
    for switch in Switches:
        B1 = switch['From']
        B2 = switch['To']
        Line_in_lines = [x for x in Lines if (x['Bus1'] == B1 and x['Bus2'] == B2) or (x['Bus1'] == B2 and x['Bus2'] == B1)]   
        if not Line_in_lines:
            NewLine = {}
            NewLine['Name'] = switch['Name']
            NewLine['Bus1'] = B1
            NewLine['Bus2'] = B2
            NewLine['Enabled'] = switch['Enabled'] 
            NewLine['isSwitch'] = True
            #phasecode = [x['Terminal.phases'] for x in Line['ConductingEquipment.Terminals'] if x['ACDCTerminal.sequenceNumber'] == 1][0]
            NewLine['numPhases'] = Buses[index_dict(Buses,'Name',B1)]['numPhases']
            NewLine['kV'] =  Buses[index_dict(Buses,'Name',B1)]['kV']
            NewLine['Length'] = 0
            NewLine['Rpu'] = 1e-6
            NewLine['Xpu'] = 1e-6
            NewLine['phasecode'] = Buses[index_dict(Buses,'Name',B1)]['nodes']
            Lines.append(NewLine)
        else:
            pass

def plot_NX_draw_kamada_kawai(G):
    Ecolors = [G[u][v]['color'] for u,v in G.edges]
    Ewidth = [G[u][v]['width'] for u,v in G.edges]
    Ncolors = [G.nodes[u]['color'] for u in G.nodes]
    Nsize = [G.nodes[u]['size'] for u in G.nodes]
    nx.draw_kamada_kawai(  G,with_labels=1,edge_color=Ecolors   ,width = Ewidth   ,node_size=Nsize   ,node_color = Ncolors) 

def update_Protection_Devices_RAVENS(Relays,Recs,Fuses):
    Relays = [x for x in Relays if x['Enabled']]
    Recs = [x for x in Recs if x['Enabled']]
    Fuses = [x for x in Fuses if x['Enabled']]
    return Relays,Recs,Fuses         

def update_DG_sources_RAVENS(Pvs,BESS,Gens):
    Pvs = [x for x in Pvs if x['Enabled']]
    BESS = [x for x in BESS if x['Enabled']]
    Gens = [x for x in Gens if x['Enabled']]
    return Pvs,BESS,Gens

def update_system_dicts_RAVENS(Buses,Lines,XFMRs):
    Lines = [x for x in Lines if x['Enabled']]
    XFMRs = [x for x in XFMRs if x['Enabled']]
    active_buses_set = set( [x['Bus1'] for x in Lines] + [x['Bus2'] for x in Lines] + [x['Bus1'] for x in XFMRs] + [x['Bus2'] for x in XFMRs] )
    Buses = [x for x in Buses if x['Name'] in active_buses_set]
    return Lines, XFMRs, Buses

def get_fault_data_RAVENS(data,Step_num):
    FData = []
    Fault_Data = [data['AnalysisResult'][x] for x in data['AnalysisResult'] if data['AnalysisResult'][x]['Ravens.cimObjectType'] == 'FaultStudyResult']
    for fault in Fault_Data:
        fault_ID = fault['FaultStudyResult.Fault']
        fault_Name = fault_ID.split('::')[-1][1:-1]
        fault_Info = data['Fault'][fault_Name]
        fault_loc = fault_Info['Fault.FaultyEquipment'].split('::')[-1][1:-1]
        fault_phase = fault_Info['Fault.phases'].split('.')[-1]
        fault_res = fault_Info['Fault.impedance']['FaultImpedance.rLineToLine'] + fault_Info['Fault.impedance']['FaultImpedance.rGround']    
        fault_type = fault_Info['Fault.kind'].split('.')[-1]
        if(fault_type == 'threePhase'):
            f_type = 'TPH_R' + str(fault_res).replace('.','_')
        elif(fault_type == 'lineToGround'):
            f_type = 'SLG_A_R' + str(fault_res).replace('.','_')
        elif(fault_type == 'lineToLine'):
            f_type = 'BC_R' + str(fault_res).replace('.','_')
        else:
            print('Error fault type not found\n')
            
        for fault_phase in fault['OperationsResult.Voltages']:
            f_pro_device_ = fault_phase['ArVoltage.ConnectivityNode']
            f_pro_device = f_pro_device_.split('::')[-1][1:-1]
            f_phase_ = fault_phase['AnalysisResultData.phase']
            f_phase = fault_phase['AnalysisResultData.phase'].split('.')[-1]
            
            step_idx = index_dict(fault_phase['AnalysisResultData.Curve']['AnalysisResultCurve.CurveDatas'],'ArCurveData.xvalue',Step_num)
            f_Vang = fault_phase['AnalysisResultData.Curve']['AnalysisResultCurve.CurveDatas'][step_idx]['ArCurveData.DataValues']['AvVoltage.angle']
            f_Vmag = fault_phase['AnalysisResultData.Curve']['AnalysisResultCurve.CurveDatas'][step_idx]['ArCurveData.DataValues']['AvVoltage.v']

            curretns = [x for x in fault['OperationsResult.CurrentFlows'] if x['AnalysisResultData.phase'] == f_phase_ and x['ArCurrent.ConnectivityNode'] == f_pro_device_][0]
            f_Imag = curretns['AnalysisResultData.Curve']['AnalysisResultCurve.CurveDatas'][step_idx]['ArCurveData.DataValues']['AvCurrent.i']
            f_Iang = curretns['AnalysisResultData.Curve']['AnalysisResultCurve.CurveDatas'][step_idx]['ArCurveData.DataValues']['AvCurrent.angle']
            
            F_idx = next((x for x,sublist in enumerate(FData) if (sublist[0] == fault_loc and sublist[1] == f_pro_device and sublist[2] == f_type)),None)
            if F_idx is None: 
                FData.append([fault_loc,f_pro_device,f_type,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
                if(f_phase == 'A'):
                    FData[-1][3] = f_Imag
                    FData[-1][4] = f_Iang
                    
                    FData[-1][9] = f_Vmag
                    FData[-1][10] = f_Vang
                if(f_phase == 'B'):
                    FData[-1][5] = f_Imag
                    FData[-1][6] = f_Iang
                    
                    FData[-1][11] = f_Vmag
                    FData[-1][12] = f_Vang
                if(f_phase == 'C'):
                    FData[-1][7] = f_Imag
                    FData[-1][8] = f_Iang
                    
                    FData[-1][13] = f_Vmag
                    FData[-1][14] = f_Vang
            else:
                if(f_phase == 'A'):
                    FData[F_idx][3] = f_Imag
                    FData[F_idx][4] = f_Iang
                    
                    FData[F_idx][9] = f_Vmag
                    FData[F_idx][10] = f_Vang
                if(f_phase == 'B'):
                    FData[F_idx][5] = f_Imag
                    FData[F_idx][6] = f_Iang
                    
                    FData[F_idx][11] = f_Vmag
                    FData[F_idx][12] = f_Vang
                if(f_phase == 'C'):
                    FData[F_idx][7] = f_Imag
                    FData[F_idx][8] = f_Iang
                    
                    FData[F_idx][13] = f_Vmag
                    FData[F_idx][14] = f_Vang

    for fdata in FData:
        Ia = (fdata[3] * np.cos(fdata[4]*(np.pi/180)) ) + ( fdata[3] * np.sin(fdata[4]*(np.pi/180)) )*1j
        Ib = (fdata[5] * np.cos(fdata[6]*(np.pi/180)) ) + ( fdata[5] * np.sin(fdata[6]*(np.pi/180)) )*1j
        Ic = (fdata[7] * np.cos(fdata[8]*(np.pi/180)) ) + ( fdata[7] * np.sin(fdata[8]*(np.pi/180)) )*1j
        
        fdata[15] = float(abs(Ia + Ib + Ic))
        fdata[16] = float(np.angle(Ia + Ib + Ic,1))
        
    return FData

def get_pf_data_phase(dev_Voltages,dev_PQ,dev_bus,devLine,Step_num,PhaseCode):
    
    pro_dev_bus_Va_ = [x for x in dev_Voltages if x['ArVoltage.ConnectivityNode'].split('::')[-1][1:-1] == dev_bus and x['AnalysisResultData.phase'].split('.')[-1] == PhaseCode]
    pro_dev_bus_PQa_ = [x for x in dev_PQ if x['ArPowerFlow.ConductingEquipment'].split('::')[-1][1:-1] == devLine and x['AnalysisResultData.phase'].split('.')[-1] == PhaseCode]
    
    if(len(pro_dev_bus_Va_) == 1 and len(pro_dev_bus_PQa_) == 1):
        pro_dev_bus_Va = pro_dev_bus_Va_[0]
        pro_dev_bus_PQa = pro_dev_bus_PQa_[0]
        
        step_idx = index_dict(pro_dev_bus_Va['AnalysisResultData.Curve']['AnalysisResultCurve.CurveDatas'],'ArCurveData.xvalue',Step_num)
        Va_pf = pro_dev_bus_Va['AnalysisResultData.Curve']['AnalysisResultCurve.CurveDatas'][step_idx]['ArCurveData.DataValues']['AvVoltage.v']
        Pa_pf = pro_dev_bus_PQa['AnalysisResultData.Curve']['AnalysisResultCurve.CurveDatas'][step_idx]['ArCurveData.DataValues']['AvPowerFlow.p']
        Qa_pf = pro_dev_bus_PQa['AnalysisResultData.Curve']['AnalysisResultCurve.CurveDatas'][step_idx]['ArCurveData.DataValues']['AvPowerFlow.q']
        PQa_pf = Pa_pf + (1j*Qa_pf)
        
    else:
        print('Protection device pre-fault voltage not found for ' + devLine)
    
    return Va_pf,PQa_pf

def abc2012(Vabc):
    a = 1*np.exp(120*(np.pi/180)*1j)
    #A = np.array([[1,1,1],[1,a**2,a],[1,a,a**2]])
    #Ainv = np.linalg.inv(A)
    
    V0 = (1/3) * (Vabc[0]+Vabc[1]+Vabc[2])
    V1 = (1/3) * (Vabc[0] + a*Vabc[1] + a**2*Vabc[2]) ;
    V2 = (1/3) * (Vabc[0] + a**2*Vabc[1] + a*Vabc[2]) ;
    return [V0,V1,V2]

def get_Device_Data_RAVENS(data,SysInfo,Step_num):
    
    Buses = SysInfo['Buses']
    devLines = [x['MonitoredObj'].split('line.')[1] for x in SysInfo['Relays']]+[x['MonitoredObj'].split('line.')[1] for x in SysInfo['Recs']]
    devNames = [x['Name'] for x in SysInfo['Relays']]+[x['Name'] for x in SysInfo['Recs']]
    dev_BusV = [Buses[index_dict(Buses,'Name',x['Bus1'])]['kV']*1e3 for x in SysInfo['Relays'] ]+[Buses[index_dict(Buses,'Name',x['Bus1'])]['kV']*1e3 for x in SysInfo['Recs']]
    dev_buses = [x['Bus1'] for x in SysInfo['Relays']] +[x['Bus1'] for x in SysInfo['Recs']]
    
    Device_Data_CSV = [dict.fromkeys(['RelayName','Ia_mag','Ia_ang','Va_mag','Va_ang','P','Q','switchState','PVs','In']) for number in range(len(devNames))]
    
    dev_Voltages = data['AnalysisResult']['OptimalPowerFlow']['OperationsResult.Voltages']
    dev_PQ = data['AnalysisResult']['OptimalPowerFlow']['OperationsResult.PowerFlows']
    
    for ii in range(len(devNames)):
        Va_pf,PQa_pf = get_pf_data_phase(dev_Voltages,dev_PQ,dev_buses[ii],devLines[ii],Step_num,'A')
        Vb_pf,PQb_pf = get_pf_data_phase(dev_Voltages,dev_PQ,dev_buses[ii],devLines[ii],Step_num,'B')
        Vc_pf,PQc_pf = get_pf_data_phase(dev_Voltages,dev_PQ,dev_buses[ii],devLines[ii],Step_num,'C')
    
        Va = (Va_pf * np.exp(1j*0*np.pi/180))
        Vb = (Vb_pf * np.exp(1j*-120*np.pi/180))
        Vc = (Vc_pf * np.exp(1j*-120*np.pi/180))
    
        Ia_pf = np.conj (PQa_pf/Va)
        Ib_pf = np.conj (PQb_pf/Vb)
        Ic_pf = np.conj (PQc_pf/Vc)
        
        Device_Data_CSV[ii]['RelayName'] = devNames[ii]
        Device_Data_CSV[ii]['Va_mag'] = float(min(Va_pf,Vb_pf,Vc_pf))
        Device_Data_CSV[ii]['Va_ang'] = 0
        
        Imax_ind = [abs(Ia_pf),abs(Ib_pf),abs(Ic_pf)].index(max([abs(Ia_pf),abs(Ib_pf),abs(Ic_pf)]))
        
        Device_Data_CSV[ii]['Ia_mag'] = float([abs(Ia_pf),abs(Ib_pf),abs(Ic_pf)][Imax_ind])
        Device_Data_CSV[ii]['Ia_ang'] = float([np.angle(Ia_pf,1),np.angle(Ib_pf,1),np.angle(Ic_pf,1)][Imax_ind])
        
        Device_Data_CSV[ii]['Vabc'] = [Va,Vb,Vc]
        Device_Data_CSV[ii]['Iabc'] = [Ia_pf,Ib_pf,Ic_pf]
        Device_Data_CSV[ii]['V012'] = abc2012(Device_Data_CSV[ii]['Vabc'])
        Device_Data_CSV[ii]['I012'] = abc2012(Device_Data_CSV[ii]['Iabc'])
        Device_Data_CSV[ii]['Z012'] = np.divide(Device_Data_CSV[ii]['V012'],Device_Data_CSV[ii]['I012'])
        Device_Data_CSV[ii]['In'] = abs(Ia_pf+Ib_pf+Ic_pf)
        
        Sa = Va * np.conj( Ia_pf )
        Sb = Vb * np.conj( Ib_pf )
        Sc = Vc * np.conj( Ic_pf )
        
        Device_Data_CSV[ii]['P'] = float(np.real(Sa+Sb+Sc))
        Device_Data_CSV[ii]['Q'] = float(np.imag(Sa+Sb+Sc))
        Device_Data_CSV[ii]['PVs'] = 1
        Device_Data_CSV[ii]['switchState'] = True
    return Device_Data_CSV

def get_system_Data_RAVENS(data,Step_num,include_Fuses = False,HCE_HC=0):
    
    Buses = get_Buses_RAVENS(data,HCE_HC=HCE_HC)
    Switches = get_switches_RAVENS(data,Buses,HCE_HC=HCE_HC) 
    Lines = get_Lines_RAVENS(data,Buses)
    XFMRs = get_XFMRs_RAVENS(data,HCE_HC=HCE_HC)
    
    # convert Generation data
    Pvs = get_Pvs_RAVENS(data)
    BESS = get_BESS_RAVENS(data)
    Gens = []
    Switches,Lines,XFMRs,Pvs,BESS = update_states_RAVENS(data,Switches,Lines,XFMRs,Pvs,BESS,Step_num)
    add_Switches_as_lines(Lines,Switches,Buses)
    
    # convert protection data
    Relays = get_Relays_RAVENS(data,Switches)
    Recs = get_Recs_RAVENS(data,Switches)
    Fuses = get_Fuses_RAVENS(data,Switches)
    # Trim dicts 
    Pvs,BESS,Gens = update_DG_sources_RAVENS(Pvs,BESS,Gens)
    Relays,Recs,Fuses = update_Protection_Devices_RAVENS(Relays,Recs,Fuses)
    Lines, XFMRs, Buses = update_system_dicts_RAVENS(Buses,Lines,XFMRs)
    
    
    SysInfo= {}   
    SysInfo['Relays'] = Relays
    SysInfo['Recs'] = Recs
    if(include_Fuses):
        SysInfo['Fuses'] = Fuses
    else:
        SysInfo['Fuses'] = [] 
    SysInfo['Lines'] = Lines
    SysInfo['XFMRs'] = XFMRs
    SysInfo['Buses'] = Buses
    SysInfo['Pvs'] = Pvs
    SysInfo['BESS'] = BESS
    SysInfo['Gens'] = Gens

    return SysInfo
