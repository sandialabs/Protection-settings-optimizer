

import networkx as nx
import numpy as np


def index_dict(lst, key, value):
    if(type(value)==str):
        for i, dic in enumerate(lst):
            if dic[key].lower() == value.lower():
                return i
        return None
    else:
        for i, dic in enumerate(lst):
            if dic[key] == value:
                return i
        return None

def genGraph(Lines,XFMRs,Buses,Relays,Recs,Fuses,Pvs,BESS,Gens,sW_Names,sW_status):
    nEdges = len(Lines)+len(XFMRs)
    Edges = [dict.fromkeys(['From','To','Name','isSwitch','isRelay','RelayName','isRecloser','RecloserName','isFuse','FuseName','color','width','Enabled']) for number in range(nEdges)]
    
    # Generate Graph()
    kk = 0
#    bus_names = [Bus['Name'] for Bus in Buses]
    for ii in range(len(Lines)):
        if (Lines[ii]['Enabled']):
            Edges[kk]['Enabled'] = 1
            Edges[kk]['Name'] = Lines[ii]['Name']
            Edges[kk]['From'] = Lines[ii]['Bus1'].split('.')[0]
            Edges[kk]['To'] = Lines[ii]['Bus2'].split('.')[0]
            Edges[kk]['isSwitch']  = Lines[ii]['isSwitch']=='True'
            Edges[kk]['numPhases'] = Lines[ii]['numPhases']
            Edges[kk]['color'] = 'k'
            Edges[kk]['width'] = 2
            kk = kk+1
            
    for ii in range(len(XFMRs)):
        if (XFMRs[ii]['Enabled']):
            Edges[kk]['Enabled'] = 1
            Edges[kk]['Name'] = XFMRs[ii]['Name']
            Edges[kk]['From'] = XFMRs[ii]['Bus1'].split('.')[0]
            Edges[kk]['To'] = XFMRs[ii]['Bus2'].split('.')[0]
            Edges[kk]['isSwitch']  = False
            Edges[kk]['numPhases']  = XFMRs[ii]['numPhases']
            Edges[kk]['color'] = 'k'
            Edges[kk]['width'] = 2
            kk = kk+1
            
    del Edges[kk:]
    
    for ii in range(len(sW_Names)):
        Edge_ind = index_dict(Edges,'Name',sW_Names[ii])
        if(Edge_ind == None):
            pass #print('Err')
        else:
            Edges[Edge_ind]['status'] = sW_status[ii]
    
    pos = dict.fromkeys([bus['Name'] for bus in Buses])
    G = nx.Graph()
    for ii in range(len(Buses)):
        G.add_node(Buses[ii]['Name'],Name=Buses[ii]['Name'])
        G.nodes[Buses[ii]['Name']]['nPhases'] = Buses[ii]['numPhases']
        G.nodes[Buses[ii]['Name']]['isPV'] = False
        G.nodes[Buses[ii]['Name']]['isBESS'] = False
        G.nodes[Buses[ii]['Name']]['isSource'] = False
        G.nodes[Buses[ii]['Name']]['pvName'] = False
        G.nodes[Buses[ii]['Name']]['bessName'] = False
        G.nodes[Buses[ii]['Name']]['sourceName'] = False
        G.nodes[Buses[ii]['Name']]['color'] = 'y'
        G.nodes[Buses[ii]['Name']]['size'] = 100
        
        pos[Buses[ii]['Name']] = np.array([Buses[ii]['X'],Buses[ii]['Y']])
        if(np.isnan(pos[Buses[ii]['Name']][0])):
            pos[Buses[ii]['Name']][0] = 0
        if(np.isnan(pos[Buses[ii]['Name']][1])):
            pos[Buses[ii]['Name']][1] = 1500        
    
    
    for ii in range(len(Edges)):        
        if(Edges[ii]['Enabled']==1):
            G.add_edge(Edges[ii]['From'],Edges[ii]['To'])
            nx.set_edge_attributes(G,{(Edges[ii]['From'],Edges[ii]['To']): Edges[ii]})
    
    # Extract realy locations 
    if(len(Relays) == 0):
        print('No Relays')
    else:
        for ii in range(len(Relays)):
            line_ind = index_dict(Edges,'Name',Relays[ii]['MonitoredObj'].split('.')[1])
            if(line_ind != None and Edges[line_ind]['Enabled']==1):
                if(Relays[ii]['Enabled']):
                    G[Relays[ii]['Bus1']][Relays[ii]['Bus2']]['isRelay'] = True
                    G[Relays[ii]['Bus1']][Relays[ii]['Bus2']]['RelayName'] = Relays[ii]['Name']
                    G[Relays[ii]['Bus1']][Relays[ii]['Bus2']]['color'] = 'g'
                    G[Relays[ii]['Bus1']][Relays[ii]['Bus2']]['width'] = 5
                else:
                    pass
                    #G[Relays[ii]['Bus1']][Relays[ii]['Bus2']] = False
                            
    # Extract Rec locations
    if(len(Recs) == 0):
        print('No Reclosers')
    else:
        for ii in range(len(Recs)):
            line_ind = index_dict(Edges,'Name',Recs[ii]['MonitoredObj'].split('.')[1])
            if(line_ind != None and Edges[line_ind]['Enabled']==1):
                if(Recs[ii]['Enabled']):
                    G[Recs[ii]['Bus1']][Recs[ii]['Bus2']]['isRecloser'] = True
                    G[Recs[ii]['Bus1']][Recs[ii]['Bus2']]['RecloserName'] = Recs[ii]['Name']
                    G[Recs[ii]['Bus1']][Recs[ii]['Bus2']]['color'] = 'b'
                    G[Recs[ii]['Bus1']][Recs[ii]['Bus2']]['width'] = 5
                else:
                    pass
                    #G[Recs[ii]['Bus1']][Recs[ii]['Bus2']] = False           
    
    # Extract Fuse Locations           
    if(len(Fuses) == 0):
        print('No Fuses')
    else:
        for ii in range(len(Fuses)):
            line_ind = index_dict(Edges,'Name',Fuses[ii]['MonitoredObj'].split('.')[1])
            if(line_ind != None and Edges[line_ind]['Enabled']==1):
                if(Fuses[ii]['Enabled']):
                    G[Fuses[ii]['Bus1']][Fuses[ii]['Bus2']]['isFuse'] = True
                    G[Fuses[ii]['Bus1']][Fuses[ii]['Bus2']]['FuseName'] = Fuses[ii]['Name']
                    G[Fuses[ii]['Bus1']][Fuses[ii]['Bus2']]['color'] = 'm'
                    G[Fuses[ii]['Bus1']][Fuses[ii]['Bus2']]['width'] = 5
                else:
                    pass
                    #G[Fuses[ii]['Bus1']][Fuses[ii]['Bus2']] = False                  
    
    # add Attrs at Nodes
    # Extract PVs
    if(len(Pvs)==0):
        print('No PV Systems active')
    else:
        for ii in range(len(Pvs)):
            if(Pvs[ii]['Enabled']):
                Bus_ind = Pvs[ii]['Bus'].split('.')[0]
                G.nodes[Bus_ind]['isPV'] = True
                G.nodes[Bus_ind]['pvName'] = Pvs[ii]['Name'] 
                G.nodes[Bus_ind]['color'] = 'g'
                G.nodes[Bus_ind]['size'] = 150
    
    # Extract BESS
    if(len(BESS)==0):
        print('No BESS active')
    else:
        for ii in range(len(BESS)):
            if(BESS[ii]['Enabled']):
                Bus_ind = BESS[ii]['Bus'].split('.')[0]
                G.nodes[Bus_ind]['isBESS'] = True
                G.nodes[Bus_ind]['bessName'] = BESS[ii]['Name'] 
                G.nodes[Bus_ind]['color'] = 'g'
                G.nodes[Bus_ind]['size'] = 150
    
    # Extract genes
    if(len(Gens)==0):
        print('No Gens active')
    else:
        for ii in range(len(Gens)):
            if(Gens[ii]['Enabled']):
                Bus_ind = Gens[ii]['Bus'].split('.')[0]
                G.nodes[Bus_ind]['isSource'] = True
                G.nodes[Bus_ind]['sourceName'] = Gens[ii]['Name'] 
                G.nodes[Bus_ind]['color'] = 'g'
                G.nodes[Bus_ind]['size'] = 150
                
    #Find Group
    Group = [list(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]
    for ii in range(len(Group)):
        for jj in range(len(Group[ii])):
            G.nodes[Group[ii][jj]]['Group'] = ii
    
    return (G,pos,Edges)