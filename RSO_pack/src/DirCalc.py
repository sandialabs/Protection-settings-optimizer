# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 13:18:12 2022

@author: maste
"""
import numpy as np
import networkx as nx
from .Read_CSV_Functions import pdef
from .GenNxGraph import index_dict

def c2p(C):
    Mag = np.abs(C)
    Ang = np.angle(C)*(180/np.pi)
    return Mag,Ang

def abc2012(Vabc):
    a = 1*np.exp(120*(np.pi/180)*1j)
    #A = np.array([[1,1,1],[1,a**2,a],[1,a,a**2]])
    #Ainv = np.linalg.inv(A)
    
    V0 = (1/3)*(Vabc[0]+Vabc[1]+Vabc[2])
    V1 = (1/3) * (Vabc[0] + a*Vabc[1] + a**2*Vabc[2]) ;
    V2 = (1/3) * (Vabc[0] + a**2*Vabc[1] + a*Vabc[2]) ;
    return [V0,V1,V2]

def calc_Z2(Z1ANG,BBB,a2):
    
    if('Vabc' not in BBB.keys()):
        VaF1 = pdef(BBB['Va_mag'],BBB['Va_ang'])*(4.16e3/np.sqrt(3))
        VbF1 = pdef(BBB['Vb_mag'],BBB['Vb_ang'])*(4.16e3/np.sqrt(3))
        VcF1 = pdef(BBB['Vc_mag'],BBB['Vc_ang'])*(4.16e3/np.sqrt(3))
        
        IaF1 = pdef(BBB['Ia_mag'],BBB['Ia_ang'])
        IbF1 = pdef(BBB['Ib_mag'],BBB['Ib_ang'])
        IcF1 = pdef(BBB['Ic_mag'],BBB['Ic_ang'])
        
        BBB['Vabc'] = [VaF1,VbF1,VcF1]
        BBB['Iabc'] = [IaF1,IbF1,IcF1]
    
    V012_F1 = abc2012(BBB['Vabc'])
    I012_F1 = abc2012(BBB['Iabc'])
    
    if(np.divide(abs(I012_F1[2]),abs(I012_F1[1]))>a2):
        Z2 = np.real( V012_F1[2]*( I012_F1[2]* np.conj(pdef(1,Z1ANG)) ))/ (np.abs(I012_F1[2])**2)
    else:
        Z2 = np.nan
    return Z2

def calcZ1Z0ANG(Rname,Device_Data_CSV,Fault_Data_CSV,G,B1,B2,Show_Warnings):
    Rind = index_dict(Device_Data_CSV,'RelayName',Rname)
    #RRind = index_dict(Relays,'Name',Rname)
    #B1 = Relays[RRind]['Bus1']
    #B2 = Relays[RRind]['Bus2']
    
    # % Record Z2 Data
    Z12T = 0.1
    Z1ANG = np.angle(Device_Data_CSV[Rind]['Z012'][1])*(180/np.pi)
    Z1MAG = abs(Device_Data_CSV[Rind]['Z012'][1])
    Z0MAG = abs(Device_Data_CSV[Rind]['Z012'][0])
    RFaults = [x for x in Fault_Data_CSV if x['Relay'].lower() == Rname.lower()]
    for F in RFaults:
        # if I1 exists (ie has none 0 current)
        if(abs(F['I012'][1]) > 0.001 or abs(F['I012'][0])>0.001):
            Z1F =  F['V012'][1]/F['I012'][1]
            Z0F =  F['V012'][0]/F['I012'][0]
            F['Z1'] = np.angle(Z1F,True) % 360
            F['Z0'] = np.angle(Z0F,True) % 360
            
            # Use Z2 if I2 has enough current (ie > threshold set by Z12T)
            F['sel'] = (abs(F['I012'][2])/abs(F['I012'][1]))>Z12T
        
            # Find Dir actual direction based on graph 
            if('_closein' in F['busNumber']):
                Fpoint = F['busNumber'].split('_closein')[1]
            else:
                Fpoint = F['busNumber']
            # Find Relay
            Fpath = nx.shortest_path(G,B1,Fpoint)
            if(B2 in Fpath):
                F['Dir_act'] = 1
            else:
                F['Dir_act'] = -1
        # if I is not > 0.001 (fault current is 0 ignore, cant trip or damage)
        else:
            F['Dir_act'] = 0
            F['sel'] = False
            F['Z0'] = np.nan
            F['Z1'] = np.nan
            
    # %% Calculate Z1 and Z0 for all forward and reverse faults
    AF =[x['Z1'] for x in RFaults if x['Dir_act']==1 and x['sel']==False]
    AR =[x['Z1'] for x in RFaults if x['Dir_act']==-1 and x['sel']==False]
    
    A0F =[x['Z0'] for x in RFaults if x['Dir_act']==1 and ('SLG' in x['FaultType'])]
    A0R =[x['Z0'] for x in RFaults if x['Dir_act']==-1 and ('SLG' in x['FaultType'])]
    
    #Find min max of angle Z1
    if(len(AF)>0 and len(AR)>0):
        AminmaxFR = [max(AF),min(AF),max(AR),min(AR)]
    elif(len(AF)>0 and len(AR)==0):
        AminmaxFR = [max(AF),min(AF),(max(AF)+180) % 360,(min(AF)+180) % 360]
    elif(len(AR)>0 and len(AF)==0):
        AminmaxFR = [(max(AR)+180) % 360,(min(AR)+180) % 360,max(AR),min(AR)]
    else:
        AminmaxFR = []
    # find min max of angles for Z0
    if(len(A0F)>0 and len(A0R)>0):
        A0minmaxFR = [max(A0F),min(A0F),max(A0R),min(A0R)]
    elif(len(A0F)>0 and len(A0R)==0):
        A0minmaxFR = [max(A0F),min(A0F),(max(A0F)+180) % 360,(min(A0F)+180) % 360]
    elif(len(A0R)>0 and len(A0F)==0):
        A0minmaxFR = [(max(A0R)+180) % 360,(min(A0R)+180) % 360,max(A0R),min(A0R)]
    else:
        A0minmaxFR = []
    # Find Z1Ang
    if(len(AminmaxFR)>0):
        if(abs(AminmaxFR[0]-AminmaxFR[2]) <= abs(AminmaxFR[1]-AminmaxFR[3])):
            x = AminmaxFR[0] + abs(AminmaxFR[0]-AminmaxFR[2])/2 - 90
        else:
            x = AminmaxFR[1] - abs(AminmaxFR[1]-AminmaxFR[3])/2 + 90
    else:
        x = 45;
    # Find Z0ANG
    if(len(A0minmaxFR)>0):
        if(abs(A0minmaxFR[0]-A0minmaxFR[2]) <= abs(A0minmaxFR[1]-A0minmaxFR[3])):
            x0= A0minmaxFR[0] + abs(A0minmaxFR[0]-A0minmaxFR[2])/2 - 90
        else:
            x0= A0minmaxFR[1] - abs(A0minmaxFR[1]-A0minmaxFR[3])/2 + 90
    else:
        x = 45;
        
    # angle is within 360
    x = x % 360
    x0 = x0 % 360
    
    # if not is 1 quad
    if(x > 90 and x<180): # Q2
        # Error 
        x = 90
        if(Show_Warnings==1):
            print('Error: DIR MTA in Q2 for Positive seq for '+Rname)
    elif(x>180 and x<270): #Q3
        # Dir Revrse 
        if(Show_Warnings==1):
            print('Error: DIR MTA in Q3 for Positive seq for '+Rname)
        x  =  (x-180)
    elif(x>270):
        if(Show_Warnings==1):
            print('Error: DIR MTA in Q2 for Positive seq for '+Rname)
        x = 90
        
    # if not is 1 quad
    if(x0 > 90 and x0<180): # Q2
        # Error 
        x0 = 90
        if(Show_Warnings==1):
            print('Error: DIR MTA in Q2 for Zero seq for '+Rname)
    elif(x0>180 and x0<270): #Q3
        # Dir Revrse 
        x0  =  (x0-180)
        if(Show_Warnings==1):
            print('Error: DIR MTA in Q3 for Zero seq for '+Rname)
    elif(x0>270):
        x0 = 90
        if(Show_Warnings==1):
            print('Error: DIR MTA in Q2 for Zero seq for '+Rname)
            
    # condition for relay settigns limits
    if(x<=90 and x>=0):
        if(x<1):
            x = 1
    if(x0<=90 and x0>=0):
        if(x0<1):
            x0 = 1
    
    Z1ANG = x
    Z0ANG = x0

    # %%
    for F in RFaults:
        Z2 = calc_Z2(Z1ANG,F,Z12T)
        if(np.isnan(Z2)):
            Z2 = calc_Z2(Z1ANG,F,0.0)
            F['Dir'] = None
            F['Z2'] = Z2
            F['Z2FT'] = None
            F['Z2RT'] = None
        else:
            F['Z2'] = Z2
    
    if(len([x['Z2'] for x in RFaults if x['Dir_act'] == -1 and x['sel']==True]) == 0 and len([x['Z2'] for x in RFaults if x['Dir_act'] == 1 and x['sel']==True]) > 0):
        Z2TH = max([x['Z2'] for x in RFaults if x['Dir_act'] == 1 and x['sel']==True])*1.5
    elif(len([x['Z2'] for x in RFaults if x['Dir_act'] == 1 and x['sel']==True]) == 0 and len([x['Z2'] for x in RFaults if x['Dir_act'] == -1 and x['sel']==True]) > 0):
        Z2TH = min([x['Z2'] for x in RFaults if x['Dir_act'] == -1 and x['sel']==True])*1.5
    elif(len([x['Z2'] for x in RFaults if x['Dir_act'] == -1 and x['sel']==True]) == 0 and len([x['Z2'] for x in RFaults if x['Dir_act'] == 1 and x['sel']==True]) == 0):
        Z2TH = 1.0
    else:
        Z2TH = (min([x['Z2'] for x in RFaults if x['Dir_act'] == -1 and x['sel']==True]) - max([x['Z2'] for x in RFaults if x['Dir_act'] == 1 and x['sel']==True])) + min([x['Z2'] for x in RFaults if x['Dir_act'] == -1 and x['sel']==True])
    
    
    Z2F = Z2TH
    Z2R = Z2TH
    
    Z1 = (Z1MAG,Z1ANG)
    Z0 = (Z0MAG,Z0ANG)
    Z2 = (Z2F,Z2R)
    return Z1,Z0,Z2,Z12T