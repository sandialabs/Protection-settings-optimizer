# import Steady state data

import csv
import numpy as np

def  read_Relays_CSV_Data(File):
    Data = []
    with open(File, mode='r') as inp:
        reader = csv.DictReader(inp)
        Data = list(reader)

    Device_Data_CSV = [dict.fromkeys(['RelayName','Ia_mag','Ia_ang','Va_mag','Va_ang','P','Q','switchState','PVs','In']) for number in range(len(Data))]
    pha  = [None]*12   
    for ii in range(len(Data)):
        
        if(ii<6): 
           Device_Data_CSV[ii]['RelayName'] = 'R'+str(ii+1) 
        else:
           Device_Data_CSV[ii]['RelayName']= 'RTL'+str(ii-5)
        
        
        p = Data[ii]['phasors']
        ph = p[1:-1].split(',')
        for jj in range(len(ph)):
            pha[jj] = float(ph[jj].replace(')','').replace('(',''))
        
        #Device_Data_CSV[ii]['RelayName'] = Data[ii]['RelayName']
        Device_Data_CSV[ii]['switchState'] = Data[ii]['digital'].lower() in ['[1]','true','close','closed']
        Device_Data_CSV[ii]['Va_mag'] = min(pha[0],pha[2],pha[4])
        Device_Data_CSV[ii]['Va_ang'] = pha[1]
        Device_Data_CSV[ii]['Ia_mag'] = max(pha[6],pha[8],pha[10])
        Device_Data_CSV[ii]['Ia_ang'] = pha[7]
        
        Device_Data_CSV[ii]['Vabc'] = [pdef(pha[0],pha[1]*(180/np.pi)),pdef(pha[2],pha[3]*(180/np.pi)),pdef(pha[4],pha[5]*(180/np.pi))]
        Device_Data_CSV[ii]['Iabc'] = [pdef(pha[6],pha[7]*(180/np.pi)),pdef(pha[8],pha[9]*(180/np.pi)),pdef(pha[10],pha[11]*(180/np.pi))]
        Device_Data_CSV[ii]['V012'] = abc2012(Device_Data_CSV[ii]['Vabc'])
        Device_Data_CSV[ii]['I012'] = abc2012(Device_Data_CSV[ii]['Iabc'])
        Device_Data_CSV[ii]['Z012'] = np.divide(Device_Data_CSV[ii]['V012'],Device_Data_CSV[ii]['I012'])
        
        Device_Data_CSV[ii]['In'] = abs(pdef(pha[6],pha[7]*(180/np.pi))+
                                        pdef(pha[8],pha[9]*(180/np.pi))+
                                        pdef(pha[10],pha[11]*(180/np.pi)))
        
        Sa = pdef(pha[0],pha[1]*(180/np.pi)) * pdef(pha[6],pha[7]*(180/np.pi))
        Sb = pdef(pha[2],pha[3]*(180/np.pi)) * pdef(pha[8],pha[9]*(180/np.pi))
        Sc = pdef(pha[4],pha[5]*(180/np.pi)) * pdef(pha[10],pha[11]*(180/np.pi))
        
        Device_Data_CSV[ii]['P'] = np.real(Sa+Sb+Sc)
        Device_Data_CSV[ii]['Q'] = np.imag(Sa+Sb+Sc)
        Device_Data_CSV[ii]['PVs'] = Data[ii]['PVStatus']
        
        #za = pdef(pha[0],pha[1]*(180/np.pi))/pdef(pha[6],pha[7]*(180/np.pi))
        #zb = pdef(pha[2],pha[3]*(180/np.pi))/pdef(pha[8],pha[9]*(180/np.pi))
        #zc = pdef(pha[4],pha[5]*(180/np.pi))/pdef(pha[10],pha[11]*(180/np.pi))
        
        #ZL1_mag = min(np.abs([za,zb,zc]))
        #ZL1_ang = np.mean(np.angle([za,zb,zc]))*(180/np.pi)
        #Device_Data_CSV[ii]['ZL1'] =  pdef(ZL1_mag,ZL1_ang)
    return  Device_Data_CSV


def read_Fault_CSV_Data(File):
    Data = []
    fieldnames = ['busNumber','Relay','FaultType','Ia_mag','Ia_ang','Ib_mag','Ib_ang','Ic_mag','Ic_ang','Va_mag','Va_ang','Vb_mag','Vb_ang','Vc_mag','Vc_ang','x3Iz_mag','Iz_ang']
    with open(File, mode='r') as inp:
        reader = csv.DictReader(inp,fieldnames=fieldnames,dialect='excel')
        Data = list(reader)
    
    
    Fault_Data_CSV = [dict.fromkeys(['busNumber','Relay','FaultType','Ia_mag','Ia_ang','Ib_mag','Ib_ang','Ic_mag','Ic_ang','Va_mag','Va_ang','Vb_mag','Vb_ang','Vc_mag','Vc_ang','In_mag','In_ang','P','Q','Z012']) for number in range(len(Data))]    
    for ii in range(len(Data)):
         Fault_Data_CSV[ii]['busNumber'] = Data[ii]['busNumber'].lstrip()
         Fault_Data_CSV[ii]['Relay'] = Data[ii]['Relay']
         Fault_Data_CSV[ii]['FaultType'] = Data[ii]['FaultType']
         Fault_Data_CSV[ii]['Ia_mag'] = float(Data[ii]['Ia_mag'])
         Fault_Data_CSV[ii]['Ia_ang'] = float(Data[ii]['Ia_ang'])
         Fault_Data_CSV[ii]['Ib_mag'] = float(Data[ii]['Ib_mag'])
         Fault_Data_CSV[ii]['Ib_ang'] = float(Data[ii]['Ib_ang'])
         Fault_Data_CSV[ii]['Ic_mag'] = float(Data[ii]['Ic_mag'])
         Fault_Data_CSV[ii]['Ic_ang'] = float(Data[ii]['Ic_ang'])
         Fault_Data_CSV[ii]['Va_mag'] = float(Data[ii]['Va_mag'])
         Fault_Data_CSV[ii]['Va_ang'] = float(Data[ii]['Va_ang'])
         Fault_Data_CSV[ii]['Vb_mag'] = float(Data[ii]['Vb_mag'])
         Fault_Data_CSV[ii]['Vb_ang'] = float(Data[ii]['Vb_ang'])
         Fault_Data_CSV[ii]['Vc_mag'] = float(Data[ii]['Vc_mag'])
         Fault_Data_CSV[ii]['Vc_ang'] = float(Data[ii]['Vc_ang'])
         Fault_Data_CSV[ii]['In_mag'] = float(Data[ii]['x3Iz_mag'])
         Fault_Data_CSV[ii]['In_ang'] = float(Data[ii]['Iz_ang'])
         
         Va = pdef(Fault_Data_CSV[ii]['Va_mag'],Fault_Data_CSV[ii]['Va_ang'])*(4.16e3/np.sqrt(3))
         Vb = pdef(Fault_Data_CSV[ii]['Vb_mag'],Fault_Data_CSV[ii]['Vb_ang'])*(4.16e3/np.sqrt(3))
         Vc = pdef(Fault_Data_CSV[ii]['Vc_mag'],Fault_Data_CSV[ii]['Vc_ang'])*(4.16e3/np.sqrt(3))
         
         Ia = pdef(Fault_Data_CSV[ii]['Ia_mag'],Fault_Data_CSV[ii]['Ia_ang'])
         Ib = pdef(Fault_Data_CSV[ii]['Ib_mag'],Fault_Data_CSV[ii]['Ib_ang'])
         Ic = pdef(Fault_Data_CSV[ii]['Ic_mag'],Fault_Data_CSV[ii]['Ic_ang'])
         
         V012 = abc2012([Va,Vb,Vc])
         I012 = abc2012([Ia,Ib,Ic])
         
         Fault_Data_CSV[ii]['Vabc'] = [Va,Vb,Vc]
         Fault_Data_CSV[ii]['Iabc'] = [Ia,Ib,Ic]
         Fault_Data_CSV[ii]['V012'] = V012
         Fault_Data_CSV[ii]['I012'] = I012
         
         Sa = Va * np.conj(Ia)
         Sb = Vb * np.conj(Ib)
         Sc = Vc * np.conj(Ic)                                            
         
         Fault_Data_CSV[ii]['P'] = np.real(Sa+Sb+Sc)
         Fault_Data_CSV[ii]['Q'] = np.imag(Sa+Sb+Sc)
         Fault_Data_CSV[ii]['Z012'] = [np.divide(V012[0],I012[0],where=I012[0]!=0),
                                       np.divide(V012[1],I012[1],where=I012[1]!=0),
                                       np.divide(V012[2],I012[2],where=I012[2]!=0)]
     
    return Fault_Data_CSV
        

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