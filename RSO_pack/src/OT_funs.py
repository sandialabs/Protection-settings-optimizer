import numpy as np

def OCOT(M,TDS,TCC,V):
    if(M <= 1.0):
        return float("inf")
    elif(M>30):
        M = 30
        
    TDS = TDS/10
    
    if(TCC>6 and V<=0.95):
        TCC = TCC-6
    elif(TCC>6 and V>0.95):
        T=float("inf")
        return T
        
    if(TCC == 1):   # U1 - MI
        T = TDS * (0.02260 + (0.01040/(M**0.02 - 1)))
    elif(TCC == 2): # U2 - I
        T = TDS * (0.18000 + (5.95000/(M**2.00 - 1)))
    elif(TCC == 3): # U3 - VI
        T = TDS * (0.09630 + (3.88000/(M**2.00 - 1)))
    elif(TCC == 4): # U4 - EI
        T = TDS * (0.03520 + (5.67000/(M**2.00 - 1)))
    elif(TCC == 5): # U5 - STI
        T = TDS * (0.00262 + (0.00342/(M**2.00 - 1)))
    elif(TCC == 6):
        T = max(TDS,0.02)
        
    # elif(TCC == 6): # C1 - SI
    #     T = TDS * (0.14/(M**0.02 - 1))
    # elif(TCC == 7): # C2 - VI
    #     T = TDS * (13.5/(M**1.00 - 1))
    # elif(TCC == 8): # C3 - EI
    #     T = TDS * (80.0/(M**2.00 - 1))
    # elif(TCC == 9): # C4 - LTI
    #     T = TDS * (120/(M**1.00 - 1))
    # elif(TCC == 10):# C5 - STI
    #     T = TDS * (0.05/(M**0.04 - 1))
    else:
        T = float("inf")
    
    if(T<0.02):
        T = 0.02+0.03
    
    return T+0.03

def mhoDist(V,I,Z1,Z2,T):
    Zmeas = np.divide(V,I)
    r = [abs(Z1)/2,abs(Z2)/2]
    zR = [Z1,Z2]
    
    xmeas = np.real(Zmeas)
    ymeas = np.imag(Zmeas)
    
    cx = np.real(zR)/2
    cy = np.imag(zR)/2
    
    r2=[0]*len(T)
    tripzone= [0]*len(T)
    
    for ii in range(len(V)):
        for jj in range(len(r)):
            r2[ii] = np.sqrt( (xmeas[ii]-cx[jj])**2 + (ymeas[ii]-cy[ii])**2 )
            if( r2[ii]<r[jj]):
                tripzone[ii] = jj
                break
            else:
                tripzone[ii] = -1
    
    for kk in range(len(V)):
        if(tripzone[ii]<0):
            Top = np.inf
        else:
            Top = T[tripzone[kk]]
    return Top

def OCIT(I,IT):
    if(abs(I)>IT and IT!=0):
        return 0.02
    else:
        return float("inf")
    
def OCTCC_Name(TCC):
    
    if(TCC>6):
        TCC = TCC-6
    
    if(TCC == 1):   # U1 - MI
        T = 'U1: moderately inverse (OC)'
    elif(TCC == 2): # U2 - I
        T = 'U2: inverse (OC)'
    elif(TCC == 3): # U3 - VI
        T = 'U3: very inverse (OC)'
    elif(TCC == 4): # U4 - EI
        T = 'U4: extremely inverse (OC)'
    elif(TCC == 5): # U5 - STI
        T = 'U5: short-time inverse (OC)'
    elif(TCC == 6):
        T = 'DT: Discrete time (OC)'
    # elif(TCC == 6): # C1 - SI
    #     T = 'C1'
    # elif(TCC == 7): # C2 - VI
    #     T = 'C2'
    # elif(TCC == 8): # C3 - EI
    #     T = 'C3'
    # elif(TCC == 9): # C4 - LTI
    #     T = 'C4'
    # elif(TCC == 10):# C5 - STI
    #     T = 'C5'
    else:
        T = 'Unknown'
    return T

def OCTCC_Num(TCC,VROC):
    
    if(TCC == 'U1: moderately inverse (OC)'):   # U1 - MI
        T = 1
    elif(TCC == 'U2: inverse (OC)'): # U2 - I
        T = 2
    elif(TCC == 'U3: very inverse (OC)'): # U3 - VI
        T = 3
    elif(TCC == 'U4: extremely inverse (OC)'): # U4 - EI
        T = 4
    elif(TCC == 'U5: short-time inverse (OC)'): # U5 - STI
        T = 5
    elif(TCC == 'DT: Discrete time (OC)'):
        T = 6
    # elif(TCC == 6): # C1 - SI
    #     T = 'C1'
    # elif(TCC == 7): # C2 - VI
    #     T = 'C2'
    # elif(TCC == 8): # C3 - EI
    #     T = 'C3'
    # elif(TCC == 9): # C4 - LTI
    #     T = 'C4'
    # elif(TCC == 10):# C5 - STI
    #     T = 'C5'
    else:
        T = 'Unknown'
        
    if(VROC):
        T = T+6
    return T