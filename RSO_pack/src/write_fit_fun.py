# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 15:38:52 2022

@author: maste
"""

def write_GA_fit_fun(CtrlDev,ProDevRef,Pairs,M2,fundir,fileName,CTI,OTmax):
    
    fileID = open(fundir+'\\'+fileName+'.py','w+')
    
    fileID.write('\n') 
    fileID.write('from RSO_pack import OCOT\n') 
    fileID.write('from RSO_pack import OCIT\n') 
    fileID.write('\n') 
    fileID.write('def fitness_func(solution, solution_idx):')
    fileID.write('\n')
    
    for ii in range(len(ProDevRef)):
        if(ProDevRef[ii]['Type'] == 'Relay_TOC'):
            
            fileID.write('\tTDS'+str(ii)+' = solution['+str(CtrlDev[ii]-4)+']\n')
            fileID.write('\tTDSg'+str(ii)+' = solution['+str(CtrlDev[ii]-3)+']\n')
            fileID.write('\trelayType'+str(ii)+' = solution['+str(CtrlDev[ii]-2)+']\n')
            fileID.write('\trelayTypeg'+str(ii)+' = solution['+str(CtrlDev[ii]-1)+']\n')
            
        elif(ProDevRef[ii]['Type'] == 'Relay_DT'):
            
            fileID.write('\tTDS'+str(ii)+' = solution['+str(CtrlDev[ii]-4)+']\n')
            fileID.write('\tTDSg'+str(ii)+' = solution['+str(CtrlDev[ii]-3)+']\n')
            fileID.write('\trelayType'+str(ii)+' = solution['+str(CtrlDev[ii]-2)+']\n')
            fileID.write('\trelayTypeg'+str(ii)+' = solution['+str(CtrlDev[ii]-1)+']\n')
        elif('Rec_' in ProDevRef[ii]['Type']):
            fileID.write('\tTDS'+str(ii)+' = solution['+str(CtrlDev[ii]-4)+']\n')
            fileID.write('\tTDSg'+str(ii)+' = solution['+str(CtrlDev[ii]-3)+']\n')
            fileID.write('\trelayType'+str(ii)+' = solution['+str(CtrlDev[ii]-2)+']\n')
            fileID.write('\trelayTypeg'+str(ii)+' = solution['+str(CtrlDev[ii]-1)+']\n')
        else:
            pass
        
    fileID.write('\n')
    fileID.write('\tTtot=[0]*8*'+str(len(ProDevRef))+'\n')
    
    cc = 0
    for ii in range(len(ProDevRef)):
        if('Relay_' in ProDevRef[ii]['Type'] or 'Rec_' in ProDevRef[ii]['Type']):
            if(ProDevRef[ii]['Vmax3ph']<0.95 and ProDevRef[ii]['Imax3ph']/ProDevRef[ii]['Ip'] > 1):
                fileID.write('\tTtot['+str(cc)+'] = OCOT('+str(ProDevRef[ii]['Imax3ph']/ProDevRef[ii]['Ip'])+',TDS'+str(ii)+',relayType'+str(ii)+','+str(ProDevRef[ii]['Vmax3ph'])+')\n')
                cc=cc+1 # ('+str(cc)+')
            if(ProDevRef[ii]['Vmin3ph']<0.95 and ProDevRef[ii]['Imin3ph']/ProDevRef[ii]['Ip'] > 1):
                fileID.write('\tTtot['+str(cc)+'] = OCOT('+str(ProDevRef[ii]['Imin3ph']/ProDevRef[ii]['Ip'])+',TDS'+str(ii)+',relayType'+str(ii)+','+str(ProDevRef[ii]['Vmin3ph'])+')\n')
                cc=cc+1
            if(ProDevRef[ii]['VmaxLL']<0.95 and ProDevRef[ii]['ImaxLL']/ProDevRef[ii]['Ip'] > 1):
                fileID.write('\tTtot['+str(cc)+'] = OCOT('+str(ProDevRef[ii]['ImaxLL']/ProDevRef[ii]['Ip'])+',TDS'+str(ii)+',relayType'+str(ii)+','+str(ProDevRef[ii]['VmaxLL'])+')\n')
                cc=cc+1                
            if(ProDevRef[ii]['VminLL']<0.95 and ProDevRef[ii]['IminLL']/ProDevRef[ii]['Ip'] > 1):
                fileID.write('\tTtot['+str(cc)+'] = OCOT('+str(ProDevRef[ii]['IminLL']/ProDevRef[ii]['Ip'])+',TDS'+str(ii)+',relayType'+str(ii)+','+str(ProDevRef[ii]['VminLL'])+')\n')
                cc=cc+1                
            if(ProDevRef[ii]['VmaxSLG']<0.95):
                if(ProDevRef[ii]['ImaxSLG']/ProDevRef[ii]['Ip'] > 1):
                    fileID.write('\tTtot['+str(cc)+'] = OCOT('+str(ProDevRef[ii]['ImaxSLG']/ProDevRef[ii]['Ip'])+',TDS'+str(ii)+',relayType'+str(ii)+','+str(ProDevRef[ii]['VmaxSLG'])+')\n')
                    cc=cc+1
                if(ProDevRef[ii]['IgmaxSLG']/ProDevRef[ii]['Inp'] > 1):
                    fileID.write('\tTtot['+str(cc)+'] = OCOT('+str(ProDevRef[ii]['IgmaxSLG']/ProDevRef[ii]['Inp'])+',TDSg'+str(ii)+',relayTypeg'+str(ii)+','+str(ProDevRef[ii]['VmaxSLG'])+')\n')
                    cc=cc+1   
            if(ProDevRef[ii]['VminSLG']<0.95):
                if(ProDevRef[ii]['IminSLG']/ProDevRef[ii]['Ip'] > 1):
                    fileID.write('\tTtot['+str(cc)+'] = OCOT('+str(ProDevRef[ii]['IminSLG']/ProDevRef[ii]['Ip'])+',TDS'+str(ii)+',relayType'+str(ii)+','+str(ProDevRef[ii]['VminSLG'])+')\n')
                    cc=cc+1
                if(ProDevRef[ii]['IgminSLG']/ProDevRef[ii]['Inp'] > 1):
                    fileID.write('\tTtot['+str(cc)+'] = OCOT('+str(ProDevRef[ii]['IgminSLG']/ProDevRef[ii]['Inp'])+',TDSg'+str(ii)+',relayTypeg'+str(ii)+','+str(ProDevRef[ii]['VminSLG'])+')\n')
                    cc=cc+1
     
    fileID.write('\tdel Ttot['+str(cc-1)+':]\n')
    fileID.write('\ty1=sum(Ttot)\n')
    fileID.write('\ty1_1 = 0 if max(Ttot)<='+str(OTmax)+' else max(Ttot)*1000 \n')
    #fileID.write('\ty1=(sum(Ttot)/len(Ttot))*10 + max(Ttot)*10 \n')
    fileID.write('\n')
    
    # %% Write Cons
    fileID.write('\n')
    fileID.write('\tgac = [None]*6*'+str(len(Pairs))+'*6\n')
    pp = 0
    for ii in range(len(Pairs)):
        Farray = [x for x in M2 if (x[1]==Pairs[ii][6] and x[20] == 1)] 
        for jj in range(len(Farray)):
            Bacdat = [x for x in M2 if (x[1]==Pairs[ii][7] and x[0]==Farray[jj][0])]
            P_ind = [x['Oind'] for x in ProDevRef].index(Pairs[ii][6])
            B_ind = [x['Oind'] for x in ProDevRef].index(Pairs[ii][7])
            fileID.write('# cons for Pair:'+str(ii)+' Pind='+str(P_ind)+' Bind='+str(B_ind)+'\n')
            FI = [2,3,8,9,14,15]
            VI = [4,5,10,11,16,17]
            NI = [6,7,12,13,18,19]
            for FF in range(0,6):
                if(Farray[jj][FI[FF]] == None or Bacdat[0][FI[FF]] == None):
                    continue
                else:
                    fileID.write('# cons for jj:'+str(jj)+' FF:'+str(FF)+'\n') 
                    Test_Val = [Farray[jj][FI[FF]]/ProDevRef[P_ind]['Ip'],Farray[jj][VI[FF]],Bacdat[0][FI[FF]]/ProDevRef[B_ind]['Ip'],Bacdat[0][VI[FF]]] 
                     
                    if( (Farray[jj][FI[FF]]/ProDevRef[P_ind]['Ip']>=1.2) and
                        (Bacdat[0][FI[FF]]/ProDevRef[B_ind]['Ip']>=1.2) ):
                            
                        if(Pairs[ii][4] == 'Relay_TOC'):
                            fileID.write('\tT1 =  OCOT('+str(Farray[jj][FI[FF]]/ProDevRef[P_ind]['Ip'])+',TDS'+str(P_ind)+',relayType'+str(P_ind)+','+str(Farray[jj][VI[FF]])+')\n')
                            fileID.write('\tT1g = OCOT('+str(Farray[jj][NI[FF]]/ProDevRef[P_ind]['Inp'])+',TDSg'+str(P_ind)+',relayTypeg'+str(P_ind)+','+str(Farray[jj][VI[FF]])+')\n')
                            fileID.write('\tT1_IT = OCIT('+str(Farray[jj][FI[FF]])+','+str(ProDevRef[P_ind]['IT'])+')\n')
                            fileID.write('\tT1_ITg = OCIT('+str(Farray[jj][NI[FF]])+','+str(ProDevRef[P_ind]['ITg'])+')\n')
                            fileID.write('\tTpri = [x for x in [T1,T1g,T1_IT,T1_ITg] if x>0]\n')
                        elif(Pairs[ii][4] == 'Rec'):
                            FP_ind = P_ind
                            SP_ind = P_ind+1
                            fileID.write('\tT1F =  OCOT('+str(Farray[jj][FI[FF]]/ProDevRef[FP_ind]['Ip'])+',TDS'+str(FP_ind)+',relayType'+str(FP_ind)+','+str(Farray[jj][VI[FF]])+')\n')
                            fileID.write('\tT1gF = OCOT('+str(Farray[jj][NI[FF]]/ProDevRef[FP_ind]['Inp'])+',TDSg'+str(FP_ind)+',relayTypeg'+str(FP_ind)+','+str(Farray[jj][VI[FF]])+')\n')
                            fileID.write('\tT1_IT = OCIT('+str(Farray[jj][FI[FF]])+','+str(ProDevRef[FP_ind]['IT'])+')\n')
                            fileID.write('\tT1_ITg = OCIT('+str(Farray[jj][NI[FF]])+','+str(ProDevRef[FP_ind]['ITg'])+')\n')
                            fileID.write('\tTpriF = [x for x in [T1F,T1gF,T1_IT,T1_ITg] if x>0]\n')  
                            
                            fileID.write('\tT1S =  OCOT('+str(Farray[jj][FI[FF]]/ProDevRef[SP_ind]['Ip'])+',TDS'+str(SP_ind)+',relayType'+str(SP_ind)+','+str(Farray[jj][VI[FF]])+')\n')
                            fileID.write('\tT1gS = OCOT('+str(Farray[jj][NI[FF]]/ProDevRef[SP_ind]['Inp'])+',TDSg'+str(SP_ind)+',relayTypeg'+str(SP_ind)+','+str(Farray[jj][VI[FF]])+')\n')
                            fileID.write('\tTpriS = [x for x in [T1S,T1gS,T1_IT,T1_ITg] if x>0]\n')  
                            
                            
                        if(Pairs[ii][5] == 'Relay_TOC'):
                            fileID.write('\tT2 = OCOT('+str(Bacdat[0][FI[FF]]/ProDevRef[B_ind]['Ip'])+',TDS'+str(B_ind)+',relayType'+str(B_ind)+','+str(Bacdat[0][VI[FF]])+')\n' )
                            fileID.write('\tT2g = OCOT('+str(Bacdat[0][NI[FF]]/ProDevRef[B_ind]['Inp'])+',TDSg'+str(B_ind)+',relayTypeg'+str(B_ind)+','+str(Bacdat[0][VI[FF]])+')\n' )
                            fileID.write('\tT2_IT = OCIT('+str(Bacdat[0][FI[FF]])+','+str(ProDevRef[B_ind]['IT'])+')\n')
                            fileID.write('\tT2_ITg = OCIT('+str(Bacdat[0][NI[FF]])+','+str(ProDevRef[B_ind]['ITg'])+')\n')
                            fileID.write('\tTbac = [x for x in [T2,T2g,T2_IT,T2_ITg] if x>0]\n')
                        elif(Pairs[ii][5] == 'Rec'):
                            FB_ind = B_ind
                            SB_ind = B_ind+1
                            fileID.write('\tT2F = OCOT('+str(Bacdat[0][FI[FF]]/ProDevRef[FB_ind]['Ip'])+',TDS'+str(FB_ind)+',relayType'+str(FB_ind)+','+str(Bacdat[0][VI[FF]])+')\n' )
                            fileID.write('\tT2gF = OCOT('+str(Bacdat[0][NI[FF]]/ProDevRef[FB_ind]['Inp'])+',TDSg'+str(FB_ind)+',relayTypeg'+str(FB_ind)+','+str(Bacdat[0][VI[FF]])+')\n' )
                            fileID.write('\tT2_IT = OCIT('+str(Bacdat[0][FI[FF]])+','+str(ProDevRef[FB_ind]['IT'])+')\n')
                            fileID.write('\tT2_ITg = OCIT('+str(Bacdat[0][NI[FF]])+','+str(ProDevRef[FB_ind]['ITg'])+')\n')
                            fileID.write('\tTbacF = [x for x in [T2F,T2gF,T2_IT,T2_ITg] if x>0]\n')
                            
                            fileID.write('\tT2S = OCOT('+str(Bacdat[0][FI[FF]]/ProDevRef[SB_ind]['Ip'])+',TDS'+str(SB_ind)+',relayType'+str(SB_ind)+','+str(Bacdat[0][VI[FF]])+')\n' )
                            fileID.write('\tT2gS = OCOT('+str(Bacdat[0][NI[FF]]/ProDevRef[SB_ind]['Inp'])+',TDSg'+str(SB_ind)+',relayTypeg'+str(SB_ind)+','+str(Bacdat[0][VI[FF]])+')\n' )
                            fileID.write('\tTbacS = [x for x in [T2S,T2gS,T2_IT,T2_ITg] if x>0]\n')
                        
                        
                        if(Pairs[ii][4] == 'Relay_TOC' and Pairs[ii][5] == 'Relay_TOC'):
                            fileID.write('\tgac['+str(pp)+'] = min(Tbac) - (min(Tpri) + '+str(CTI)+')\n')
                            pp = pp + 1
                        elif(Pairs[ii][4] == 'Relay_TOC' and Pairs[ii][5] == 'Rec'):
                            fileID.write('\tgac['+str(pp)+'] = min(TbacF) - (min(Tpri) + '+str(CTI)+')\n')
                            fileID.write('\tgac['+str(pp+1)+'] = min(TbacS) - (min(TbacF) + '+str(0.2)+')\n')
                            pp = pp + 2
                        elif(Pairs[ii][4] == 'Rec' and Pairs[ii][5] == 'Relay_TOC'):
                            fileID.write('\tgac['+str(pp)+'] = min(Tbac) - (min(TpriS) + '+str(CTI)+')\n')
                            fileID.write('\tgac['+str(pp+1)+'] = min(TpriS) - (min(TpriF) + '+str(0.2)+')\n')
                            pp = pp + 2
                        
                        fileID.write('\n')
                    else:
                        fileID.write('# Test:'+str(Test_Val)+'\n')
                        fileID.write('# Farr:'+str(Farray[jj])+'\n')
                        fileID.write('# Barr:'+str(Bacdat[0])+'\n')


    fileID.write('\tdel gac['+str(pp-1)+':]\n')
    fileID.write('\n')
    fileID.write('\tPen = 0\n')
    fileID.write('\tfor i in gac:\n')
    fileID.write('\t\tif(i>=0):\n')
    fileID.write('\t\t\tPen =  Pen + i*0.01\n')
    fileID.write('\t\telse:\n')
    fileID.write('\t\t\tPen =Pen + (-100000*i)\n')
    fileID.write('\n')        
    fileID.write('\ty2=min(gac) if min(gac)>=0 else -min(gac)*1000\n')
    fileID.write('\n')
    
    #fileID.write('\tprint(Ttot)\n')
    #fileID.write('\tprint(\'\\n\')\n')
    fileID.write('\tZ = -(y1+y1_1+y2+Pen)\n')
    fileID.write('\treturn Z')
    fileID.close()
    return fileID



# %% CTI and Constraint Calculation
def write_Con(CtrlDev,ProDevRef,Pairs,M2,fundir,fileName,CTI):
    
    fileID = open(fundir+'\\'+fileName+'_con.py','w+')
    
    fileID.write('\n') 
    fileID.write('from RSO_pack import OCOT\n') 
    fileID.write('from RSO_pack import OCIT\n') 
    fileID.write('\n') 
    fileID.write('def Con_func(solution, solution_idx):')
    fileID.write('\n')
    
    for ii in range(len(ProDevRef)):
        if(ProDevRef[ii]['Type'] == 'Relay_TOC'):
            
            fileID.write('\tTDS'+str(ii)+' = solution['+str(CtrlDev[ii]-4)+']\n')
            fileID.write('\tTDSg'+str(ii)+' = solution['+str(CtrlDev[ii]-3)+']\n')
            fileID.write('\trelayType'+str(ii)+' = solution['+str(CtrlDev[ii]-2)+']\n')
            fileID.write('\trelayTypeg'+str(ii)+' = solution['+str(CtrlDev[ii]-1)+']\n')
            
        elif(ProDevRef[ii]['Type'] == 'Relay_DT'):
            
            fileID.write('\tTDS'+str(ii)+' = solution['+str(CtrlDev[ii]-4)+']\n')
            fileID.write('\tTDSg'+str(ii)+' = solution['+str(CtrlDev[ii]-3)+']\n')
            fileID.write('\trelayType'+str(ii)+' = solution['+str(CtrlDev[ii]-2)+']\n')
            fileID.write('\trelayTypeg'+str(ii)+' = solution['+str(CtrlDev[ii]-1)+']\n')
        elif('Rec_' in ProDevRef[ii]['Type']):
            fileID.write('\tTDS'+str(ii)+' = solution['+str(CtrlDev[ii]-4)+']\n')
            fileID.write('\tTDSg'+str(ii)+' = solution['+str(CtrlDev[ii]-3)+']\n')
            fileID.write('\trelayType'+str(ii)+' = solution['+str(CtrlDev[ii]-2)+']\n')
            fileID.write('\trelayTypeg'+str(ii)+' = solution['+str(CtrlDev[ii]-1)+']\n')
        else:
            pass
        
    fileID.write('\n')
    fileID.write('\tTtot=[0]*8*'+str(len(ProDevRef))+'\n')
    
    cc = 0
    for ii in range(len(ProDevRef)):
        if('Relay_' in ProDevRef[ii]['Type'] or 'Rec_' in ProDevRef[ii]['Type']):
            if(ProDevRef[ii]['Vmax3ph']<0.95 and ProDevRef[ii]['Imax3ph']/ProDevRef[ii]['Ip']>1):
                fileID.write('\tTtot['+str(cc)+'] = OCOT('+str(ProDevRef[ii]['Imax3ph']/ProDevRef[ii]['Ip'])+',TDS'+str(ii)+',relayType'+str(ii)+','+str(ProDevRef[ii]['Vmax3ph'])+')\n')
                cc=cc+1 # ('+str(cc)+')
            if(ProDevRef[ii]['Vmin3ph']<0.95 and ProDevRef[ii]['Imin3ph']/ProDevRef[ii]['Ip']>1):
                fileID.write('\tTtot['+str(cc)+'] = OCOT('+str(ProDevRef[ii]['Imin3ph']/ProDevRef[ii]['Ip'])+',TDS'+str(ii)+',relayType'+str(ii)+','+str(ProDevRef[ii]['Vmin3ph'])+')\n')
                cc=cc+1
            if(ProDevRef[ii]['VmaxLL']<0.95 and ProDevRef[ii]['ImaxLL']/ProDevRef[ii]['Ip'] > 1):
                fileID.write('\tTtot['+str(cc)+'] = OCOT('+str(ProDevRef[ii]['ImaxLL']/ProDevRef[ii]['Ip'])+',TDS'+str(ii)+',relayType'+str(ii)+','+str(ProDevRef[ii]['VmaxLL'])+')\n')
                cc=cc+1                
            if(ProDevRef[ii]['VminLL']<0.95 and ProDevRef[ii]['IminLL']/ProDevRef[ii]['Ip'] > 1):
                fileID.write('\tTtot['+str(cc)+'] = OCOT('+str(ProDevRef[ii]['IminLL']/ProDevRef[ii]['Ip'])+',TDS'+str(ii)+',relayType'+str(ii)+','+str(ProDevRef[ii]['VminLL'])+')\n')
                cc=cc+1                
            if(ProDevRef[ii]['VmaxSLG']<0.95):
                if(ProDevRef[ii]['ImaxSLG']/ProDevRef[ii]['Ip'] > 1):
                    fileID.write('\tTtot['+str(cc)+'] = OCOT('+str(ProDevRef[ii]['ImaxSLG']/ProDevRef[ii]['Ip'])+',TDS'+str(ii)+',relayType'+str(ii)+','+str(ProDevRef[ii]['VmaxSLG'])+')\n')
                    cc=cc+1
                if(ProDevRef[ii]['IgmaxSLG']/ProDevRef[ii]['Inp'] > 1):
                    fileID.write('\tTtot['+str(cc)+'] = OCOT('+str(ProDevRef[ii]['IgmaxSLG']/ProDevRef[ii]['Inp'])+',TDSg'+str(ii)+',relayTypeg'+str(ii)+','+str(ProDevRef[ii]['VmaxSLG'])+')\n')
                    cc=cc+1   
            if(ProDevRef[ii]['VminSLG']<0.95):
                if(ProDevRef[ii]['IminSLG']/ProDevRef[ii]['Ip'] > 1):
                    fileID.write('\tTtot['+str(cc)+'] = OCOT('+str(ProDevRef[ii]['IminSLG']/ProDevRef[ii]['Ip'])+',TDS'+str(ii)+',relayType'+str(ii)+','+str(ProDevRef[ii]['VminSLG'])+')\n')
                    cc=cc+1
                if(ProDevRef[ii]['IgminSLG']/ProDevRef[ii]['Inp'] > 1):
                    fileID.write('\tTtot['+str(cc)+'] = OCOT('+str(ProDevRef[ii]['IgminSLG']/ProDevRef[ii]['Inp'])+',TDSg'+str(ii)+',relayTypeg'+str(ii)+','+str(ProDevRef[ii]['VminSLG'])+')\n')
                    cc=cc+1
                    
    fileID.write('\tdel Ttot['+str(cc-1)+':]\n')
    fileID.write('\n')
    
    # Write Cons
    fileID.write('\n')
    fileID.write('\tgac = [None]*6*'+str(len(Pairs))+'*6\n')
    pp = 0
    for ii in range(len(Pairs)):
        Farray = [x for x in M2 if (x[1]==Pairs[ii][6] and x[20] == 1)] 
        for jj in range(len(Farray)):
            Bacdat = [x for x in M2 if (x[1]==Pairs[ii][7] and x[0]==Farray[jj][0])]
            P_ind = [x['Oind'] for x in ProDevRef].index(Pairs[ii][6])
            B_ind = [x['Oind'] for x in ProDevRef].index(Pairs[ii][7])
            fileID.write('# cons for Pair:'+str(ii)+' Pind='+str(P_ind)+' Bind='+str(B_ind)+'\n')
            FI = [2,3,8,9,14,15]
            VI = [4,5,10,11,16,17]
            NI = [6,7,12,13,18,19]
            for FF in range(0,6):
                if(Farray[jj][FI[FF]] == None or Bacdat[0][FI[FF]] == None):
                    continue
                else:
                    fileID.write('# cons for jj:'+str(jj)+' FF:'+str(FF)+'\n') 
                    Test_Val = [Farray[jj][FI[FF]]/ProDevRef[P_ind]['Ip'],Farray[jj][VI[FF]],Bacdat[0][FI[FF]]/ProDevRef[B_ind]['Ip'],Bacdat[0][VI[FF]]] 
                     
                    if( (Farray[jj][FI[FF]]/ProDevRef[P_ind]['Ip']>1.2 and  Farray[jj][VI[FF]] < 0.95) and
                        (Bacdat[0][FI[FF]]/ProDevRef[B_ind]['Ip']>1.2 and Bacdat[0][VI[FF]] < 0.95 ) ):
                            
                        if(Pairs[ii][4] == 'Relay_TOC'):
                            fileID.write('\tT1 =  OCOT('+str(Farray[jj][FI[FF]]/ProDevRef[P_ind]['Ip'])+',TDS'+str(P_ind)+',relayType'+str(P_ind)+','+str(Farray[jj][VI[FF]])+')\n')
                            fileID.write('\tT1g = OCOT('+str(Farray[jj][NI[FF]]/ProDevRef[P_ind]['Inp'])+',TDSg'+str(P_ind)+',relayTypeg'+str(P_ind)+','+str(Farray[jj][VI[FF]])+')\n')
                            fileID.write('\tT1_IT = OCIT('+str(Farray[jj][FI[FF]])+','+str(ProDevRef[P_ind]['IT'])+')\n')
                            fileID.write('\tT1_ITg = OCIT('+str(Farray[jj][NI[FF]])+','+str(ProDevRef[P_ind]['ITg'])+')\n')
                            fileID.write('\tTpri = [x for x in [T1,T1g,T1_IT,T1_ITg] if x>0]\n')
                        elif(Pairs[ii][4] == 'Rec'):
                            FP_ind = P_ind
                            SP_ind = P_ind+1
                            fileID.write('\tT1F =  OCOT('+str(Farray[jj][FI[FF]]/ProDevRef[FP_ind]['Ip'])+',TDS'+str(FP_ind)+',relayType'+str(FP_ind)+','+str(Farray[jj][VI[FF]])+')\n')
                            fileID.write('\tT1gF = OCOT('+str(Farray[jj][NI[FF]]/ProDevRef[FP_ind]['Inp'])+',TDSg'+str(FP_ind)+',relayTypeg'+str(FP_ind)+','+str(Farray[jj][VI[FF]])+')\n')
                            fileID.write('\tT1_IT = OCIT('+str(Farray[jj][FI[FF]])+','+str(ProDevRef[FP_ind]['IT'])+')\n')
                            fileID.write('\tT1_ITg = OCIT('+str(Farray[jj][NI[FF]])+','+str(ProDevRef[FP_ind]['ITg'])+')\n')
                            fileID.write('\tTpriF = [x for x in [T1F,T1gF,T1_IT,T1_ITg] if x>0]\n')  
                            
                            fileID.write('\tT1S =  OCOT('+str(Farray[jj][FI[FF]]/ProDevRef[SP_ind]['Ip'])+',TDS'+str(SP_ind)+',relayType'+str(SP_ind)+','+str(Farray[jj][VI[FF]])+')\n')
                            fileID.write('\tT1gS = OCOT('+str(Farray[jj][NI[FF]]/ProDevRef[SP_ind]['Inp'])+',TDSg'+str(SP_ind)+',relayTypeg'+str(SP_ind)+','+str(Farray[jj][VI[FF]])+')\n')
                            fileID.write('\tTpriS = [x for x in [T1S,T1gS,T1_IT,T1_ITg] if x>0]\n')  
                            
                            
                        if(Pairs[ii][5] == 'Relay_TOC'):
                            fileID.write('\tT2 = OCOT('+str(Bacdat[0][FI[FF]]/ProDevRef[B_ind]['Ip'])+',TDS'+str(B_ind)+',relayType'+str(B_ind)+','+str(Bacdat[0][VI[FF]])+')\n' )
                            fileID.write('\tT2g = OCOT('+str(Bacdat[0][NI[FF]]/ProDevRef[B_ind]['Inp'])+',TDSg'+str(B_ind)+',relayTypeg'+str(B_ind)+','+str(Bacdat[0][VI[FF]])+')\n' )
                            fileID.write('\tT2_IT = OCIT('+str(Bacdat[0][FI[FF]])+','+str(ProDevRef[B_ind]['IT'])+')\n')
                            fileID.write('\tT2_ITg = OCIT('+str(Bacdat[0][NI[FF]])+','+str(ProDevRef[B_ind]['ITg'])+')\n')
                            fileID.write('\tTbac = [x for x in [T2,T2g,T2_IT,T2_ITg] if x>0]\n')
                        elif(Pairs[ii][5] == 'Rec'):
                            FB_ind = B_ind
                            SB_ind = B_ind+1
                            fileID.write('\tT2F = OCOT('+str(Bacdat[0][FI[FF]]/ProDevRef[FB_ind]['Ip'])+',TDS'+str(FB_ind)+',relayType'+str(FB_ind)+','+str(Bacdat[0][VI[FF]])+')\n' )
                            fileID.write('\tT2gF = OCOT('+str(Bacdat[0][NI[FF]]/ProDevRef[FB_ind]['Inp'])+',TDSg'+str(FB_ind)+',relayTypeg'+str(FB_ind)+','+str(Bacdat[0][VI[FF]])+')\n' )
                            fileID.write('\tT2_IT = OCIT('+str(Bacdat[0][FI[FF]])+','+str(ProDevRef[FB_ind]['IT'])+')\n')
                            fileID.write('\tT2_ITg = OCIT('+str(Bacdat[0][NI[FF]])+','+str(ProDevRef[FB_ind]['ITg'])+')\n')
                            fileID.write('\tTbacF = [x for x in [T2F,T2gF,T2_IT,T2_ITg] if x>0]\n')
                            
                            fileID.write('\tT2S = OCOT('+str(Bacdat[0][FI[FF]]/ProDevRef[SB_ind]['Ip'])+',TDS'+str(SB_ind)+',relayType'+str(SB_ind)+','+str(Bacdat[0][VI[FF]])+')\n' )
                            fileID.write('\tT2gS = OCOT('+str(Bacdat[0][NI[FF]]/ProDevRef[SB_ind]['Inp'])+',TDSg'+str(SB_ind)+',relayTypeg'+str(SB_ind)+','+str(Bacdat[0][VI[FF]])+')\n' )
                            fileID.write('\tTbacS = [x for x in [T2S,T2gS,T2_IT,T2_ITg] if x>0]\n')
                        
                        
                        if(Pairs[ii][4] == 'Relay_TOC' and Pairs[ii][5] == 'Relay_TOC'):
                            fileID.write('\tgac['+str(pp)+'] = min(Tbac) - (min(Tpri) + '+str(CTI)+')\n')
                            pp = pp + 1
                        elif(Pairs[ii][4] == 'Relay_TOC' and Pairs[ii][5] == 'Rec'):
                            fileID.write('\tgac['+str(pp)+'] = min(TbacF) - (min(Tpri) + '+str(CTI)+')\n')
                            fileID.write('\tgac['+str(pp+1)+'] = min(TbacS) - (min(TbacF) + '+str(0.2)+')\n')
                            pp = pp + 2
                        elif(Pairs[ii][4] == 'Rec' and Pairs[ii][5] == 'Relay_TOC'):
                            fileID.write('\tgac['+str(pp)+'] = min(Tbac) - (min(TpriS) + '+str(CTI)+')\n')
                            fileID.write('\tgac['+str(pp+1)+'] = min(TpriS) - (min(TpriF) + '+str(0.2)+')\n')
                            pp = pp + 2
                        
                        fileID.write('\n')
                    else:
                        fileID.write('# Test:'+str(Test_Val)+'\n')
                        fileID.write('# Farr:'+str(Farray[jj])+'\n')
                        fileID.write('# Barr:'+str(Bacdat[0])+'\n')


    fileID.write('\tdel gac['+str(pp-1)+':]\n')
    fileID.write('\n')
    fileID.write('\tPen = [0]*len(gac)\n')
    #fileID.write('\tPen = 0\n')
    fileID.write('\tfor i in range(len(gac)):\n')
    fileID.write('\t\tif(i>=0):\n')
    fileID.write('\t\t\tPen[i] =  i*0.01\n')
    fileID.write('\t\telse:\n')
    fileID.write('\t\t\tPen[i] = (i*-10000)\n')

    fileID.write('\treturn (Pen,gac,Ttot)')
    fileID.close()
    return fileID

    