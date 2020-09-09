#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 12:10:10 2020

@author: mdzik
"""

import CLB.CLBXMLWriter as CLBXML
import CLB.CLBXMLHandler
import CLB.VTIFile
# import bearded_octo_wookie.jednostki2 as j2
import os
import numpy as np
import scipy.optimize as so
import scipy.integrate as sint
import glob
import re
import CLB.VTIFile
import pandas as pd
import matplotlib.pyplot as plt

idx = 0
#R0 = 2.2  # Basic Reproduction Number - the number of secondary infections each infected individual produces.
#T_rec = 5.3  # days to recovery
T_rec = 1.
initial_phi = 0.10  
lambda_ph = 0.1

L2 = list()
for lbdt in 2.**-np.arange(1,9):
    tc = 100
    def getXML(**kwars):
            global idx
            
            idx = idx + 1
            prefix = '/tmp/id%03d/'%idx
            if 'clear' in kwars and kwars['clear']:
                os.system('rm -r %s'%prefix)
    
            os.system('mkdir %s'%prefix)
        
            CLBc = CLBXML.CLBConfigWriter( )
            fname = prefix+"run"
            CLBc.addGeomParam('nx', 16)
            CLBc.addGeomParam('ny', 16)
            
            CLBc.addTRT_SOI()
            CLBc.addBox()
            
            # CLBc.addSmoothing()
            # CLBc.addBox()
            
            params = {
            		"diffusivity_phi":0.1666666, #unimportant, no space variablitiy
            		"lambda":lambda_ph*lbdt,
                    "magic_parameter": 0.1666666,
            		"Init_PhaseField":initial_phi ,	
            		"phase_field_smoothing_coeff":0.1,
            }
            
            CLBc.addModelParams(params)
                     
      
            current = 0
            #for stop in np.logspace(0, np.log10(tc/lbdt), 100):    
            for stop in np.linspace(0, tc/lbdt, 101)[1:]:    
                
                CLBc.addSolve(iterations=stop-current)
                CLBc.addVTK()
                current = stop
           
            #CLBc.addSolve(iterations=tc/lbdt, vtk=50)
            CLBc.write(fname+'.xml')
    
    
            return prefix
        
        
    
    d0 = getXML(clear=True)
    
    wdir = d0 +'output'
    
    # os.system("cd %s && ~/projekty/TCLB/tools/sirun.sh d2q9_Allen_Cahn_SOI   ./run.xml >/dev/null"%d0)

    # os.system("cd ~/GITHUB/LBM/TCLB && CLB/d2q9_Allen_Cahn_SOI/main example/experimental/d2q9_Allen_Cahn_SOI.xml")
    os.system("cd %s && ~/GITHUB/LBM/TCLB/CLB/d2q9_Allen_Cahn_SOI/main ./run.xml >/dev/null"%d0)
    fname_base = "run_"    
    fconfig =  wdir + '/run_config_P00_00000000.xml'
    d = wdir
    if not os.path.isfile(fconfig):
        raise Exception("Not such case: " + fconfig)
         
        
    CLBc, CLBcf, CLBCn = CLB.CLBXMLHandler.parseConfig(fconfig,time=1E8)
    
    
    
    tmp = glob.glob(d + '/%sVTK_P00*.pvti'%fname_base)
    
    tmp = np.sort([ int(re.findall('[0-9]+',s)[-1]) for s in tmp])
    
    
    data = list()
    
    for f in tmp:
    
        fvti = d + '/%sVTK_P00_%08d.pvti'% (fname_base, f)
        vti = CLB.VTIFile.VTIFile(fvti, True)
        
        row = {}
        
        for fn in vti.getNames():
            row[fn] = np.average(vti.get(fn))
        row['Time'] = f*lbdt
    
        data.append(row)
        
    data = pd.DataFrame.from_records(data)
    
    ####
    from scipy.integrate import solve_ivp
    
    
    def AllenCahn(t, z, lambda_ph):
        return [lambda_ph*z*(1-z**2)]
    
    
    # CONSTANTS
    
    
    # INITIAL CONDItIONS
    IC = np.array([initial_phi])
    
    sol = solve_ivp(AllenCahn,
                    [0,data.Time.iloc[-1]],
                    IC,
                    method='RK45',
                    args=[lambda_ph],
                    dense_output=True)
    
    t = data.Time.to_numpy()
    z = sol.sol(t).T
    
    plt.figure()
    plt.title('RK45 / LBM for LBM_dt = %.1e'%lbdt)
    plt.plot(t,z, 'r-', lw=2, label='RK45')
    plt.plot(t,data.PhaseField, 'ko', lw=2, label='LBM')

    
    plt.xlabel('Time')
    plt.xlim(0, tc)
    plt.ylim(initial_phi, 1.01)
    plt.legend()
    plt.grid(which='both')
    plt.savefig( 'AC_LBM_0D__%.1e.png'%lbdt, dpi=200)
    
    L2.append(
        {
            'LBMdt': lbdt,
            'LBM': data.PhaseField,  
            'LBM_time': data.Time,
            'RK45': z, 
            'RK45_time': t
        }
      )
    print("dt done")





#L2 = pd.DataFrame.from_records(L2)

plt.figure()
reference = L2[-1]['LBM']

for i in range(len(L2)-1):
    t = L2[i]['LBM_time']
    # plt.loglog(L2[i]['LBMdt'],  np.sqrt(sint.trapz(( L2[i]['LBM'] - reference)**2, t)), 'ko')
    L2[i]['err'] = np.sqrt(sint.trapz(( L2[i]['LBM'] - reference)**2, t))


L2dr = pd.DataFrame.from_records(L2)

plt.loglog(L2dr.LBMdt, L2dr.err, 'ko', 'LBM')

dt = np.logspace(np.log10(L2[0]['LBMdt']), np.log10(L2[-1]['LBMdt']),100)
y = dt**2
y = y / y[0] * L2[0]['err']

plt.loglog(dt,y, label=r'${x^2}$')
plt.grid(which='both')
plt.xlabel('$\Delta t$')
plt.ylabel('$L_2(\phi(t,dt), \phi(t,dt_{min})$')

plt.legend()

plt.savefig( 'AC_LBM_0D_conv.png', dpi=200)
