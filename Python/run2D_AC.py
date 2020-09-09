#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 12:10:10 2020

@author: mdzik
"""

import CLB.CLBXMLWriter as CLBXML
import CLB.CLBXMLHandler
import CLB.VTIFile
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

domain_size=128
box_size=50

blocks_up = list()
blocks_down = list()

np.random.seed(seed=1)

nrand = 1
# for xi in range(nrand):
#     x = np.random.randint(0,domain_size)
#     y = np.random.randint(0,domain_size)
#     blocks_up.append((x,y))
#     # CLBc.addSmoothing()


for xi in range(nrand):
    x = np.random.randint(0,domain_size)
    y = np.random.randint(0,domain_size)
    blocks_down.append((x,y))


# for lbdt in 2.**-np.arange(1,8):
# for lbdt in 1./np.linspace(1, 10, 10): # (start, stop, num)
# for lbdt in  1./np.logspace(0, 1, 10, base=10)   : # (start, stop, num)

for lbdt in 1./np.logspace(0, 3, 20, base=10): # (start, stop, num)
    tc = 100
    def getXML(**kwars):
            print(f"running case: lbdt={lbdt}")
            global idx
            
            idx = idx + 1
            prefix = '/tmp/id%03d/'%idx
            if 'clear' in kwars and kwars['clear']:
                print(f"removing {prefix}")
                os.system('rm -r %s'%prefix)
    
            os.system('mkdir %s'%prefix)
        
            CLBc = CLBXML.CLBConfigWriter( )
            fname = prefix+"run"
            CLBc.addGeomParam('nx', domain_size)
            CLBc.addGeomParam('ny', domain_size)
            
            
            CLBc.addTRT_SOI()
            # CLBc.addSRT_SOI()
            CLBc.addBox()
            
            CLBc.addZoneBlock(name='up')
           
            for x,y in blocks_up:
                CLBc.addBox(dx=x,dy=y,nx=box_size,ny=box_size)
                # CLBc.addSmoothing()
            
            CLBc.addZoneBlock(name='down')

            for x,y in blocks_down:
                CLBc.addBox(dx=x,dy=y,nx=box_size,ny=box_size)
        
            params = {
            		"diffusivity_phi":0.1666666*lbdt,
                    "magic_parameter": 0.25,
            		"lambda":lambda_ph*lbdt,
            		"Init_PhaseField":0 ,	
            		"phase_field_smoothing_coeff":0.1,
            }
            
            CLBc.addModelParams(params)

            CLBc.addModelParam("Init_PhaseField",0.9,'up' )
            CLBc.addModelParam("Init_PhaseField",-0.9,'down' )
                     
            current = 0
            #for stop in np.logspace(0, np.log10(tc/lbdt), 100):    
            # for stop in np.linspace(1, tc/lbdt, 101): # watch out for fractions    
            #     CLBc.addSolve(iterations=stop-current)
            #     CLBc.addVTK()
            #     current = stop

            
            #CLBc.addSolve(iterations=tc/lbdt, vtk=50)

            CLBc.addSolve(iterations=tc/lbdt)
            CLBc.addVTK()
            
            CLBc.write(fname+'.xml')
            return prefix
        
    
    d0 = getXML(clear=True)
    wdir = d0 + '/output'
    
    # os.system("cd %s && ~/projekty/TCLB/tools/sirun.sh d2q9_Allen_Cahn_SOI   ./run.xml >/dev/null"%d0)
    os.system(f"cd {d0} && ~/GITHUB/LBM/TCLB/CLB/d2q9_Allen_Cahn_SOI/main ./run.xml >/dev/null")
    
    fname_base = "run_"    
    fconfig =  wdir + '/run_config_P00_00000000.xml'
    d = wdir
    if not os.path.isfile(fconfig):
        raise Exception("Not such case: " + fconfig)
         
        
    CLBc, CLBcf, CLBCn = CLB.CLBXMLHandler.parseConfig(fconfig,time=1E8)

    tmps = glob.glob(wdir + '/%sVTK_P00*.pvti'%fname_base)
    tmps = np.sort([int(re.findall('[0-9]+',s)[-1]) for s in tmps])
    
    data = list()
    for tmp in tmps:
        fvti = wdir + '/%sVTK_P00_%08d.pvti'% (fname_base, tmp)
        vti = CLB.VTIFile.VTIFile(fvti, True)
        
        row = {}
        for field_name in vti.getNames():
            row[field_name] = vti.get(field_name)[50,50]
        row['Time'] = tmp*lbdt
    
        data.append(row)
        
    data = pd.DataFrame.from_records(data)
    
    
    L2.append(
        {
            'LBMdt': lbdt,
            'LBM_point': data.PhaseField,  
            'LBM_time': data.Time,
            'LBM_final': vti.get('PhaseField'),        
        }
      )
    
    
reference = L2[-1]['LBM_final']
# plt.figure()
# for i in range(len(L2)-1):
#     t = L2[i]['LBM_time']
#     # plt.loglog(L2[i]['LBMdt'],  np.sqrt(sint.trapz(( L2[i]['LBM'] - reference)**2, t)), 'ko')
#     plt.plot(t, L2[i]['LBM_point'])

def calc_L2(anal, num):
    # Eq. 4.57
    return np.sqrt(np.sum((anal - num) * (anal - num)) / np.sum(anal * anal))    

plot_dir = 'AC_plots'
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

for i in range(len(L2)-1):
    # t = L2[i]['LBM_time']
    # plt.loglog(L2[i]['LBMdt'],  np.sqrt(sint.trapz(( L2[i]['LBM'] - reference)**2, t)), 'ko')
    # L2[i]['err'] = np.sqrt(sint.trapz(( L2[i]['LBM_point'] - reference)**2, t))

    L2[i]['err_field'] = np.sqrt((L2[i]['LBM_final'] - reference)**2)
    L2[i]['err_L2'] = calc_L2(reference, L2[i]['LBM_final'])

    plt.figure()
    plt.imshow(L2[i]['LBM_final'])
    plt.savefig(f'{plot_dir}/AC_LBM_2D_final_{L2[i][r"LBMdt"]:.1e}.png', dpi=300)

    plt.figure()
    plt.imshow(L2[i]['err_field'])
    plt.savefig(f'{plot_dir}/AC_LBM_2D_err_field_{L2[i][r"LBMdt"]:.1e}.png', dpi=300)

    # plt.show()
    

plt.figure()

L2dr = pd.DataFrame.from_records(L2)

plt.loglog(L2dr.LBMdt, L2dr.err_L2, 'ko', 'LBM_point')

dt = np.logspace(np.log10(L2[0]['LBMdt']), np.log10(L2[-1]['LBMdt']),100)

y = dt**1
y = y / y[0] * L2[0]['err_L2']
plt.loglog(dt,y, label=r'${x}$')

y2 = dt**2
y2 = y2 / y2[0] * L2[0]['err_L2']
plt.loglog(dt,y2, label=r'${x^2}$')

plt.grid(which='both')
plt.xlabel('$\Delta t$')
plt.ylabel('$L_2(\phi(t,dt), \phi(t,dt_{min})$')

plt.legend()

plt.savefig('AC_LBM_2D_conv.png', dpi=200)
