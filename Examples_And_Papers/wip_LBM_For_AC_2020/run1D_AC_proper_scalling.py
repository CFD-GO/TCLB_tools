#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 12:10:10 2020

@author: mdzik
"""

import CLB.CLBXMLWriter as CLBXML
import CLB.CLBXMLHandler

#assert(CLB.CLBXMLWriter.version >= 100)

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

def accoustic_scalling(dtn):
    return dtn

def diffusive_scalling(dtn):
    return dtn*2


################################################
# CONFIG
scalling = diffusive_scalling 

lambda_ph = 1E-3 # 1E-12 for pure diffusion
tc = 50 # number of timesteps for dt=1 aka Time
domain_size0=32
nsamples = 6 # number of resolutions
ny = 3 #numebr of nodes in second dimension, min 2 for CPU, min 3 for GPU
################################################

L2 = list()
idx = 0

for dtn in np.arange(0,nsamples)  : # (start, stop, num)
    lbdt = 1./(2**dtn)
    
    domain_size = domain_size0 * 2**scalling(dtn)
    
    lbdx = domain_size/domain_size0
    
    diffusivity = 1./6. * (lbdt/lbdx**2)
    
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
            CLBc.addGeomParam('ny', ny)
            
            
            CLBc.addTRT_M_SOI()
            CLBc.addBox()
            
            params = {
            		"diffusivity_phi": diffusivity,
                    "magic_parameter": 0.25,
            		"lambda":lambda_ph*lbdt,
            		"Init_PhaseField":-1 ,	
            		"phase_field_smoothing_coeff":0.0,
            }
            
            CLBc.addModelParams(params)

            CLBc.addRunR(eval=\
"""
		x = Solver$Geometry$X 
		x = (x -0.5)/ ({domain_size}) * 2 * pi
		y = Solver$Geometry$Y		
		Solver$Fields$Init_PhaseField_External[] = exp(sin(x))
		Solver$Actions$InitFromFields() 
        
""".format(domain_size=domain_size)
            
            )

            CLBc.addVTK()
            CLBc.addSolve(iterations=tc/lbdt-1)
            CLBc.addVTK()
            
            CLBc.write(fname+'.xml')
            return prefix
        
    
    d0 = getXML(clear=True)
    wdir = d0 + '/output'
    
    # os.system("cd %s && ~/projekty/TCLB/tools/sirun.sh d2q9_Allen_Cahn_SOI   ./run.xml"%d0)
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
        for field_name in ['PhaseField']:
            row[field_name] = vti.get(field_name)[0,::2**scalling(dtn)]
           # plt.plot(row[field_name])
           # plt.plot(-row[field_name][::-1], '-')

        row['Time'] = tmp*lbdt
    
        data.append(row)
        
    data = pd.DataFrame.from_records(data)
    
    L2.append(
        {            
            'dtN': dtn,
            'LBMdt': lbdt,
            'LBM_point': data.PhaseField,  
            'LBM_time': data.Time,
            'LBM_final': vti.get('PhaseField'),        
        }
      )
    
    
reference = L2[-1]['LBM_final'][:,::2**scalling(dtn)]
# plt.figure()
# for i in range(len(L2)-1):
#     t = L2[i]['LBM_time']
#     # plt.loglog(L2[i]['LBMdt'],  np.sqrt(sint.trapz(( L2[i]['LBM'] - reference)**2, t)), 'ko')
#     plt.plot(t, L2[i]['LBM_point'])

def calc_L2(anal, num):
    # Eq. 4.57
    return np.sqrt(np.sum((anal - num) * (anal - num)) / np.sum(anal * anal))    
    #return  np.max(np.abs(anal - num))

plot_dir = 'AC_plots'
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

plt.figure()
for i in range(len(L2)):
    # t = L2[i]['LBM_time']
    # plt.loglog(L2[i]['LBMdt'],  np.sqrt(sint.trapz(( L2[i]['LBM'] - reference)**2, t)), 'ko')
    # L2[i]['err'] = np.sqrt(sint.trapz(( L2[i]['LBM_point'] - reference)**2, t))

    final = L2[i]['LBM_final'][:,::2**scalling(L2[i]['dtN'])]
    
    L2[i]['err_field'] = np.sqrt(( final - reference)**2)
    L2[i]['err_L2'] = calc_L2(reference, final)

    plt.plot(np.linspace(0,1, L2[i]['LBM_final'].shape[1]), L2[i]['LBM_final'][0,:])
  #  plt.savefig(f'{plot_dir}/AC_LBM_2D_final_{L2[i][r"LBMdt"]:.1e}.png', dpi=300)
 #   plt.close(fig)

plt.grid(which='both')

plt.figure()
fig = plt.gcf()  # get current figure
L2dr = pd.DataFrame.from_records(L2[:-1])

plt.plot(L2dr.LBMdt, L2dr.err_L2, 'ko', 'LBM_point')

dt = np.logspace(np.log10(L2[0]['LBMdt']), np.log10(L2[-1]['LBMdt']),100)


y = dt**6
y = y / y[0] * L2[0]['err_L2']
plt.loglog(dt,y, label=r'${x^6}$')

y = dt**4
y = y / y[0] * L2[0]['err_L2']
plt.loglog(dt,y, label=r'${x^4}$')

y2 = dt**2
y2 = y2 / y2[0] * L2[0]['err_L2']
plt.loglog(dt,y2, label=r'${x^2}$')

plt.grid(which='both')
plt.xlabel('$\epsilon$')
plt.ylabel('$L_2(\phi(t,dt), \phi(t,dt_{min})$')

plt.legend()

plt.savefig('AC_LBM_2D_conv.png', dpi=200)
#plt.close(fig)

print("DONE.")

# #####################################################

# x = np.linspace(0,2*np.pi,1024)
# y0  = np.cos(x)

# import scipy.fftpack as sfp

# plt.plot(x,-np.sin(x))
# plt.plot(x,sfp.diff(y0,order=2))
