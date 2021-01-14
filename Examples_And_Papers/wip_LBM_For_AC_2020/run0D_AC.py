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

################################################
# CONFIG
initial_phi = 0.10  
lambda_ph0 = 1E-2
tc = 100
magic_parameter = 1./12. # to control even relaxation rate in TRT model
diffusivity0 = 1./6 #unimportant, no space variablitiy

################################################


# single branch solution: positive IC and Lambda
def AC0D(t, lambda_phi, phi_0):
    return np.sqrt(-1/(np.exp(-2*lambda_phi*t) - 1 - np.exp(-2*lambda_phi*t)/phi_0**2))


plot_dir = f'AC_plots_0D'
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    
L2 = list()
idx = 0
data = list()
for lbdt in 2.**-np.arange(1,11):
    
    lambda_ph = lambda_ph0*lbdt
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
        
        CLBc.addTRT_M_SOI()
        #CLBc.addSRT_M_SOI()
        #CLBc.addMRT()
        CLBc.addBox()
        
        # CLBc.addSmoothing()
        # CLBc.addBox()
        
        params = {
        		"diffusivity_phi": diffusivity0, 
        		"lambda":lambda_ph,
            "magic_parameter": magic_parameter,
        		"Init_PhaseField":initial_phi - lambda_ph*initial_phi*(1-initial_phi**2)/2 ,	
            #"Init_PhaseField":initial_phi,
        		#"phase_field_smoothing_coeff":0.0,
        }
        
        CLBc.addModelParams(params)
                 
  
        current = 0
        #for stop in np.logspace(0, np.log10(tc/lbdt), 100):    
        for stop in np.linspace(0, tc/lbdt, 101)[1:]:    
            
            CLBc.addSolve(iterations=stop-current)
            CLBc.addVTK()
            current = stop
        
        CLBc.addSolve(iterations=1)
        CLBc.addVTK()   
        #CLBc.addSolve(iterations=tc/lbdt, vtk=50)
        CLBc.write(fname+'.xml')


        return prefix
        
        
    
    d0 = getXML(clear=True)
    
    wdir = d0 +'output'
    
    os.system("cd %s && ~/projekty/TCLB/tools/sirun.sh d2q9_Allen_Cahn_SOI   ./run.xml >/dev/null"%d0)
    #os.system("cd %s && ~/projekty/TCLB/tools/sirun.sh d2q9_AllenCahn_BGK   ./run.xml >/dev/null "%d0)

    # os.system("cd ~/GITHUB/LBM/TCLB && CLB/d2q9_Allen_Cahn_SOI/main example/experimental/d2q9_Allen_Cahn_SOI.xml")
    #os.system("cd %s && ~/GITHUB/LBM/TCLB/CLB/d2q9_Allen_Cahn_SOI/main ./run.xml >/dev/null"%d0)
    fname_base = "run_"    
    fconfig =  wdir + '/run_config_P00_00000000.xml'
    d = wdir
    if not os.path.isfile(fconfig):
        raise Exception("Not such case: " + fconfig)
         
        
    CLBc, CLBcf, CLBCn = CLB.CLBXMLHandler.parseConfig(fconfig,time=1E8)
    
    
    
    tmp = glob.glob(d + '/%sVTK_P00*.pvti'%fname_base)
    
    tmp = np.sort([ int(re.findall('[0-9]+',s)[-1]) for s in tmp])
    
    
    
    
    for f in tmp:
    
        fvti = d + '/%sVTK_P00_%08d.pvti'% (fname_base, f)
        vti = CLB.VTIFile.VTIFile(fvti, True)
        
        row = {}
        
        for fn in vti.getNames():
            row[fn] = np.average(vti.get(fn))
        row['Time'] = f*lbdt
        row['lbdt_pow'] = np.log(lbdt)/np.log(2)
        data.append(row)
        

data = pd.DataFrame.from_records(data)
data['dt'] = 2**data.lbdt_pow

local_error = list()
for lbdt, ldata in data.groupby(by='lbdt_pow'):
    plt.plot(ldata.Time,ldata.PhaseField, 'o', lw=2, label='LBM dt=%d'%lbdt, markevery=10)
    
    err = np.abs(AC0D(ldata.dt.iloc[0], lambda_ph0, ldata.PhaseField.iloc[-2]) - ldata.PhaseField.iloc[-1])
    
    local_error.append({'dt': 2**lbdt, 'err': err})
    
local_error = pd.DataFrame.from_records(local_error)

plt.plot(ldata.Time, AC0D(ldata.Time, lambda_ph0, initial_phi)   ) 
plt.legend()


plt.figure()
plt.title('Local timestep error')
y2 = local_error.dt**3
y2 = y2 / local_error.dt.max()**3 * local_error.err.max()

plt.loglog(local_error.dt, local_error.err, 'o')
plt.loglog(local_error.dt,y2, label=r'$x^3$',base=10)

plt.legend()
plt.xlabel(r'$\Delta t$')
plt.ylabel(r'$\left|\phi(\Delta t) - \Phi(\Delta t)_{LBM} \right|$')
plt.grid(which='both')

plt.tight_layout()  # otherwise the right y-label is slightly clipped

plt.savefig('0d-local.pdf', bbox_inches='tight', dpi=200)
plt.savefig('0d-local.png', bbox_inches='tight', dpi=200) # for preview


data['PointError_sq'] = (AC0D(data.Time,  lambda_ph0, initial_phi) - data.PhaseField)**2
data['PointError'] = np.abs(AC0D(data.Time,  lambda_ph0, initial_phi) - data.PhaseField)

plt.figure()

for lbdt, ldata in data.groupby(by='lbdt_pow'):
    plt.plot(ldata.Time,ldata.PointError, 'o', lw=2, label='LBM dt=%d'%lbdt)
plt.legend()



L2 = data.groupby(by='lbdt_pow')[['PointError_sq']].sum().reset_index()
L2['dt'] = 2**L2.lbdt_pow
L2['L2'] = np.sqrt(ldata.PointError_sq)
 

# plt.figure()
y2 = L2.dt**2
y2 = y2 / 1E-2 * L2.L2.min()
# plt.loglog(L2.dt,y2)

# plt.loglog(L2.dt, L2.L2*np.sqrt(L2.dt), base=2)
# plt.grid(which='both')





plt.figure()

plt.title('Accumulated error')

ldata = data.query('Time == 100')
plt.loglog(ldata.dt,ldata.PointError, 'o-',label='T=100')

ldata = data.query('Time == 10')
plt.loglog(ldata.dt,ldata.PointError, 'o-',label='T=10')

ldata = data.query('Time == 1')
plt.loglog(ldata.dt,ldata.PointError, 'o-',label='T=1', base=2)

plt.grid(which='both')

plt.loglog(L2.dt,y2,label=r'$x^2$')



plt.legend()
plt.xlabel(r'$\Delta t$')
plt.ylabel(r'$\left|\phi(T) - \Phi(T)_{LBM} \right|$')


plt.legend()


plt.tight_layout()  # otherwise the right y-label is slightly clipped

plt.savefig('0d-accumulated.pdf', bbox_inches='tight', dpi=200)
plt.savefig('0d-accumulated.png', bbox_inches='tight', dpi=200) # for preview
# plt.show()


# plt.figure()

# plt.title('Self. conv at fixed time')

# ldata = data.query('Time == 100')

# ref = ldata.query('lbdt_pow == %f'%float(ldata.lbdt_pow.min()))

# plt.loglog(ldata.dt,np.absolute(ldata.PhaseField - float(ref.PhaseField)), 'o-', base=2)

# plt.loglog(L2.dt,y2)
# plt.grid(which='both')

