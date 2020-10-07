#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 12:10:10 2020

@author: mdzik
"""

import CLBXMLWriter as CLBXML
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
lambda_ph = 1E-14

L2 = list()

domain_size=32

blocks_up = list()
blocks_down = list()

np.random.seed(seed=1)

# nrand = 1 # doesn't converge nice for nrand 100
# for xi in range(nrand):
#     x = np.random.randint(0,domain_size)
#     y = np.random.randint(0,domain_size)
#     blocks_up.append((x,y))
#     # CLBc.addSmoothing()


# for xi in range(nrand):
#     x = np.random.randint(0,domain_size)
#     y = np.random.randint(0,domain_size)
#     blocks_down.append((x,y))

# for lbdt in 2.**-np.arange(1,8):
# for lbdt in 1./np.linspace(1, 10, 10): # (start, stop, num)

# x1 = 1./np.logspace(0, 1, 10, base=10)   
# x2 = 1./np.linspace(1, 10, 10)   

# all_together = np.concatenate([x1,x2]) 
# for lbdt in all_together: # (start, stop, num)

for dtn in np.arange(1,9)  : # (start, stop, num)
    lbdt = 1./(2**dtn)

    tc = 400
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
            CLBc.addGeomParam('ny', 2)
            
            
            CLBc.addTRT_SOI()
            CLBc.addBox()
            # CLBc.addSRT_SOI()

            
            params = {
            		"diffusivity_phi": 1./6.*lbdt,
                "magic_parameter": 1./6.,
            		"lambda":10E-15,#lambda_ph*lbdt,
            		"Init_PhaseField":-1 ,	
            		"phase_field_smoothing_coeff":0.0,
            }
            
            CLBc.addModelParams(params)

           # CLBc.addModelParam("Init_PhaseField",1,'up' )
           # CLBc.addModelParam("Init_PhaseField",-0.9,'down' )
                     
            current = 0
            #for stop in np.logspace(0, np.log10(tc/lbdt), 100):    
            # for stop in np.linspace(1, tc/lbdt, 51): # watch out for fractions    
            #     CLBc.addSolve(iterations=stop-current)
            #     CLBc.addVTK()
            #     current = stop

            
            #CLBc.addSolve(iterations=tc/lbdt, vtk=50)
            CLBc.addRunR(eval=\
"""
		x = Solver$Geometry$X 
		x = x / ({domain_size}) * 2 * pi
		y = Solver$Geometry$Y		
		Solver$Fields$Init_PhaseField_External[] = sin(x)
		Solver$Actions$InitFromFields() 
        
""".format(domain_size=domain_size)
            )
            CLBc.addVTK()
            CLBc.addSolve(iterations=tc/lbdt, vtk=tc/lbdt/4)
            
            CLBc.write(fname+'.xml')
            return prefix
        
    
    d0 = getXML(clear=True)
    wdir = d0 + '/output'
    
    os.system("cd %s && ~/projekty/TCLB/tools/sirun.sh d2q9_Allen_Cahn_SOI   ./run.xml"%d0)
    #os.system(f"cd {d0} && ~/GITHUB/LBM/TCLB/CLB/d2q9_Allen_Cahn_SOI/main ./run.xml >/dev/null")
    
    fname_base = "run_"    
    fconfig =  wdir + '/run_config_P00_00000000.xml'
    d = wdir
    if not os.path.isfile(fconfig):
        raise Exception("Not such case: " + fconfig)
         
        
    CLBc, CLBcf, CLBCn = CLB.CLBXMLHandler.parseConfig(fconfig,time=1E8)

    tmps = glob.glob(wdir + '/%sVTK_P00*.pvti'%fname_base)
    tmps = np.sort([int(re.findall('[0-9]+',s)[-1]) for s in tmps])
    
    data = list()
    x = (np.arange(0,domain_size) + 0.5)  / float(domain_size)*2*np.pi
    C1 = - 4. * np.pi**2 * 1./6. / float(domain_size)**2
    
    for tmp in tmps:
        fvti = wdir + '/%sVTK_P00_%08d.pvti'% (fname_base, tmp)
        vti = CLB.VTIFile.VTIFile(fvti, True)
        
        row = {}
        for field_name in ['PhaseField']:
            row[field_name] = vti.get(field_name)[0,:]
            phi = row[field_name]
         #   plt.plot(row[field_name]+row[field_name][::-1], '-')
            
        row['Time'] = (tmp)*lbdt
        row['TimeFromSinExp'] = np.average(np.log(phi/np.sin(x)) / C1)
        print( row['Time'],  row['TimeFromSinExp'] )
        data.append(row)
        sdfsfd

    #plt.plot(x,np.exp(C1*tc)*np.sin(x/np.max(x)*2*np.pi), 'k')
    plt.plot(x,phi-np.exp(C1*tc)*np.sin(x),'-')
    
    
    
    data = pd.DataFrame.from_records(data)
    
    L2.append(
        {
            'LBMdt': lbdt,
            'LBM_point': data.PhaseField,  
            'LBM_time': data.Time,
            'LBM_final': vti.get('PhaseField'),        
        }
      )
    
    
    
    
x = (np.arange(0,domain_size) + 0.5)  / float(domain_size)*2*np.pi
C1 = - 4. * np.pi**2 * 1./6. / float(domain_size)**2

#reference = np.exp(C1*tc)*np.sin(x)

#plt.plot(reference-final)

reference = L2[-1]['LBM_final'][0,:]
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

    final = L2[i]['LBM_final'][0,:]
    
    L2[i]['err_field'] = np.sqrt(( final - reference)**2)
    L2[i]['err_L2'] = calc_L2(reference, final)

    plt.plot(final-reference, '-')
    plt.plot(reference, 'k')
  #  plt.savefig(f'{plot_dir}/AC_LBM_2D_final_{L2[i][r"LBMdt"]:.1e}.png', dpi=300)
 #   plt.close(fig)



plt.figure()
fig = plt.gcf()  # get current figure
L2dr = pd.DataFrame.from_records(L2[:-1])

plt.plot(L2dr.LBMdt, L2dr.err_L2, 'ko', 'LBM_point')

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
#plt.close(fig)



# #####################################################

# x = np.linspace(0,2*np.pi,1024)
# y0  = np.cos(x)

# import scipy.fftpack as sfp

# plt.plot(x,-np.sin(x))
# plt.plot(x,sfp.diff(y0,order=2))
