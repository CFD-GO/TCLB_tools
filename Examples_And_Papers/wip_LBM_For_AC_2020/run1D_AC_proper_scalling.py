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

def accoustic_scalling(n):
    return n

def diffusive_scalling(n):
    return n*2


################################################
# CONFIG
scalling = accoustic_scalling 
# scalling = diffusive_scalling 
diffusivity0 = 1./6.    # initial diffusivity
lambda_ph0 = 1E-10       # 1E-12 for pure diffusion, source term strength
domain_size0 = 32      # initial size of the domain
# domain_size0 = 16 # TODO nie dziala dla 16 operands could not be broadcast together with shapes (3,32) (3,16)
tc = 50                 # number of timesteps for dt=1 aka Time
nsamples = 6            # number of resolutions
ny = 3                  # number of nodes in second dimension, min 2 for CPU, min 3 for GPU
magic_parameter = 0.25  # to control even relaxation rate in TRT model

Da0 = (lambda_ph0 *  domain_size0**2) / diffusivity0  # Damkohler similarity number
print(f"Initial Damköhler number: {Da0:.2e}")
################################################

L2 = list()
idx = 0
for n in np.arange(0,nsamples): # (start, stop, step=1)

    # MD version
    lbdt = 1./(2**n)
    domain_size = domain_size0 * 2**scalling(n)
    lbdx = domain_size/domain_size0
    diffusivity = diffusivity0 * (lbdt/lbdx**2)
    lambda_ph = lambda_ph0*lbdt
    # # end od MD version
    
    # domain_size = domain_size0 * 2**n
    # lbdt = 1./(2**scalling(n))
    # lbdx = 1./2**n
    # diffusivity = diffusivity0 * lbdt / lbdx**2
    # lambda_ph = lambda_ph0 * lbdt

    Da = (lambda_ph *  domain_size**2) / diffusivity  
    print(f"running case {n}/{nsamples} Da = {Da:.2e} ; lbdt={lbdt}, lbdx={lbdx},  diffusivity={diffusivity:.2e} lambda_ph={lambda_ph:.2e}")

    def getXML(**kwars):

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
        # CLBc.addSRT_M_SOI()
        CLBc.addBox()
        
        params = {
                "diffusivity_phi": diffusivity,
                "magic_parameter": magic_parameter,
                "lambda": lambda_ph,
                "Init_PhaseField": -1,
                "phase_field_smoothing_coeff": 0.0,
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
        CLBc.addSolve(iterations=tc/lbdt-1) # why -1?
        CLBc.addVTK()
        CLBc.write(fname+'.xml')
        return prefix
        
    
    d0 = getXML(clear=True)
    wdir = d0 + '/output'
    
    # os.system("cd %s && ~/projekty/TCLB/tools/sirun.sh d2q9_Allen_Cahn_SOI   ./run.xml"%d0)
    os.system(f"cd {d0} && ~/GITHUB/LBM/TCLB/CLB/d2q9_SourceTerm_SOI_AllenCahn/main ./run.xml >/dev/null")
    
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
            # row[field_name] = vti.get(field_name)[50,50] # get point like (50, 50)
            row[field_name] = vti.get(field_name)[:,::2**scalling(n)]  # get each n-th point according to scalling

        row['Time'] = tmp*lbdt
        data.append(row)
        
    data = pd.DataFrame.from_records(data)
    
    L2.append(
        {            
            'n': n,
            'LBMdt': lbdt,
            'LBM_time': data.Time,
            'LBM_field_slice': data.PhaseField, # slices from all iterations
            'LBM_field_all': vti.get('PhaseField'), # not sliced, last iteration   
        }
      )
    # xxxreference0 = data.PhaseField[0]
    # xxxreference1 = data.PhaseField[1]
    # xxxreference2 = L2[-1]['LBM_field_slice'][0]
    # xxxreference3 = L2[-1]['LBM_field_slice'][1]
    # xxxreference4 = L2[-1]['LBM_field_all'][:,::2**scalling(n)]
    # xxxreference5 = vti.get(field_name)[:,::2**scalling(n)]


reference = L2[-1]['LBM_field_all'][:,::2**scalling(n)]

# plt.figure()
# for i in range(len(L2)-1):
#     t = L2[i]['LBM_time']
#     # plt.loglog(L2[i]['LBMdt'],  np.sqrt(sint.trapz(( L2[i]['LBM'] - reference)**2, t)), 'ko')
#     plt.plot(t, L2[i]['LBM_point'])

def calc_L2(anal, num):
    # Eq. 4.57
    return np.sqrt(np.sum((anal - num) * (anal - num)) / np.sum(anal * anal))    
    #return  np.max(np.abs(anal - num))

plot_dir = 'AC_plots_1D'
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

plt.figure()
for i in range(len(L2)):
    # t = L2[i]['LBM_time']
    # plt.loglog(L2[i]['LBMdt'],  np.sqrt(sint.trapz(( L2[i]['LBM'] - reference)**2, t)), 'ko')
    # L2[i]['err'] = np.sqrt(sint.trapz(( L2[i]['LBM_point'] - reference)**2, t))

    final = L2[i]['LBM_field_all'][:,::2**scalling(L2[i]['n'])]

    L2[i]['err_field'] = np.sqrt(( final - reference)**2)
    L2[i]['err_L2'] = calc_L2(reference, final)

    plt.plot(np.linspace(0,1, L2[i]['LBM_field_all'].shape[1]), L2[i]['LBM_field_all'][0,:])
    plt.savefig(f'{plot_dir}/AC_LBM_1D_LBM_field_all_{L2[i][r"LBMdt"]:.1e}.png', dpi=300)
     
fig = plt.gcf()  # get current figure
plt.close(fig)

plt.grid(which='both')
plt.figure()
fig = plt.gcf()  # get current figure

L2dr = pd.DataFrame.from_records(L2[:-1])

plt.plot(L2dr.LBMdt, L2dr.err_L2, 'ko')


dt = np.logspace(np.log10(L2[0]['LBMdt']), np.log10(L2[-1]['LBMdt']),100)

#y = dt**1
#y = y / y[0] * L2[0]['err_L2']
#plt.loglog(dt,y, label=r'${x}$')

y = dt**2
y = y / y[0] * L2[0]['err_L2']
plt.loglog(dt,y, label=r'${x^2}$')

y = dt**4
y = y / y[0] * L2[0]['err_L2']
plt.loglog(dt,y, label=r'${x^4}$')

y = dt**6
y = y / y[0] * L2[0]['err_L2']
plt.loglog(dt,y, label=r'${x^6}$')



plt.grid(which='both')
plt.xlabel('$\epsilon$')
plt.ylabel('$L_2(\phi(t,dt), \phi(t,dt_{min})$')
plt.legend()


fig_name = f"AC_LBM_1D_conv_{scalling.__name__}_diffusivity0_{diffusivity0:.2e}_lambda_ph0_{lambda_ph0:.2e}.png"
plt.savefig(os.path.join(plot_dir,fig_name), dpi=200)
plt.close(fig)

print(f"DONE. saved: {fig_name}")


# #####################################################

# x = np.linspace(0,2*np.pi,1024)
# y0  = np.cos(x)

# import scipy.fftpack as sfp

# plt.plot(x,-np.sin(x))
# plt.plot(x,sfp.diff(y0,order=2))
