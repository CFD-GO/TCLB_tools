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
import math
import matplotlib.pyplot as plt

from numpy.testing import assert_almost_equal

def reactive_scalling(n):
    return 0

def acoustic_scalling(n):
    return n

def diffusive_scalling(n):
    return n*2


def eat_dots_for_texmaker(value):
    # s_value = str(value)
    s_value = f"{value:.2e}"
    s_value = re.sub(r"\.", 'o', s_value)
    return s_value

################################################
# CONFIG
tc = 64                     # number of timesteps for dt=1 aka Time
domain_size0 = 32           # initial size of the domain
# box_size0 = int(domain_size0/8)  # depreciated: prefer RunR # initial size of the box # 

# magic_parameter = 1./6     # best for pure diffusion   # to control even relaxation rate in TRT model
magic_parameter = 1./12      # best for pure advection   # to control even relaxation rate in TRT model

np.random.seed(seed=1)           # same random sequence each run
nrand = int(domain_size0/4)      # number of boxes, doesn't converge nice for nrand 100

initialDa = 1E3
initialPe = 1E3

# diffusivity0 = 1./6. * 2* 1E-5  # initial diffusivity
# lambda_ph0 = 1E-3               # worse than 2nd order, set 1E-12 for pure diffusion
# nsamples = 6                    # number of resolutions

# scalling = reactive_scalling 
# scalling = acoustic_scalling 
scalling = diffusive_scalling 


diffusivity0 = 1./6. * 1E-2  # initial diffusivity
lambda_ph0 = initialDa*diffusivity0/domain_size0**2 
nsamples = 5                 # number of resolutions

Ux0 = initialPe*diffusivity0/domain_size0    # macroscopic advection velocity


# check Da, Pe
Da0 = (lambda_ph0 *  domain_size0**2) / diffusivity0  # Damkohler similarity number
Pe0 = Ux0*domain_size0 / diffusivity0  # Peclet similarity number (similar to Reynolds, but refers to advection-diffusion eq)
print(f"Initial Damköhler number: {Da0:.2e} \t Peclet number: {Pe0:.2e}")

dfObj = pd.DataFrame()
################################################

# blocks_up0 = list()
# blocks_down0 = list()
# for xi in range(nrand):
#     x = np.random.randint(0,domain_size0)
#     y = np.random.randint(0,domain_size0)
#     blocks_up0.append((x,y))
 

# for xi in range(nrand):
#     x = np.random.randint(0,domain_size0)
#     y = np.random.randint(0,domain_size0)
#     blocks_down0.append((x,y))

idx = 0
L2 = list()

for n in np.arange(0,nsamples): # (start, stop, step=1)
    
    domain_size = domain_size0 * 2**n

    # begin of acoustic digression  
    # todo: for acoustic scaling dt = (U_LB/U_si) * dx  -> or dx/dt = const... = e?
    # we are going to stay on the same characteristic... but e used to be 1?
    coeff = 1
     
    # there are infinitely many ways to match Da, let
    # lbdt = coeff*1./(2**scalling(n))
    # lbdx = 1./2**n
    # diffusivity = diffusivity0 * lbdt / lbdx**2
    # lambda_ph = lambda_ph0 * lbdt
    
    # alternatively
    lbdt = 1./(2**scalling(n))
    lbdx = 1./2**n
    diffusivity = coeff * diffusivity0 * lbdt / lbdx**2
    lambda_ph = coeff *  lambda_ph0 * lbdt
    # end of acoustic digression
    
    Da = (lambda_ph *  domain_size**2) / diffusivity  
    
    Ux = Pe0*diffusivity/domain_size
    
    # box_size = box_size0 * 2**n
    # blocks_up   = [ (point[0] * 2**n, point[1] * 2**n)  for point in blocks_up0]
    # blocks_down = [ (point[0] * 2**n, point[1] * 2**n)  for point in blocks_down0]
    
    n_iterations = tc/lbdt
    print(f"running case {n}/{nsamples} Da = {Da:.2e}; lbdt={lbdt}, lbdx={lbdx}, Ux={Ux:.2e} diffusivity={diffusivity:.2e} lambda_ph={lambda_ph:.2e} ")
    assert_almost_equal(Da, Da0, decimal=4)
    assert_almost_equal(Ux*domain_size/diffusivity, Pe0, decimal=4)
    assert_almost_equal(math.modf(n_iterations)[0] , 0, decimal=4) # check decimal places
    
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
        CLBc.addGeomParam('ny', domain_size)
        
        # CLBc.addSmoothing()
        CLBc.addTRT_M_SOI()
        # CLBc.addSRT_M_SOI()
        CLBc.addBox()
        
        # CLBc.addZoneBlock(name='up')
        # for x,y in blocks_up:
        #     CLBc.addBox(dx=x,dy=y,nx=box_size,ny=box_size)
 
        # CLBc.addZoneBlock(name='down')
        # for x,y in blocks_down:
        #     CLBc.addBox(dx=x,dy=y,nx=box_size,ny=box_size)


        CLBc.addModelParams({
                "diffusivity_phi":diffusivity,
                "magic_parameter": magic_parameter,
                "lambda":lambda_ph,
                "Init_PhaseField": 0,
                "Init_UX" : Ux
                })

        # CLBc.addModelParam("Init_PhaseField",0.9,'up' )
        # CLBc.addModelParam("Init_PhaseField",-0.9,'down' )
                    
        CLBc.addRunR(eval=\
            """
                    x = Solver$Geometry$X 
                    x = (x - 0.5)/ ({domain_size}) * 2 * pi
                    y = Solver$Geometry$Y
                    y = (y - 0.5)/ ({domain_size}) * 4 * pi
                    Solver$Fields$Init_PhaseField_External[] = exp(sin(x)) - 2*exp(sin(y)) # to benchmark diffusion & source term
                    # Solver$Fields$Init_PhaseField_External[] = 0.1 # to benchmark the source term only
                    Solver$Actions$InitFromFields()        
            """.format(domain_size=domain_size))
            
        #CLBc.addSolve(iterations=n_iterations, vtk=50)
        CLBc.addVTK()
        CLBc.addSolve(iterations=n_iterations, vtk=int(n_iterations/10))
        CLBc.addVTK()
        
        CLBc.write(fname+'.xml')
        return prefix
    
    
    d0 = getXML(clear=True)
    wdir = d0 + '/output'
    
    # os.system("cd %s && ~/projekty/TCLB/tools/sirun.sh d2q9_Allen_Cahn_SOI   ./run.xml >/dev/null"%d0)
    os.system(f"cd {d0} && ~/GITHUB/LBM/TCLB/CLB/d2q9_Allen_Cahn_SOI/main ./run.xml > log.txt")
    
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
            row[field_name] = vti.get(field_name)[::2**n,::2**n] # get each n-th point according to scalling
        
        row['Time'] = tmp*lbdt
    
        data.append(row)
        
    data = pd.DataFrame.from_records(data)
    
    
    L2.append(
        {
            'n': n,
            'LBMdx': lbdx,
            'LBMdt': lbdt,
            'SI_time': data.Time,
            'LBM_field_slice': data.PhaseField,        # slices from all iterations
            'LBM_field_all': vti.get('PhaseField'),    # not sliced, last iteration 
        }
      )
    
    tmp_df = pd.DataFrame({
           'L': domain_size,
           'n': int(n_iterations),
           # 'log2(LBMdx)': np.log2(lbdx),
           # 'log2(LBMdt)': np.log2(lbdt),
           r'$time_{SI}$': int(data.Time.iloc[-1]),
           r'$\Lambda$': magic_parameter,
           # 'Da': int(Da),
           # 'Pe': int(Pe0),
           r'U': Ux,
           'M':diffusivity,
           r'$\lambda$':lambda_ph,
         }, index=[idx]) 
    
    dfObj = dfObj.append(tmp_df) 
    
reference = L2[-1]['LBM_field_all'][::2**n,::2**n]


def calc_L2(anal, num):
    # Eq. 4.57
    return np.sqrt(np.sum((anal - num) * (anal - num)) / np.sum(anal * anal))    
    #return  np.max(np.abs(anal - num))
    
    
plot_dir = f'AC_plots_2D_{scalling.__name__}'
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)


for i in range(len(L2)-1):
    final = L2[i]['LBM_field_all'][::2**L2[i]['n'],::2**L2[i]['n']]

    L2[i]['err_field'] = np.sqrt((final - reference)**2)
    L2[i]['err_L2'] = calc_L2(reference, final)

    # plt.figure()
    # fig = plt.gcf()  # get current figure
    # plt.imshow(L2[i]['LBM_field_all'])
    # plt.savefig(f'{plot_dir}/AC_LBM_2D_field_all_{i}_{L2[i][r"LBMdt"]:.1e}_diffusivity0_{diffusivity0:.2e}_lambda_ph0_{lambda_ph0:.2e}.png', dpi=300)
    # plt.close(fig)

    # plt.figure()
    # fig = plt.gcf()  # get current figure
    # plt.imshow(L2[i]['err_field'])
    # plt.savefig(f'{plot_dir}/AC_LBM_2D_err_field_{i}_{L2[i][r"LBMdt"]:.1e}_diffusivity0_{diffusivity0:.2e}_lambda_ph0_{lambda_ph0:.2e}.png', dpi=300)
    # plt.close(fig)
    # plt.show()
    


##########################

# fig_name = os.path.join(plot_dir, f"{scalling.__name__}_2D_conv_Da_{eat_dots_for_texmaker(Da)}_diffusivity0_{diffusivity0:.2e}_lambda_ph0_{lambda_ph0:.2e}")


fig_name = os.path.join(plot_dir, 
                        f"{scalling.__name__}"
                        f"_2D_conv_Da_{eat_dots_for_texmaker(Da)}"
                        f"_diffusivity0_{eat_dots_for_texmaker(diffusivity)}"
                        f"_lambda_ph0_{eat_dots_for_texmaker(lambda_ph0)}")

plt.rcParams.update({'font.size': 26})


fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))    
# fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(16, 10))

ax1 = axs
# ax2 = axs[1]

# prepare data
L2dr = pd.DataFrame.from_records(L2[:-1])

dt = np.logspace(np.log10(L2[0]['LBMdt']), np.log10(L2[-1]['LBMdt']),100)
dx = np.logspace(np.log10(L2[0]['LBMdx']), np.log10(L2[-1]['LBMdx']),100)

# y = np.sqrt(dx)
# y = y / y[0] * L2[0]['err_L2']
# ax1.loglog(dx,y, label=r'${x^{1/2}}$', linewidth=2)

y = dx**1
y = y / y[0] * L2[0]['err_L2']
ax1.loglog(dx,y, label=r'${x}$', linewidth=2)

y = dx**2
y = y / y[0] * L2[0]['err_L2']
ax1.loglog(dx,y, label=r'${x^2}$', linewidth=2)



ax1.loglog(L2dr.LBMdx, L2dr.err_L2, linestyle="", color='black', marker='x', markersize=12)
ax1.set_xscale('log', base=2)
ax1.legend()
ax1.grid(True)
ax1.set(xlabel=r'$\epsilon_x$', ylabel=r'$L_2(\phi(\delta x), \phi(\delta x_{min})$')

    

# y = np.sqrt(dt)
# y = y / y[0] * L2[0]['err_L2']
# ax2.loglog(dx,y, label=r'${t^{1/2}}$')

# y = dt**1
# y = y / y[0] * L2[0]['err_L2']
# ax2.loglog(dx,y, label=r'${t}$')

# y = dt**2
# y = y / y[0] * L2[0]['err_L2']
# ax2.loglog(dx,y, label=r'${t^2}$')



# ax2.loglog(L2dr.LBMdt, L2dr.err_L2,  linestyle="", color='black', marker='o', markersize=10)
# ax2.set_xscale('log', base=2)
# ax2.legend()
# ax2.grid(True)
# ax2.set(xlabel=r'$\epsilon_t$', ylabel=r'$L_2(\phi(\delta t), \phi(\delta t_{min})$')

fig.tight_layout()  # otherwise the right y-label is slightly clipped

plt.pause(1e-9)  # there is a race condition somewhere in the matplotlib code.
fig.savefig(fig_name + '.png', bbox_inches='tight', dpi=200)
# plt.show()
plt.close(fig)  # close the figure




print(f"saved: {fig_name}")

print(f"SUMMARY:")
print(dfObj.to_latex(index=False, escape=False))

with pd.ExcelWriter(f"{fig_name}.xlsx") as writer:  # doctest: +SKIP
    # dfObj = dfObj.sort_values(by=['is3D', 'Collision_Kernel', 'D', 'BC_order'])
    dfObj.to_excel(writer, sheet_name='EnhancedTable', index=False)
    # legend_df.to_excel(writer, sheet_name='EnhancedTable', startrow=len(new_df) + 2, index=False)


print(f"DONE.")
# According to chapter 10.7.2, eq 10.48, p429 from 'The Lattice Boltzmann Method: Principles and Practice'
# by T. Krüger, H. Kusumaatmaja, A. Kuzmin, O. Shardt, G. Silva, E.M. Viggen
# There are certain values of magic_parameter that show distinctive properties:
# • magic_parameter 1./12 = 0.08(3) cancels the third-order spatial error, leading to optimal results for pure advection problems.
# • magic_parameter 1./6 = 0.1(6) cancels the fourth-order spatial error, providing the most accurate results for the pure diffusion equation.
# • magic_parameter 3./16 = 0.1875 results in the boundary wall location implemented via bounce-back for the Poiseuille flow exactly in the middle between horizontal walls and fluid nodes.
# • magic_parameter 1./4 = 0.25 provides the most stable simulations. 
