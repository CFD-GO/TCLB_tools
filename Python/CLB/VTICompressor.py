# -*- coding: utf-8 -*-

import vtk
import sys
import os


from multiprocessing import Pool
NP = 1

save = sys.argv[1].split(",")

if save == [ sys.argv[1] ]:
    save = ['Rho', 'U', 'BOUNDARY']
    startWith = 1
else:
    startWith = 2

print save
    
def revti(name):
    in_vtifname = name
    out_vtifname = "/tmp/"+os.path.basename( os.path.splitext(name)[0] )+'.vti'
    
    reader = vtk.vtkXMLPImageDataReader()
    reader.SetFileName(in_vtifname)
    reader.Update()
    data = reader.GetOutput()
    cd = data.GetCellData()
    
    
    
    
    
    k = 0
    bc = cd.GetNumberOfArrays()
    while cd.GetNumberOfArrays() > len(save) or k > 10:
        k = k + 1
                
        for i in range(cd.GetNumberOfArrays()):
            if not ( cd.GetArrayName(i) in save ):
#                print "Removing:", cd.GetArrayName(i)
                cd.RemoveArray(cd.GetArrayName(i))
                break
#            else:
#                print "Keeping:", cd.GetArrayName(i)   
            
            
            
    print out_vtifname +" : ", bc, '/', cd.GetNumberOfArrays()
    
    
    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(out_vtifname)
    
    writer.SetInputData(data)
    writer.Write()


#pool = Pool(NP)
#pool.map( revti, sys.argv[startWith:])

for f in sys.argv[startWith:]:
    revti(f)