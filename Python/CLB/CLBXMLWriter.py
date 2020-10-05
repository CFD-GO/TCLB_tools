# -*- coding: utf-8 -*-

import lxml.etree as ET

## DECORATORS ##

def geometryElement(func):
    def fin(self, *args, **kwargs):
        nargs = func(self, *args, **kwargs)

        if len(nargs) < 1 or \
        not  '_xml_node_name' in nargs  \
        :
            raise BaseException("NOT ENOUGHT ARGS FOR GEOM ELEMENT")

        n = ET.SubElement(self.current_geometry,nargs['_xml_node_name'])
        for k in nargs:
            if not k == '_xml_node_name':
                n.set(k, str(kwargs[k]))
        return n
    return fin

def rootElement(func):
    def fin(self, *args, **kwargs):
        nargs = func(self, *args, **kwargs)

        if len(nargs) < 1 or \
        not '_xml_node_name' in nargs.keys() \
        :
            raise BaseException("NOT ENOUGHT ARGS FOR ROOT ELEMENT")

        if 'parent' in nargs.keys():
            n = ET.SubElement(nargs['parent'],nargs['_xml_node_name'])
        else:
            n = ET.SubElement(self.root,nargs['_xml_node_name'])
        for k in nargs:
            if not k in ('_xml_node_name', 'parent'):
                n.set(k, str(kwargs[k]))
        return n
    return fin


def addCDATA(name):
    def nif(func):
        def fin(self, *args, **kwargs):
            n = func(self, *args, **kwargs)
            n.text = ET.CDATA(kwargs['eval'])
            #n.append(cdata)
        return fin
    return nif


def BCElement(func):
    def fin(self, *args, **kwargs):
        nargs = func(self, *args, **kwargs)

        if len(nargs) < 1 or \
        not '_xml_node_name' in nargs \
        :
            raise BaseException("NOT ENOUGHT ARGS FOR BC ELEMENT")

        n = ET.SubElement(self.geometry,nargs['_xml_node_name'])
        for k in nargs:
            if not k == '_xml_node_name':
                n.set(k, str(kwargs[k]))

        self.current_geometry = n
        return self.current_geometry
    return fin


def defaultArg(name, value):
    def nif(func):
        def fin(self, *args, **kwargs):
            if not name in kwargs:
                kwargs[name] = value
            return func(self, *args, **kwargs)
        return fin
    return nif

def requireArg(name):
    def nif(func):
        def fin(self, *args, **kwargs):
            if not name in kwargs:
                raise BaseException("Argument "+str(name)+ " is required")
            return func(self, *args, **kwargs)
        return fin
    return nif



def addSimpleGeomElements(nameList):
    def nif(cls):
        def inf(name,fname,cls):
            @geometryElement
            def fin(self, **kwargs):
                kwargs['_xml_node_name'] = name
                return kwargs
            setattr(cls, "add"+fname, fin)                
        for n in nameList:
            if type(n) is tuple:
                inf(n[0], n[1],cls)
            else:
                inf(n,n,cls)
        return cls
    return nif

def addSimpleBCElements(nameList):
    def nif(cls):
        def inf(name,fname,cls):
            @BCElement
            def fin(self, **kwargs):
                kwargs['_xml_node_name'] = name
                return kwargs
            setattr(cls, "add"+fname, fin)                
        for n in nameList:
            if type(n) is tuple:
                inf(n[0], n[1],cls)
            else:
                inf(n,n,cls)
        return cls
    return nif

def addRootElements(nameList):
    def nif(cls):
        def inf(name,csl):
            @rootElement
            def fin(self, **kwargs):
                kwargs['_xml_node_name'] = name
                return kwargs
#            cls.__dict__["add"+name] = fin
            setattr(cls, "add"+name, fin)                
        for n in nameList:
            inf(n,cls)
        return cls
    return nif

    
def _set_by_kw(kw, name, default):
    if name in kw:
        return kw[name]
    else:
        return default
        
        
@addSimpleGeomElements([
    'Box',
    'Sphere',
    'HalfSphere',
    'OffgridSphere',
    'OffgridPipe',
    ('Outlet','OutletElement'),
    ('Inlet','InletElement'),
    ])
@addSimpleBCElements([
    'MRT',
    'TRT_SOI',    
    'Smoothing',
    'Cumulant',    
    'ESymmetry',
    'NSymmetry',
    'MovingWall',
    'SSymmetry',
    'None',
    'EPressure',
    'WPressure',
    'EVelocity',
    'WVelocity',    
    'SVelocity',
    'NVelocity',     
    ('Outlet','OutletDef'),
    ('Inlet','InletDef'),
    ('SolidBoundary1','SolidBoundary1Def'),
    ('SolidBoundary2','SolidBoundary2Def'),
    ('SolidBoundary3','SolidBoundary3Def'),    
    ])
@addRootElements(
[   
 'EvalIf',
 'CallPython',
 'VTK',
 'RunR',
 'Log',
 'Andersen'
 ]        
)
class CLBConfigWriter:

    def __init__(self, output="output/", sign=''):
        self.root = ET.Element('CLBConfig')
        if not sign == '':
            self.root.append(ET.Comment(sign))
        self.root.append(ET.Comment("Created using CLBConfigWriter"))

        self.geometry = ET.SubElement(self.root,'Geometry')
        self.model = ET.SubElement(self.root, 'Model')

        self.root.set("version", "2.0")
        self.root.set("output", output)
        self.geometry.set("predef", "none")
        self.geometry.set("model", "MRT")

        self.current_geometry = self.geometry


    def newModel(self):
        self.model = ET.SubElement(self.root, 'Model')

    def dump(self):
        self.indent(self.root)
        return ET.tostring(self.root,  pretty_print=True)

    def indent(self, elem, level=0):
        i = "\n" + level*"  "
        if len(elem):
            if not elem.text or not elem.text.strip():
                elem.text = i + "  "
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
            for elem in elem:
                self.indent(elem, level+1)
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                elem.tail = i


    def write(self, filename):
        tree = ET.ElementTree(element=self.root) 
        tree.write(filename,  pretty_print=True)

    def addModelParam(self, name, value,zone=''):
        n = ET.SubElement(self.model,'Param')
        try: 
            n.set("name", str(name))
            n.set("value", str('%.16f'%value))
            if zone != '':
                n.set("zone", str(zone))

        except TypeError as e:
            print(name, '==>', value)
            raise e
    
    def addModelParams(self, params_dict,zone=''):
        for n in params_dict:
            self.addModelParam(n, params_dict[n],zone)
            
    def addModelParams_old(self, params_dict):
        for n in params_dict:
            self.addModelParam_old(n, params_dict[n])      


    def addModelParam_old(self, name, value):
        n = ET.SubElement(self.model,'Params')
        try: 
            n.set(str(name), str('%.16f'%value))
        except TypeError as e:
            print(name, '==>', value)
            raise e
        
    def addRootParams(self, params_dict):
        for n in params_dict:
            self.addParamRoot(n, params_dict[n])      

    def addRootParam(self, name, value):
        n = ET.SubElement(self.root,'Params')
        n.set(str(name), str('%.16f'%value))
        
    def addParamRoot(self, name, value):
        self.addRootParam(name,value)
        
    def addGeomParam(self, name, value):
        self.geometry.set(str(name),  str('%.16f'%value))

    def setCG(self, cg):
        self.current_geometry = cg

    def addInit(self):
        ET.SubElement(self.root, 'Init')

#    def addSaveVTK(self,  **kwargs):
#        
#        vtk_what = _set_by_kw(kwargs, 'vtk_fields', 0)
#        
#        
#        n = ET.SubElement(self.root, 'VTK')
#        if vtk_what > 0:
#            n.set('what', str(vtk_what))       
    
    def addSaveCheckpoint(self, fname,comp=False, iterations=0):
        if not comp:
            n = ET.SubElement(self.root, 'SaveMemoryDump')
        else:
            n = ET.SubElement(self.root, 'SaveBinary')
            n.set('comp', comp)
        if iterations > 0:
            n.set('Iterations', str(iterations))
        n.set('filename', str(fname))

    def addLoadCheckpoint(self, fname,comp=False):
        if not comp:
            n = ET.SubElement(self.root, 'LoadMemoryDump')
        else:
            n = ET.SubElement(self.root, 'LoadBinary')
            n.set('comp', comp)

        n.set('filename', str(fname))


    def addSolve(self, **kwargs):
        
        iterations = _set_by_kw(kwargs, 'iterations', 1)
        vtk = _set_by_kw(kwargs, 'vtk', 0)
        vtk_what = _set_by_kw(kwargs, 'vtk_fields', '')
        log = _set_by_kw(kwargs, 'log', 0)
        failcheck = _set_by_kw(kwargs, 'failcheck', 0)

        failcheck_nx = _set_by_kw(kwargs, 'failcheck_nx', 1)
        failcheck_ny = _set_by_kw(kwargs, 'failcheck_ny', 1)
        failcheck_nz = _set_by_kw(kwargs, 'failcheck_nz', 1)

        failcheck_dx = _set_by_kw(kwargs, 'failcheck_dx', 1)
        failcheck_dy = _set_by_kw(kwargs, 'failcheck_dy', 1)
        failcheck_dz = _set_by_kw(kwargs, 'failcheck_dz', 1)
        failcheck_fields = _set_by_kw(kwargs, 'failcheck_fields', "")
        
        #self.model = self.root
        n = ET.SubElement(self.root, 'Solve')
        n.set('Iterations', str(iterations))
        if vtk > 0:
            n2 = ET.SubElement(n, 'VTK')
            n2.set('Iterations', str(vtk))
            if len(vtk_what) > 0:
                n2.set('what', str(vtk_what))       
                        
        if log > 0:
            n3 = ET.SubElement(n, 'Log')
            n3.set('Iterations', str(log))
        if failcheck > 0:
            n4 = ET.SubElement(n, 'Failcheck')
            n4.set('Iterations', str(failcheck))
            n4.set('nx', str(failcheck_nx))
            n4.set('ny', str(failcheck_ny))
            n4.set('nz', str(failcheck_nz))

            n4.set('dx', str(failcheck_dx))
            n4.set('dy', str(failcheck_dy))
            n4.set('dz', str(failcheck_dz))

            n4.set('what', str(failcheck_fields))


            ET.SubElement(n4, 'VTK')     
        return n
##############
# ELEMENT METHODE
#############
    @requireArg('name')
    @BCElement
    def addZoneBlock(self, **kwargs):
        kwargs['_xml_node_name'] = 'None'
        return kwargs

    @defaultArg('mask','ALL')
    @BCElement
    def addWall(self, **kwargs):
        kwargs['_xml_node_name'] = 'Wall'
        return kwargs

    @requireArg('file')
    @geometryElement
    def addText(self, **kwargs):
        kwargs['_xml_node_name'] = 'Text'
        return kwargs

    @requireArg('eval')
    @addCDATA('eval')
    @geometryElement
    def addPythonInline(self, **kwargs):
        kwargs['_xml_node_name'] = 'PythonInline'
        del kwargs['eval']
        return kwargs
##############
#  END ELEMENT FUNCTIONS, END CLASS
#############


if __name__ == "__main__":
    CLBc = CLBConfigWriter()
    CLBc.addGeomParam('ny', 256)
    CLBc.addGeomParam('nx', 160)
    
    CLBc.addModelParam('v',1)
    
    CLBc.newModel()
    CLBc.addModelParam('v',2)
    CLBc.addMRT()
    CLBc.addBox()
    
    CLBc.addZoneBlock(name='zwet')
    
    CLBc.addBox(dy=90, fy=-90)
    
    CLBc.addEPressure()
    CLBc.addInletElement()
    
    CLBc.addInletDef()
    CLBc.addBox(nx=1)

    CLBc.addSolve(iterations=1, vtk=1)
    CLBc.addSolve(iterations=100, vtk=50, failcheck=1)
    
    print(CLBc.dump())
   # CLBc.write('/home/michal/tach-17/mnt/fhgfs/users/mdzikowski/yang-laplace-sphere-matrix/test.xml')
   