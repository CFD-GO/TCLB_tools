# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 15:55:33 2015

@author: mdzikowski
"""

import xml.sax
import re
import numpy as np

def toFloat(val):
    try:
        return float(val)
    except ValueError:
        return np.nan
        
class CLBXMLHandler(xml.sax.ContentHandler):
    
    def __init__(self, config_ref, mp, time):
        self.config = config_ref
        self.mp = mp
        self.iterations = 0
        self.time =  time
        self.rewind = False
        if self.mp:
            self.startElement = self.startElement_mp
        else:
            self.startElement = self.startElement_sp            
        
    def startElement_sp(self, name, attrs):
        
        if name == "Solve":
            iters = attrs.items()
            if self.iterations >= self.time:
                self.rewind = True
            for (k,v) in iters:
                if k == 'Iterations':
                    try:
                        self.iterations = self.iterations + int(v)
                    except ValueError:
                        self.iterations = self.iterations + 0
            #print self.iterations
                
        if self.rewind:
                return
        
        if name == "Params":
            a = dict()
            for (k,v) in attrs.items():
                if k == 'gauge':
                    g = re.findall('[-,\.,e,0-9]+', v)
                    if len(g) == 1:
                        a['gauge'] = toFloat(g[0])
                    else:
                        print("No gauge value: ", a, v)
                        a['float'] = np.nan                        
                else:
                    a['name'] = k
                    a['value'] = v
                    g = re.findall('[-,\.,e,0-9]+', v)
                    if len(g) >= 1:
                        a['float'] = toFloat(g[0])                    
                    else:
                        print("No value: ", k, v)
                        a['float'] = np.nan
            
            ###catch overwrite
            if (a['name'] in self.config) and ('gauge' in self.config[a['name']]) :
                tmp = self.config[a['name']]
                for k in a.keys():
                    if not tmp.has_key(k):
                        tmp[k] = a[k]                        
                a = tmp
                print("Safe-Overwriting ", a['name'])
                
            self.config[a['name']] = a
            self.config[a['name']]['time'] = self.iterations
            
        if name == "Param":
            a = dict()

            v = attrs['value']
            a['name'] = attrs['name']

            if 'zone' in attrs.keys():
                a['name'] = a['name'] + '-'  + attrs['zone']
                
            g = re.findall('[-,\.,e,0-9]+', v)
            if len(g) >= 1:
                a['float'] = toFloat(g[0])                    
            else:
                print("No value: ", k, v)
                a['float'] = np.nan
            
            ###catch overwrite
            if (a['name'] in self.config) and ('gauge' in self.config[a['name']]) :
                tmp = self.config[a['name']]
                for k in a.keys():
                    if not tmp.has_key(k):
                        tmp[k] = a[k]                        
                a = tmp
                print("Safe-Overwriting ", a['name'])
                
            self.config[a['name']] = a
            self.config[a['name']]['time'] = self.iterations
            
        if name == "Geometry":
            
            for (k,v) in attrs.items():
                
                if k in ('nx', 'ny', 'nz'):
                    a = dict()
                    a['name'] = k
                    a['value'] = v
                    g = re.findall('[-,\.,e,0-9]+', v)
                    if len(g) == 1:
                        a['float'] = toFloat(g[0])         
                    else:
                        a['float'] = np.nan                            
                    self.config[a['name']] = a

            
                    
    def startElement_mp(self, name, attrs):
        if name == "Params":
            
            for (k,v) in attrs.items():
                a = dict()
                if k == 'gauge':
                    g = re.findall('[-,\.,e,0-9]+', v)
                    if len(g) == 1:
                        a['gauge'] = toFloat(g[0])
                    else:
                        print("No gauge value: ", k, v)
                        a['float'] = np.nan                              
                else:
                    a['name'] = k
                    a['value'] = v
                    g = re.findall('[-,\.,e,0-9]+', v)
                    if len(g) >= 1:
                        a['float'] = toFloat(g[0])                    
                    else:
                        print("No value: ", k, v)
                        a['float'] = np.nan
                self.config[a['name']] = a
        if name == "Geometry":
            
            for (k,v) in attrs.items():
                
                if k in ('nx', 'ny'):
                    a = dict()
                    a['name'] = k
                    a['value'] = v
                    g = re.findall('[-,\.,e,0-9]+', v)
                    if len(g) == 1:
                        a['float'] = toFloat(g[0])         
                    else:
                        a['float'] = np.nan                            
                    self.config[a['name']] = a        

            
def parseConfig(fconfig, **kwargs):
    CLBc = dict()
    parser = xml.sax.make_parser()
    if 'multiparams' in kwargs and kwargs['multiparams']:
        mp = True
    else:
        mp = False
        
    if 'time' in kwargs:
        time = kwargs['time']
    else:
        time = 0
    parser.setContentHandler(CLBXMLHandler(CLBc, mp, time))
    
    
    try:
        parser.parse(open(fconfig,"r"))
    except xml.sax._exceptions.SAXParseException:
        # this may be due <Run> element at the end
        parser.reset()
        f = file(fconfig,"r")
        temp = list()
        for l in f:
            temp.append(str(l))  
        t = '\n'.join(temp[:-3])
        parser.feed(t)
        
        
    CLBcf = dict()
    CLBcg = dict()
    for c in  CLBc:
        CLBcf[c] = CLBc[c]['float']
        if 'gauge' in CLBc[c]:
            CLBcg[c] = CLBc[c]['gauge']
       #     print c, "  gauge = ", CLBc[c]['gauge']        
       # print c, " = ", CLBc[c]['float']
        
    return CLBc, CLBcf, CLBcg

