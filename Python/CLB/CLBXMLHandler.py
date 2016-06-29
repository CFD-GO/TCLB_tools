# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 15:55:33 2015

@author: mdzikowski
"""

import xml.sax
import re
import numpy as np
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
                    self.iterations = self.iterations + int(v)
            
                
        if self.rewind:
                return
        
        if name == "Params":
            a = dict()
            for (k,v) in attrs.items():
                if k == 'gauge':
                    g = re.findall('[-,\.,e,0-9]+', v)
                    if len(g) == 1:
                        a['gauge'] = float(g[0])
                else:
                    a['name'] = k
                    a['value'] = v
                    g = re.findall('[-,\.,e,0-9]+', v)
                    if len(g) >= 1:
                        a['float'] = float(g[0])                    
                    else:
                        print "No value: ", k, v
                        a['float'] = np.nan
            self.config[a['name']] = a
            self.config[a['name']]['time'] = self.iterations
        if name == "Geometry":
            
            for (k,v) in attrs.items():
                
                if k in ('nx', 'ny'):
                    a = dict()
                    a['name'] = k
                    a['value'] = v
                    g = re.findall('[-,\.,e,0-9]+', v)
                    if len(g) == 1:
                        a['float'] = float(g[0])         
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
                        a['gauge'] = float(g[0])
                else:
                    a['name'] = k
                    a['value'] = v
                    g = re.findall('[-,\.,e,0-9]+', v)
                    if len(g) >= 1:
                        a['float'] = float(g[0])                    
                    else:
                        print "No value: ", k, v
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
                        a['float'] = float(g[0])         
                    else:
                        a['float'] = np.nan                            
                    self.config[a['name']] = a        
            
def parseConfig(fconfig, **kwargs):
    CLBc = dict()
    parser = xml.sax.make_parser()
    if kwargs.has_key('multiparams') and kwargs['multiparams']:
        mp = True
    else:
        mp = False
        
    if kwargs.has_key('time'):
        time = kwargs['time']
    else:
        time = 0
    parser.setContentHandler(CLBXMLHandler(CLBc, mp, time))
    parser.parse(open(fconfig,"r"))

    
    CLBcf = dict()
    CLBcg = dict()
    for c in  CLBc:
        CLBcf[c] = CLBc[c]['float']
        if 'gauge' in CLBc[c]:
            CLBcg[c] = CLBc[c]['gauge']
       #     print c, "  gauge = ", CLBc[c]['gauge']        
       # print c, " = ", CLBc[c]['float']
        
    return CLBc, CLBcf, CLBcg