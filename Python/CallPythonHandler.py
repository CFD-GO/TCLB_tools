# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 09:22:24 2016

@author: mdzik
"""

import CLB.CallPythonHelper as CPH
import IPython
from IPython.config.loader import Config
from IPython.frontend.terminal.embed import InteractiveShellEmbed
def console(*args):
    
    cph = CPH.CallPythonHelper(*args)
        
    
    cfg = Config()
    cfg.InteractiveShellEmbed.prompt_in1="CLB [\\#]> "
    cfg.InteractiveShellEmbed.prompt_out="CLB [\\#]: "
    #cfg.InteractiveShellEmbed.profile=ipythonprofile
    # directly open the shell
    baner = """
    
                     (          
                (    )\ )   (   
        (       )\  (()/( ( )\  
 `  )   )\ )  (((_)  /(_)))((_) 
 /(/(  (()/(  )\___ (_)) ((_)_  
((_)_\  )(_))((/ __|| |   | _ ) 
| '_ \)| || | | (__ | |__ | _ \ 
| .__/  \_, |  \___||____||___/ 
|_|     |__/                        
    

    example: 
    #cph is already definied for you :)
    
    Fx = cph.getScalar(0)
    Fy = cph.getVector(1)
    X,Y = cph.getXY(0)
    
    """
    IPython.embed(config=cfg, banner2=baner)
    # or get shell object and open it later
    
    #shell = InteractiveShellEmbed(config=cfg, user_ns=namespace, banner2=banner)
    #shell.user_ns = {}
    #shell()


    return 0


def dump(*args):
    cph = CPH.CallPythonHelper(*args)
    cph.dump('dump_'+str(cph.time)+'.npz')    
    
    
    
if "__main__" == __name__:
    console('a')