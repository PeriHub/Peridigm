# Copyright (C) 2021 Deutsches Zentrum fuer Luft- und Raumfahrt(DLR, German Aerospace Center) <www.dlr.de>


from GIICmodel.GIICmodel import GIICmodel
#from XFEM_Bechnmark.XFEMdcb import XFEMDCB
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt

import pandas as pd

class ModelControl(object):


    def __init__(self,**kwargs):
        """doc"""
        self.returnDir = None

    def run(self,**kwargs):
        """doc"""
        
        L = 152
        L = 52
        B = 10
        h = 4.95
        nn = 19
        
        
        
        nn = 2*int(nn/2)+1
        dx=[h/nn,h/nn,h/nn]
        
        print(dx, 4.01*dx[0])
        
        gc = GIICmodel(xend = L, yend = h, zend = B, dx=dx, solvertype = 'Verlet', TwoD = True, filetype = 'xml')
        model = gc.createModel(rot=True)
        #xm = XFEMDCB(xend = L, yend = 2*h, dx=[0.08,0.08])
    def endRunOnError(self):
        pass
        
        
    def endRun(self, returnDir = None, feFilename = None, runDir = None):
       pass
    
            
        
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        

