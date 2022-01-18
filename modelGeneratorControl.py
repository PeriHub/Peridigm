# Copyright (C) 2021 Deutsches Zentrum fuer Luft- und Raumfahrt(DLR, German Aerospace Center) <www.dlr.de>


from api.app.GIICmodel.GIICmodel import GIICmodel
from api.app.DCBmodel.DCBmodel import DCBmodel
from api.app.Verification.verificationModels import VerificationModels
#from XFEM_Bechnmark.XFEMdcb import XFEMDCB
from api.app.Dogbone.Dogbone import Dogbone
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
        h = 0.019
        nn = 21
        
        
        
        nn = 2*int(nn/2)+1
        nn = 2*int(nn/2)+1
        dx=[h/nn,h/nn,h/nn]
        
        print(dx, 4.01*dx[0])
        
        gc = GIICmodel(xend = L, yend = h, zend = B, dx=dx, TwoD = True)
        model = gc.createModel()
        #xm = XFEMDCB(xend = L, yend = 2*h, dx=[0.08,0.08])
        # dx=[0.00025,0.00025,0.00025]
        # db = DCBmodel(dx = dx, TwoD = True)
        # model = db.createModel()
        # db = Dogbone(dx = dx, TwoD = False, h1=h)
        # model = db.createModel()
        #veri = VerificationModels()
        #veri.createVerificationModels()
    def endRunOnError(self):
        pass
        
        
    def endRun(self, returnDir = None, feFilename = None, runDir = None):
       pass
    
            
        
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        

