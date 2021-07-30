# Copyright (C) 2021 Deutsches Zentrum fuer Luft- und Raumfahrt(DLR, German Aerospace Center) <www.dlr.de>


from ntpath import join
from numpy import string_
from numpy.lib.shape_base import split
from GIICmodel.GIICmodel import GIICmodel
from DCBmodel.DCBmodel import DCBmodel
#from XFEM_Bechnmark.XFEMdcb import XFEMDCB
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt

import pandas as pd

from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi.responses import StreamingResponse
from fastapi.middleware.cors import CORSMiddleware

from enum import Enum

import shutil
import paramiko
import os
class ModelName(str, Enum):
    GIICmodel = "GIICmodel"
    DCBmodel = "DCBmodel"
class SolverType(str, Enum):
    Verlet = "Verlet"
    NOXQuasiStatic = "NOXQuasiStatic"
class FileType(str, Enum):
    yaml = "yaml"
    xml = "xml"

app = FastAPI()

origins = [
    "http://localhost",
    "http://localhost:8080",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
class ModelControl(object):

    @app.post("/generateModel")
    def generateModel(ModelName: ModelName, Length: float, Width: float, Height: float, Discretization: float, TwoDimensional: bool, RotatedAngles: bool, Solvertype: SolverType, Filetype: FileType):
        # L = 152
        # L = 50
        # W = 10
        # h = 4.95
        # nn = 11

        L = Length
        W = Width
        h = Height
        nn = Discretization
        dx=[h/nn,h/nn,h/nn]

        if ModelName==ModelName.GIICmodel:
            gc = GIICmodel(xend = L, yend = h, zend = W, dx=dx, solvertype = Solvertype, TwoD = TwoDimensional, filetype = Filetype, rot=RotatedAngles)
            model = gc.createModel()

        if ModelName==ModelName.DCBmodel:
            dcb = DCBmodel(xend = L, yend = h, zend = W, dx=dx, solvertype = Solvertype, TwoD = TwoDimensional, filetype = Filetype, rot=RotatedAngles)
            model = dcb.createModel()

        return {ModelName: ModelName + "-Model has been created", 'dx': dx}
        
        
    @app.get("/viewInputFile")
    def viewInputFile(ModelName: ModelName, FileType: FileType):

        return FileResponse('./Output/' + ModelName + '/'  + ModelName + '.' + FileType)

    @app.post("/writeInputFile")
    def writeInputFile(ModelName: ModelName, InputString: str, FileType: FileType):

        fid = open('./Output/' + ModelName + '/'  + ModelName + '.' + FileType ,'w')
        fid.write(InputString)
        fid.close()

        return {ModelName: ModelName + "-InputFile has been written"}

    @app.get("/getModel")
    def getModel(ModelName: ModelName):

        shutil.make_archive(ModelName, "zip", './Output/' + ModelName)

        response = FileResponse(ModelName + ".zip", media_type="application/x-zip-compressed")
        response.headers["Content-Disposition"] = "attachment; filename=" + ModelName + ".zip"
        # return StreamingResponse(iterfile(), media_type="application/x-zip-compressed")
        return response

    @app.get("/copyModelToCluster")
    def copyModelToCluster(ModelName: ModelName):

        username='hess_ja'
        server='129.247.54.37'
        remotepath = './Peridigm/apiModels/' + ModelName
        ssh = paramiko.SSHClient() 
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(server, username=username)
        sftp = ssh.open_sftp()
        try:
            sftp.chdir(remotepath)  # Test if remote_path exists
        except IOError:
            sftp.mkdir(remotepath)  # Create remote_path
            sftp.chdir(remotepath)
        for root, dirs, files in os.walk('./Output/' + ModelName):
            for name in files:
                sftp.put(os.path.join(root,name), name)

        sftp.close()
        ssh.close()

        return {ModelName: ModelName + "-Model has been copied to Cluster"}

    @app.get("/runModel")
    def runModel(ModelName: ModelName, FileType: FileType, Tasks: int):

        username='hess_ja'
        server='129.247.54.37'
        remotepath = './Peridigm/apiModels/' + ModelName
        ssh = paramiko.SSHClient() 
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(server, username=username)
        command = 'cd ' + remotepath + ' \n qperidigm -d -c ' + str(Tasks) + ' -O tgz -J ' + ModelName +' -E /home/hess_ja/PeridigmInstall/build/bin/Peridigm '+ ModelName + '.' + FileType
        ssh.exec_command(command)
        ssh.close()

        return {ModelName: ModelName + "-Model has been submitted"}
        
    def __init__(self,**kwargs):
        """doc"""
        self.returnDir = None

    def run(self,**kwargs):
        """doc"""
        
        L = 152
        L = 50
        W = 10
        h = 4.95
        nn = 31
        dx=[h/nn,h/nn,h/nn]
        
        print(dx, 1.92/dx[0])
        
        # dcb = DCBmodel()
        # model = dcb.createModel()
        gc = GIICmodel(xend = L, yend = h, zend = W, dx=dx, solvertype = 'Verlet', TwoD = False, filetype = 'xml', rot=False)
        model = gc.createModel()
        #xm = XFEMDCB(xend = L, yend = 2*h, dx=[0.08,0.08])

    def endRunOnError(self):
        pass
        
        
    def endRun(self, returnDir = None, feFilename = None, runDir = None):
       pass
    
            
        
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        

