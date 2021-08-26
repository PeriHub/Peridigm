# Copyright (C) 2021 Deutsches Zentrum fuer Luft- und Raumfahrt(DLR, German Aerospace Center) <www.dlr.de>


from ntpath import join
from numpy import string_
from numpy.lib.shape_base import split
from GIICmodel.GIICmodel import GIICmodel
from DCBmodel.DCBmodel import DCBmodel
from support.sbatchCreator  import SbatchCreator
#from XFEM_Bechnmark.XFEMdcb import XFEMDCB
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt

import pandas as pd

from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi.responses import StreamingResponse
from fastapi.middleware.cors import CORSMiddleware
import uvicorn

from enum import Enum
from pydantic import BaseModel
from typing import Dict, List

import shutil
import paramiko
import os
import csv


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

    def __init__(self,**kwargs):
        """doc"""
        self.returnDir = None

    def run(self,**kwargs):
        """doc"""
        
        L = 152
        L = 50
        W = 10
        h = 4.95
        nn = 11

        nn = 2*int(nn/2)+1
        dx=[h/nn,h/nn,h/nn]
        
        print(dx, 4.01*dx[0])
        
        # dx=[0.001,0.001,0.001]
        # db = DCBmodel(dx = dx, solvertype = 'Verlet', TwoD = True, filetype = 'xml')
        # model = db.createModel()
        gc = GIICmodel(xend = L, yend = h, zend = W, dx=dx, solvertype = 'Verlet', TwoD = False, filetype = 'yaml', rot=False)
        model = gc.createModel()
        #xm = XFEMDCB(xend = L, yend = 2*h, dx=[0.08,0.08])

    def endRunOnError(self):
        pass
        
    def endRun(self, returnDir = None, feFilename = None, runDir = None):
       pass


    @app.post("/generateModel")
    def generateModel(ModelName: ModelName, Length: float, Width: float, Height: float, Discretization: float, TwoDimensional: bool, RotatedAngles: bool, Angle0: float, Angle1: float, Param: dict):#Material: dict, Output: dict):
        # L = 152
        # L = 50
        # W = 10
        # h = 4.95
        # nn = 12

        L = Length
        W = Width
        h = Height
        nn = Discretization
        dx=[h/nn,h/nn,h/nn]

        if ModelName==ModelName.GIICmodel:
            gc = GIICmodel(xend = L, yend = h, zend = W, dx=dx, solvertype = Param['Param']['Solver']['solvertype'], finalTime = Param['Param']['Solver']['finalTime'], TwoD = TwoDimensional, filetype = Param['Param']['Solver']['filetype'], rot=RotatedAngles, angle=[Angle0,Angle1], material=Param['Param']['Material'], damage=Param['Param']['Damage'], block=Param['Param']['Block'], bc=Param['Param']['BoundaryConditions'], output=Param['Param']['Output'], solver=Param['Param']['Solver'])
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

        return {ModelName: ModelName + "-InputFile has been saved"}

    @app.get("/getModel")
    def getModel(ModelName: ModelName):

        shutil.make_archive(ModelName, "zip", './Output/' + ModelName)

        response = FileResponse(ModelName + ".zip", media_type="application/x-zip-compressed")
        response.headers["Content-Disposition"] = "attachment; filename=" + ModelName + ".zip"
        # return StreamingResponse(iterfile(), media_type="application/x-zip-compressed")
        return response

    @app.get("/getResults")
    def getResults(ModelName: ModelName, Cluster: str):

        if Cluster=='FA-Cluster':
            username='hess_ja'
            server='129.247.54.37'
            keypath = 'id_rsa_cluster'
            remotepath = './Peridigm/apiModels/' + ModelName
            ssh = paramiko.SSHClient() 
            ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            ssh.connect(server, username=username, allow_agent=False, key_filename=keypath)
            command = 'cd ' + remotepath + ' \n qperidigm -d -c ' + str(Param['Param']['Job']['Tasks']) + ' -O tgz -J ' + ModelName +' -E /home/hess_ja/PeridigmInstall/build/bin/Peridigm '+ ModelName + '.' + FileType
            ssh.exec_command(command)
            ssh.close()

            shutil.make_archive(ModelName, "zip", './Output/' + ModelName)

            response = FileResponse(ModelName + ".zip", media_type="application/x-zip-compressed")
            response.headers["Content-Disposition"] = "attachment; filename=" + ModelName + ".zip"
            # return StreamingResponse(iterfile(), media_type="application/x-zip-compressed")
            return response

        elif Cluster=='Cara':
            username='hess_ja'
            server='cara.dlr.de'
            keypath = 'id_rsa_cara'
            remotepath = './PeridigmJobs/apiModels/' + ModelName
            ssh = paramiko.SSHClient() 
            ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            ssh.connect(server, username=username, allow_agent=False, key_filename=keypath)
            sftp = ssh.open_sftp()
            sftp.get (remotepath + '/' + ModelName + '.sbatch', "a", -1)
            file.write(sbatchString)
            file.flush()
            sftp.close()

            command = 'cd ' + remotepath + ' \n sbatch '+ ModelName + '.sbatch'
            ssh.exec_command(command)
            ssh.close()

            shutil.make_archive(ModelName, "zip", './Output/' + ModelName)

            response = FileResponse(ModelName + ".zip", media_type="application/x-zip-compressed")
            response.headers["Content-Disposition"] = "attachment; filename=" + ModelName + ".zip"
            # return StreamingResponse(iterfile(), media_type="application/x-zip-compressed")
            return response

        else:
            return {Cluster: Param['Param']['Job']['cluster'] + "unknown"}
            

    @app.get("/getPointData")
    def getPointData(ModelName: ModelName):

        pointString=''
        blockIdString=''
        firstRow=True  
        with open('./Output/' + "GIICmodel" + '/'  + "GIICmodel" + '.txt', 'r') as f:
                reader = csv.reader(f)
                for row in reader:
                    if firstRow==False:
                        str1 = ''.join(row)
                        parts = str1.split(" ")
                        pointString+=parts[0]+','+parts[1]+','+parts[2]+','
                        blockIdString+=str(int(parts[3])/10)+','
                    firstRow=False
                    
        response=[pointString.rstrip(pointString[-1]),blockIdString.rstrip(blockIdString[-1])]
        return response

    @app.get("/copyModelToCluster")
    def copyModelToCluster(ModelName: ModelName, Cluster: str):

        if Cluster=='FA-Cluster':
            username='hess_ja'
            server='129.247.54.37'
            keypath = 'id_rsa_cluster'
            remotepath = './Peridigm/apiModels/' + ModelName
            ssh = paramiko.SSHClient() 
            ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            ssh.connect(server, username=username, allow_agent=False, key_filename=keypath)
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
        
        elif Cluster=='Cara':
            username='hess_ja'
            server='cara.dlr.de'
            keypath = 'id_rsa_cara'
            remotepath = './PeridigmJobs/apiModels/' + ModelName
            ssh = paramiko.SSHClient() 
            ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            ssh.connect(server, username=username, allow_agent=False, key_filename=keypath)
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

        else:
            return {Cluster: Cluster + "unknown"}


    @app.post("/runModel")
    def runModel(ModelName: ModelName, FileType: FileType, Param: dict):

        if Param['Param']['Job']['cluster']=='FA-Cluster':
            username='hess_ja'
            server='129.247.54.37'
            keypath = 'id_rsa_cluster'
            remotepath = './Peridigm/apiModels/' + ModelName
            ssh = paramiko.SSHClient() 
            ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            ssh.connect(server, username=username, allow_agent=False, key_filename=keypath)
            command = 'cd ' + remotepath + ' \n qperidigm -d -c ' + str(Param['Param']['Job']['Tasks']) + ' -O tgz -J ' + ModelName +' -E /home/hess_ja/PeridigmInstall/build/bin/Peridigm '+ ModelName + '.' + FileType
            ssh.exec_command(command)
            ssh.close()

            return {ModelName: ModelName + "-Model has been submitted"}

        elif Param['Param']['Job']['cluster']=='Cara':
            sb = SbatchCreator(filename=ModelName, output=Param['Param']['Output'], job=Param['Param']['Job'])
            sbatchString = sb.createSbatch()
            username='hess_ja'
            server='cara.dlr.de'
            keypath = 'id_rsa_cara'
            remotepath = './PeridigmJobs/apiModels/' + ModelName
            ssh = paramiko.SSHClient() 
            ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            ssh.connect(server, username=username, allow_agent=False, key_filename=keypath)
            sftp = ssh.open_sftp()
            file=sftp.file(remotepath + '/' + ModelName + '.sbatch', "a", -1)
            file.write(sbatchString)
            file.flush()
            sftp.close()

            command = 'cd ' + remotepath + ' \n sbatch '+ ModelName + '.sbatch'
            ssh.exec_command(command)
            ssh.close()

            return {ModelName: ModelName + "-Model has been submitted"}

        else:
            return {Cluster: Param['Param']['Job']['cluster'] + "unknown"}

    # uvicorn.run(app, host="0.0.0.0", port=8000)    
    # pointString=''
    # firstRow=True  
    # with open('./Output/' + "GIICmodel" + '/'  + "GIICmodel" + '.txt', 'r') as f:
    #         reader = csv.reader(f)
    #         for row in reader:
    #             if firstRow==False:
    #                 str1 = ''.join(row)
    #                 parts = str1.split(" ")
    #                 pointString+=parts[0]+','+parts[1]+','+parts[2]+','
    #             firstRow=False
                
    # pointString=pointString.rstrip(pointString[-1])

    # from pyevtk.hl import pointsToVTK
    # import meshio
    # import numpy as np 
    # npoints = 100 
    # x = np.random.rand(npoints) 
    # y = np.random.rand(npoints) 
    # z = np.random.rand(npoints) 
    # pressure = np.random.rand(npoints) 
    # temp = np.random.rand(npoints) 
    # test = pointsToVTK("./points", x, y, z, data = {"temp" : temp, "pressure" : pressure})
    # mesh = meshio.read("./points.vtu")
    # points = [
    # [0.0, 0.0],
    # [1.0, 0.0],
    # [0.0, 1.0],
    # [1.0, 1.0],
    # [2.0, 0.0],
    # [2.0, 1.0],
    # ]
    # cells = [
    #     ("triangle", [[0, 1, 2], [1, 3, 2]]),
    #     ("quad", [[1, 4, 5, 3]]),
    # ]

    # mesh = meshio.Mesh(
    #     points,
    #     cells,
    #     # Optionally provide extra data on points, cells, etc.
    #     point_data={"T": [0.3, -1.2, 0.5, 0.7, 0.0, -3.0]},
    #     # Each item in cell data must match the cells array
    #     cell_data={"a": [[0.1, 0.2], [0.4]]},
    # )
    # mesh.write("./foo.stl")
    # print(mesh)
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        

