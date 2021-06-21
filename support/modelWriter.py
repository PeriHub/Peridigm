import numpy as np
import os
from support.xmlCreator  import XMLcreator
from support.yamlCreator  import YAMLcreator

class ModelWriter(object):
    def __init__(self, filename = 'mesh'):
        
        self.filename = filename
        self.nsName   = 'ns_' + filename
        self.path = 'Output/'+filename
        if not os.path.exists('Output'):
            os.mkdir('Output')   
    def writeNodeSets(self, model, nslist):
        for idx, k in enumerate(nslist):
            points = np.where(model['k'] == k)
            string = ''
            for pt in points[0]:
                string += str(int(pt)+1) + '\n'
            self.fileWriter(self.nsName + '_' + str(idx+1) + '.txt', string)

    def fileWriter(self, filename, string):
        if not os.path.exists(self.path):
            os.mkdir(self.path)
        fid = open(self.path+'/'+filename,'w')
        fid.write(string)
        fid.close()
    def writeMesh(self, model):    
        string = '# x y z block_id volume\n'
        for idx in range(0, len(model['x'])):
            string += str(model['x'][idx]) + " " + str(model['y'][idx])+ " " + str(model['z'][idx]) + " " + str(model['k'][idx]) + " " + str(model['vol'][idx]) + "\n"
        self.fileWriter(self.filename + '.txt', string)
    def writeMeshWithAngles(self, model):    
        string = '# x y z block_id volume angle_x angle_y angle_z\n'
        for idx in range(0, len(model['x'])):
            string += str(model['x'][idx]) + " " + str(model['y'][idx])+ " " + str(model['z'][idx]) + " " + str(model['k'][idx]) + " " + str(model['vol'][idx]) + " " + str(model['angle_x'][idx]) +" " + str(model['angle_y'][idx]) +" " + str(model['angle_z'][idx]) +"\n"
        self.fileWriter(self.filename + '.txt', string)       
    def createFile(self, filetype, solvertype, bcDict,damageDict, materialDict, blockDef, bondfilters,TwoD):
        
        if filetype == 'yaml':
            yl = YAMLcreator(filename = self.filename, nsName = self.nsName, solvertype = solvertype, bc = bcDict, damageDict = damageDict, materialDict = materialDict, blockDef = blockDef, bondfilters = bondfilters, TwoD = TwoD)
            string = yl.createYAML()
            self.fileWriter(self.filename + '.yaml', string)
            
        elif filetype == 'xml':
            xl = XMLcreator(filename = self.filename, nsName = self.nsName, solvertype = solvertype, bc = bcDict, damageDict = damageDict, materialDict = materialDict, blockDef = blockDef, bondfilters = bondfilters, TwoD = TwoD)
            string = xl.createXML()
        else:
            print('Not a supported filetye: ', filetype)    
        self.fileWriter(self.filename + '.' + filetype, string)
            
        
