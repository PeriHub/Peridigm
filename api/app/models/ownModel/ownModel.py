from support.modelWriter import ModelWriter

class OwnModel(object):
    def __init__(self, dx=[0.0005,0.0005,0.0005], DiscType = 'txt', TwoD = False, horizon = 0.1, filename = 'ownModel', material = '', damage = '', block = '', bc = '', bf = '', compute = '', output = '', solver = '', username = ''):
        
        self.filename = filename
        self.scal = 1
        self.DiscType = DiscType
        self.TwoD = TwoD
        self.horizon = horizon
        self.dx   = dx
        self.blockDef = block
        self.username = username
        self.damageDict = damage
        self.computeDict = compute
        self.outputDict = output
        self.materialDict = material
        self.bondfilters = bf
        self.bcDict = bc
        self.solverDict = solver

    def createModel(self):

        writer = ModelWriter(modelClass = self)
        self.writeFILE(writer = writer)
        
        return 'Model created'

    def writeFILE(self, writer):

        for idx in range(0,len(self.blockDef)):
            self.blockDef[idx]['horizon']= self.horizon

        writer.createFile(self.blockDef)