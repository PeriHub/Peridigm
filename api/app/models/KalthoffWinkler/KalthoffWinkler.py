import numpy as np
from support.baseModels import (
    Adapt,
    Block,
    BondFilters,
    BoundaryConditions,
    Compute,
    Damage,
    Material,
    Output,
    Newton,
    Solver,
    Verlet,
)
from support.modelWriter import ModelWriter
from support.geometry import Geometry


class KalthoffWinkler(object):
    def __init__(
        self,
        xend=100,
        yend=100,
        zend=0.003,
        dx=[0.001, 0.001, 0.001],
        filename="Kalthoff-Winkler",
        TwoD=False,
        rot="False",
        angle=[0, 0],
        material="",
        damage="",
        block="",
        bc="",
        bf="",
        compute="",
        output="",
        solver="",
        username="",
        maxNodes=10000000,
        ignoreMesh=False,
    ):
        """
        definition der blocks
        k =
        1 basisplatte
        2 X Zero Node Set
        3 X Zero Node Set
        4 Load Node Set
        """

        self.filename = filename
        self.scal = 4.01
        self.DiscType = "txt"
        self.TwoD = TwoD
        self.nsList = [3, 4]
        self.dx = dx
        self.xbegin = 0.0
        self.ybegin = -yend / 2
        self.xend = xend + dx[0]
        self.yend = yend / 2 + dx[1]
        # self.xend = xend
        # self.yend = yend/2
        # self.zend = zend
        self.rot = rot
        self.blockDef = block
        self.username = username
        self.maxNodes = maxNodes
        self.ignoreMesh = ignoreMesh
        if self.TwoD:
            self.zbegin = 0
            self.zend = 0
            self.dx[2] = 1
        else:
            self.zbegin = -zend
            self.zend = zend + dx[2]

        numberOfBlocks = 4

        """ Definition of model
        """
        matNameList = ["PMMA"]
        self.materialDict = []
        self.angle = [0, 0]
        if damage == "":
            damageDict = Damage(
                id=1,
                Name="PMMADamage",
                damageModel="Critical Stretch",
                criticalStretch=10,
                criticalEnergy=0.0022170,
                interblockdamageEnergy=0.01,
                planeStress=True,
                onlyTension=False,
                detachedNodesCheck=True,
                thickness=10,
                hourglassCoefficient=1.0,
                stabilizatonType="Global Stiffness",
            )
            self.damageDict = [damageDict]
        else:
            self.damageDict = damage

        if compute == "":
            computeDict1 = Compute(
                id=1,
                Name="External_Displacement",
                variable="Displacement",
                calculationType="Minimum",
                blockName="block_3",
            )
            computeDict2 = Compute(
                id=2,
                Name="External_Force",
                variable="Force",
                calculationType="Sum",
                blockName="block_3",
            )
            self.computeDict = [computeDict1, computeDict2]
        else:
            self.computeDict = compute

        if output == "":
            outputDict1 = Output(
                id=1,
                Name="Output1",
                Displacement=True,
                Force=True,
                Damage=True,
                Velocity=True,
                Partial_Stress=False,
                External_Force=False,
                External_Displacement=False,
                Number_Of_Neighbors=True,
                Frequency=10,
                InitStep=0,
            )
            self.outputDict = [outputDict1]
        else:
            self.outputDict = output

        if material == "":
            i = 0
            for material in matNameList:
                matDict = Material(
                    id=i + 1,
                    Name=material,
                    MatType="Elastic Bond Based",
                    density=8.0e-6,
                    bulkModulus=1.2666666666666667e5,
                    shearModulus=None,
                    youngsModulus=None,
                    poissonsRatio=None,
                    tensionSeparation=False,
                    nonLinear=True,
                    planeStress=True,
                    materialSymmetry="Isotropic",
                    stabilizatonType="Global Stiffness",
                    thickness=10.0,
                    hourglassCoefficient=1.0,
                    actualHorizon=None,
                    yieldStress=None,
                    Parameter=[],
                    Properties=[],
                )
                i += 1
                self.materialDict.append(matDict)
        else:
            self.angle = angle
            self.materialDict = material

        if bf == "":
            bf1 = BondFilters(
                id=1,
                Name="bf_1",
                type="Rectangular_Plane",
                normalX=0.0,
                normalY=1.0,
                normalZ=0.0,
                lowerLeftCornerX=-0.5,
                lowerLeftCornerY=25.0,
                lowerLeftCornerZ=-0.5,
                bottomUnitVectorX=1.0,
                bottomUnitVectorY=0.0,
                bottomUnitVectorZ=0.0,
                bottomLength=50.5,
                sideLength=1.0,
                centerX=0.0,
                centerY=1.0,
                centerZ=0.0,
                radius=1.0,
                show=True,
            )
            bf2 = BondFilters(
                id=1,
                Name="bf_2",
                type="Rectangular_Plane",
                normalX=0.0,
                normalY=1.0,
                normalZ=0.0,
                lowerLeftCornerX=-0.5,
                lowerLeftCornerY=-25.0,
                lowerLeftCornerZ=-0.5,
                bottomUnitVectorX=1.0,
                bottomUnitVectorY=0.0,
                bottomUnitVectorZ=0.0,
                bottomLength=50.5,
                sideLength=1.0,
                centerX=0.0,
                centerY=1.0,
                centerZ=0.0,
                radius=1.0,
                show=True,
            )
            self.bondfilters = [bf1, bf2]
        else:
            self.bondfilters = bf

        if bc == "":
            bc1 = BoundaryConditions(
                id=1,
                Name="BC_1",
                NodeSets=None,
                boundarytype="Prescribed Displacement",
                blockId=2,
                coordinate="x",
                value="0*t",
            )
            bc2 = BoundaryConditions(
                id=2,
                Name="BC_2",
                NodeSets=None,
                boundarytype="Prescribed Displacement",
                blockId=3,
                coordinate="x",
                value="0*t",
            )
            bc3 = BoundaryConditions(
                id=3,
                Name="BC_3",
                NodeSets=None,
                boundarytype="Prescribed Displacement",
                blockId=4,
                coordinate="x",
                value="10*t",
            )
            self.bcDict = [bc1, bc2, bc3]
        else:
            self.bcDict = bc

        if solver == "":
            self.solverDict = Solver(
                verbose=False,
                initialTime=0.0,
                finalTime=0.03,
                fixedDt=None,
                solvertype="Verlet",
                safetyFactor=0.95,
                numericalDamping=0.000005,
                peridgimPreconditioner="None",
                nonlinearSolver="Line Search Based",
                numberofLoadSteps=100,
                maxSolverIterations=50,
                relativeTolerance=1e-8,
                maxAgeOfPrec=100,
                directionMethod="Newton",
                newton=Newton(),
                lineSearchMethod="Polynomial",
                verletSwitch=False,
                verlet=Verlet(),
                stopAfterDamageInitation=False,
                stopBeforeDamageInitation=False,
                adaptivetimeStepping=False,
                adapt=Adapt(),
                filetype="yaml",
            )
        else:
            self.solverDict = solver

        self.damBlock = [""] * numberOfBlocks
        self.damBlock[0] = "PMMADamage"
        self.damBlock[1] = "PMMADamage"
        self.damBlock[2] = "PMMADamage"
        self.damBlock[3] = "PMMADamage"

        self.intBlockId = [""] * numberOfBlocks
        self.matBlock = ["PMMA"] * numberOfBlocks

    def createBoundaryConditionBlock(self, x, y, k):
        k = np.where(
            np.logical_and(
                x > self.xend - self.dx[0] * 3,
                np.logical_and(y < 90 + self.dx[0] * 3, y > 90 - self.dx[0] * 3),
            ),
            2,
            k,
        )
        k = np.where(
            np.logical_and(
                x > self.xend - self.dx[0] * 3,
                np.logical_and(y < -90 + self.dx[0] * 3, y > -90 - self.dx[0] * 3),
            ),
            3,
            k,
        )
        return k

    def createLoadIntroNode(self, x, y, k):
        k = np.where(
            np.logical_and(x < self.dx[0] * 3, np.logical_and(y < 25, y > -25)), 4, k
        )
        return k

    def createModel(self):

        geo = Geometry()

        x, y, z = geo.createPoints(
            coor=[
                self.xbegin,
                self.xend,
                self.ybegin,
                self.yend,
                self.zbegin,
                self.zend,
            ],
            dx=self.dx,
        )

        if len(x) > self.maxNodes:
            return (
                "The number of nodes ("
                + str(len(x))
                + ") is larger than the allowed "
                + str(self.maxNodes)
            )

        if self.ignoreMesh and self.blockDef != "":

            writer = ModelWriter(modelClass=self)
            for idx in range(0, len(self.blockDef)):
                self.blockDef[idx].horizon = self.scal * max([self.dx[0], self.dx[1]])
            blockDef = self.blockDef

            try:
                writer.createFile(blockDef)
            except TypeError as e:
                return str(e)

        else:

            vol = np.zeros(len(x))
            k = np.ones(len(x))
            if self.rot:
                angle_x = np.zeros(len(x))
                angle_y = np.zeros(len(x))
                angle_z = np.zeros(len(x))

            k = self.createBoundaryConditionBlock(x, y, k)
            k = self.createLoadIntroNode(x, y, k)

            vol = np.full_like(x, self.dx[0] * self.dx[1] * self.dx[2])

            writer = ModelWriter(modelClass=self)

            if self.rot:
                model = np.transpose(
                    np.vstack(
                        [
                            x.ravel(),
                            y.ravel(),
                            z.ravel(),
                            k.ravel(),
                            vol.ravel(),
                            angle_x.ravel(),
                            angle_y.ravel(),
                            angle_z.ravel(),
                        ]
                    )
                )
                writer.writeMeshWithAngles(model)
            else:
                model = np.transpose(
                    np.vstack([x.ravel(), y.ravel(), z.ravel(), k.ravel(), vol.ravel()])
                )
                writer.writeMesh(model)
            writer.writeNodeSets(model)

            blockLen = int(max(k))

            writeReturn = self.writeFILE(writer=writer, blockLen=blockLen)

            if writeReturn != 0:
                return writeReturn

        return "Model created"

    def createBlockdef(self, blockLen):
        blockDict = []
        for idx in range(0, blockLen):
            blockDef = Block(
                id=1,
                Name="block_" + str(idx + 1),
                material=self.matBlock[idx],
                damageModel=self.damBlock[idx],
                horizon=self.scal * max([self.dx[0], self.dx[1]]),
                interface=self.intBlockId[idx],
                show=False,
            )
            blockDict.append(blockDef)
        # 3d tbd
        return blockDict

    def writeFILE(self, writer, blockLen):

        if self.blockDef == "":
            blockDef = self.createBlockdef(blockLen)
        else:
            for idx in range(0, len(self.blockDef)):
                self.blockDef[idx].horizon = self.scal * max([self.dx[0], self.dx[1]])
            blockDef = self.blockDef

        try:
            writer.createFile(blockDef)
        except TypeError as e:
            return str(e)
        return 0
