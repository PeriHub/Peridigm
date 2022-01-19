from pydantic import BaseModel
from typing import Optional
from typing import List
from enum import Enum

class Model(BaseModel):
    ownModel: bool
    translated: bool
    length: float
    width: float
    height: float
    height2: float
    structured: bool
    discretization: int
    horizon: float
    twoDimensional: bool
    rotatedAngles: bool
    angles: List[float]
class Properties(BaseModel):
    id: int
    Name: str
    value: Optional[float] = None
class Parameter(BaseModel):
    Name: str
    value: Optional[float] = None
class Material(BaseModel):
    id: int
    Name: str
    MatType: str
    density: float
    bulkModulus: Optional[float] = None
    shearModulus: Optional[float] = None
    youngsModulus: float
    poissonsRatio: float
    tensionSeparation: bool
    nonLinear: bool
    planeStress: bool
    materialSymmetry: str
    stabilizatonType: str
    thickness: float
    hourglassCoefficient: float
    actualHorizon: Optional[float] = None
    yieldStress: Optional[float] = None
    Parameter: List[Parameter]
    Properties: List[Properties]
class Damage(BaseModel):
    id: int
    Name: str
    damageModel: str
    criticalStretch: float
    criticalEnergy: float
    interblockdamageEnergy: float
    planeStress: bool
    onlyTension: bool
    detachedNodesCheck: bool
    thickness: float
    hourglassCoefficient: float
    stabilizatonType: str
class Block(BaseModel):
    id: int
    Name: str
    material: str
    damageModel: str
    horizon: Optional[float] = None
    interface: Optional[str] = None
    show: bool
class BoundaryConditions(BaseModel):
    id: int
    Name: str
    nodeSet: Optional[str] = None
    boundarytype: str
    blockId: int
    coordinate: str
    value: str
class Compute(BaseModel):
    id: int
    Name: str
    variable: str
    calculationType: str
    blockName: str
class Output(BaseModel):
    id: int
    Name: str
    Displacement: bool
    Force: bool
    Damage: bool
    Partial_Stress: bool
    External_Force: bool
    External_Displacement: bool
    Number_Of_Neighbors: bool
    Frequency: int
    InitStep: int
class Newton(BaseModel):
    jacobianOperator: str
    preconditioner: str
class Verlet(BaseModel):
    safetyFactor: float
    numericalDamping: float
    outputFrequency: int
class Adapt(BaseModel):
    stableStepDifference: int
    maximumBondDifference: int
    stableBondDifference: int
class FileType(str, Enum):
    yaml = "yaml"
    xml = "xml"
class Solver(BaseModel):
    verbose: bool
    initialTime: float
    finalTime: float
    solvertype: str
    safetyFactor: float
    numericalDamping: float
    peridgimPreconditioner: str
    nonlinearSolver: str
    numberofLoadSteps: int
    maxSolverIterations: int
    relativeTolerance: float
    maxAgeOfPrec: float
    directionMethod: str
    newton: Newton
    lineSearchMethod: str
    verletSwitch: bool
    verlet: Verlet
    stopAfterDamageInitation: bool
    stopBeforeDamageInitation: bool
    adaptivetimeStepping: bool
    adapt: Adapt
    filetype: FileType
class Job(BaseModel):
    cluster: str
    tasks: int
    time: str
    account: int
class RunData(BaseModel):
    job: Job
    materials: List[Material]
    outputs: List[Output]
    class Config:
        schema_extra = {
            "example": {
                "job": {
                    "cluster": 'None',
                    "tasks": 1,
                    "time": '40:00:00',
                    "account": '2263032',
                },
                "materials": [
                    {
                    "id": 1,
                    "Name": "PMMA",
                    "MatType": "Linear Elastic Correspondence",
                    "density": "1.4e5",
                    "bulkModulus": None,
                    "shearModulus": None,
                    "youngsModulus": "2.997e9",
                    "poissonsRatio": 0.3,
                    "tensionSeparation": False,
                    "nonLinear": True,
                    "planeStress": True,
                    "materialSymmetry": "Isotropic",
                    "stabilizatonType": "Global Stiffness",
                    "thickness": "0.01",
                    "hourglassCoefficient": 1,
                    "actualHorizon": None,
                    "yieldStress": None,
                    "Parameter": {
                        "C11": None,
                        "C12": None,
                        "C13": None,
                        "C14": None,
                        "C15": None,
                        "C16": None,
                        "C22": None,
                        "C23": None,
                        "C24": None,
                        "C25": None,
                        "C26": None,
                        "C33": None,
                        "C34": None,
                        "C35": None,
                        "C36": None,
                        "C44": None,
                        "C45": None,
                        "C46": None,
                        "C55": None,
                        "C56": None,
                        "C66": None
                    },
                    "Properties": [
                        {
                        "id": 1,
                        "Name": "Prop_1",
                        "value": None
                        }
                    ]
                    },
                    {
                    "id": 2,
                    "Name": "PMMAElast",
                    "MatType": "Linear Elastic Correspondence",
                    "density": "1.4e5",
                    "bulkModulus": None,
                    "shearModulus": None,
                    "youngsModulus": "2.997e9",
                    "poissonsRatio": "0.3",
                    "tensionSeparation": False,
                    "nonLinear": True,
                    "planeStress": True,
                    "materialSymmetry": "Isotropic",
                    "stabilizatonType": "Global Stiffness",
                    "thickness": "0.01",
                    "hourglassCoefficient": "1",
                    "actualHorizon": None,
                    "yieldStress": None,
                    "Parameter": {
                        "C11": None,
                        "C12": None,
                        "C13": None,
                        "C14": None,
                        "C15": None,
                        "C16": None,
                        "C22": None,
                        "C23": None,
                        "C24": None,
                        "C25": None,
                        "C26": None,
                        "C33": None,
                        "C34": None,
                        "C35": None,
                        "C36": None,
                        "C44": None,
                        "C45": None,
                        "C46": None,
                        "C55": None,
                        "C56": None,
                        "C66": None
                    },
                    "Properties": [
                        {
                        "id": 1,
                        "Name": "Prop_1",
                        "value": None
                        }
                    ]
                    }
                ],
                "outputs": [
                    {
                    "id": 1,
                    "Name": "Output1",
                    "Displacement": True,
                    "Force": True,
                    "Damage": True,
                    "Partial_Stress": True,
                    "External_Force": True,
                    "External_Displacement": True,
                    "Number_Of_Neighbors": False,
                    "Frequency": "100",
                    "InitStep": 0
                    }
                ],
            }
        }
class Data(BaseModel):
    model: Model
    materials: List[Material]
    damages: List[Damage]
    blocks: List[Block]
    boundaryConditions: List[BoundaryConditions]
    computes: List[Compute]
    outputs: List[Output]
    solver: Solver

    class Config:
        schema_extra = {
            "example": {
                "model": {
                    "ownModel": False,
                    "translated": False,
                    "length": 0.15,
                    "width": 0.001,
                    "height": 0.02,
                    "height2": 0.01,
                    "structured": True,
                    "discretization": 21,
                    "horizon": 0.01,
                    "twoDimensional": True,
                    "rotatedAngles": False,
                    "angles": [
                    0,
                    0
                    ]
                },
                "materials": [
                    {
                    "id": 1,
                    "Name": "PMMA",
                    "MatType": "Linear Elastic Correspondence",
                    "density": "1.4e5",
                    "bulkModulus": None,
                    "shearModulus": None,
                    "youngsModulus": "2.997e9",
                    "poissonsRatio": 0.3,
                    "tensionSeparation": False,
                    "nonLinear": True,
                    "planeStress": True,
                    "materialSymmetry": "Isotropic",
                    "stabilizatonType": "Global Stiffness",
                    "thickness": "0.01",
                    "hourglassCoefficient": 1,
                    "actualHorizon": None,
                    "yieldStress": None,
                    "Parameter": {
                        "C11": None,
                        "C12": None,
                        "C13": None,
                        "C14": None,
                        "C15": None,
                        "C16": None,
                        "C22": None,
                        "C23": None,
                        "C24": None,
                        "C25": None,
                        "C26": None,
                        "C33": None,
                        "C34": None,
                        "C35": None,
                        "C36": None,
                        "C44": None,
                        "C45": None,
                        "C46": None,
                        "C55": None,
                        "C56": None,
                        "C66": None
                    },
                    "Properties": [
                        {
                        "id": 1,
                        "Name": "Prop_1",
                        "value": None
                        }
                    ]
                    },
                    {
                    "id": 2,
                    "Name": "PMMAElast",
                    "MatType": "Linear Elastic Correspondence",
                    "density": "1.4e5",
                    "bulkModulus": None,
                    "shearModulus": None,
                    "youngsModulus": "2.997e9",
                    "poissonsRatio": "0.3",
                    "tensionSeparation": False,
                    "nonLinear": True,
                    "planeStress": True,
                    "materialSymmetry": "Isotropic",
                    "stabilizatonType": "Global Stiffness",
                    "thickness": "0.01",
                    "hourglassCoefficient": "1",
                    "actualHorizon": None,
                    "yieldStress": None,
                    "Parameter": {
                        "C11": None,
                        "C12": None,
                        "C13": None,
                        "C14": None,
                        "C15": None,
                        "C16": None,
                        "C22": None,
                        "C23": None,
                        "C24": None,
                        "C25": None,
                        "C26": None,
                        "C33": None,
                        "C34": None,
                        "C35": None,
                        "C36": None,
                        "C44": None,
                        "C45": None,
                        "C46": None,
                        "C55": None,
                        "C56": None,
                        "C66": None
                    },
                    "Properties": [
                        {
                        "id": 1,
                        "Name": "Prop_1",
                        "value": None
                        }
                    ]
                    }
                ],
                "damages": [
                    {
                    "id": 1,
                    "Name": "PMMADamage",
                    "damageModel": "Critical Energy Correspondence",
                    "criticalStretch": 10,
                    "criticalEnergy": "10.1",
                    "interblockdamageEnergy": "0.01",
                    "planeStress": True,
                    "onlyTension": False,
                    "detachedNodesCheck": True,
                    "thickness": 10,
                    "hourglassCoefficient": 1,
                    "stabilizatonType": "Global Stiffness"
                    }
                ],
                "blocks": [
                    {
                    "id": 1,
                    "Name": "block_1",
                    "material": "PMMAElast",
                    "damageModel": "",
                    "interface": "",
                    "show": True
                    },
                    {
                    "id": 2,
                    "Name": "block_2",
                    "material": "PMMAElast",
                    "damageModel": "",
                    "interface": "",
                    "show": True
                    },
                    {
                    "id": 3,
                    "Name": "block_3",
                    "material": "PMMA",
                    "damageModel": "PMMADamage",
                    "interface": "",
                    "show": True
                    },
                    {
                    "id": 4,
                    "Name": "block_4",
                    "material": "PMMAElast",
                    "damageModel": "",
                    "interface": "",
                    "show": True
                    },
                    {
                    "id": 5,
                    "Name": "block_5",
                    "material": "PMMAElast",
                    "damageModel": "",
                    "interface": "",
                    "show": True
                    }
                ],
                "boundaryConditions": [
                    {
                    "id": 1,
                    "Name": "BC_1",
                    "nodeSet": None,
                    "boundarytype": "Prescribed Displacement",
                    "blockId": 1,
                    "coordinate": "x",
                    "value": "0*t"
                    },
                    {
                    "id": 2,
                    "Name": "BC_2",
                    "nodeSet": None,
                    "boundarytype": "Prescribed Displacement",
                    "blockId": 5,
                    "coordinate": "x",
                    "value": "0.05*t"
                    }
                ],
                "computes": [
                    {
                    "id": 1,
                    "Name": "External_Displacement",
                    "variable": "Displacement",
                    "calculationType": "Maximum",
                    "blockName": "block_5"
                    },
                    {
                    "id": 2,
                    "Name": "External_Force",
                    "variable": "Force",
                    "calculationType": "Sum",
                    "blockName": "block_5"
                    }
                ],
                "outputs": [
                    {
                    "id": 1,
                    "Name": "Output1",
                    "Displacement": True,
                    "Force": True,
                    "Damage": True,
                    "Partial_Stress": True,
                    "External_Force": True,
                    "External_Displacement": True,
                    "Number_Of_Neighbors": False,
                    "Frequency": "100",
                    "InitStep": 0
                    }
                ],
                "solver": {
                    "verbose": False,
                    "initialTime": 0,
                    "finalTime": "0.0075",
                    "solvertype": "Verlet",
                    "safetyFactor": "0.9",
                    "numericalDamping": "0.0005",
                    "peridgimPreconditioner": "None",
                    "nonlinearSolver": "Line Search Based",
                    "numberofLoadSteps": 100,
                    "maxSolverIterations": 50,
                    "relativeTolerance": 1e-8,
                    "maxAgeOfPrec": 100,
                    "directionMethod": "Newton",
                    "newton": {
                    "jacobianOperator": "Matrix-Free",
                    "preconditioner": "None"
                    },
                    "lineSearchMethod": "Polynomial",
                    "verletSwitch": True,
                    "verlet": {
                    "safetyFactor": 0.95,
                    "numericalDamping": 0.000005,
                    "outputFrequency": 7500
                    },
                    "stopAfterDamageInitation": False,
                    "stopBeforeDamageInitation": False,
                    "adaptivetimeStepping": False,
                    "adapt": {
                    "stableStepDifference": 4,
                    "maximumBondDifference": 4,
                    "stableBondDifference": 1
                    },
                    "filetype": "yaml"
                }
            }
        }