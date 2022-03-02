import numpy as np
import os


class SbatchCreator(object):
    def __init__(
        self,
        filename="Model",
        filetype="yaml",
        remotepath="",
        output="",
        job="",
        username="",
        usermail="",
    ):
        self.tasks = job.tasks
        self.time = job.time
        self.filename = filename
        self.filetype = filetype
        self.user = "f_peridi"
        self.account = job.account
        self.mail = usermail
        self.outputDict = output
        self.remotepath = remotepath

    def createSbatch(self):
        nodes = -(-int(self.tasks) // 64)
        string = "#!/bin/bash" + "\n"
        string += "#SBATCH --nodes=" + str(nodes) + "\n"
        string += "#SBATCH --hint=multithread" + "\n"
        string += "#SBATCH --time=" + self.time + "\n"
        string += "#SBATCH --account=" + str(self.account) + "\n"
        string += "#SBATCH --output=simulation-%j.log" + "\n"
        string += "#SBATCH --error=simulation-%j.log" + "\n"
        string += "#SBATCH --open-mode=append" + "\n"

        if self.mail != "":
            string += "#SBATCH --mail-user=" + self.mail + "\n"
            string += "#SBATCH --mail-type ALL" + "\n"

        string += "module purge" + "\n"
        string += "module load GCC" + "\n"
        string += "module load OpenMPI" + "\n"
        string += "module load OpenBLAS" + "\n"
        string += "module load netCDF" + "\n"
        string += "module load ScaLAPACK" + "\n"
        string += "module load GCCcore/8.3.0" + "\n"
        string += "module load HDF5" + "\n"
        string += "module load GCCcore/10.2.0" + "\n"
        string += "module load CMake/3.18.4" + "\n"
        string += (
            "export PATH=$PATH:/home/" + self.user + "/software/trilinos/bin" + "\n"
        )

        string += (
            "srun --ntasks="
            + str(self.tasks)
            + " /home/"
            + self.user
            + "/software/peridigm/bin/Peridigm "
            + self.filename
            + ".yaml"
            + "\n"
        )

        for out in self.outputDict:
            string += (
                "python /home/"
                + self.user
                + "/peridigm/src/scripts/MergeFiles.py "
                + out.Name
                + " "
                + str(self.tasks)
                + "\n"
            )

        string += "rm *.e." + str(self.tasks) + ".*" + "\n"
        string += "rm Output*.log" + "\n"

        return string

    def createSh(self):
        # if(self.tasks==1):
        string = "#!/bin/sh" + "\n"
        string += (
            "/Peridigm/build/src/Peridigm "
            + self.filename
            + "."
            + self.filetype
            + "& echo $! > pid.txt \n"
        )
        string += "pid=`cat pid.txt` \n"
        string += "tail --pid=$pid -f /dev/null \n"
        string += "rm pid.txt \n"
        # else:
        #     string='#!/bin/sh' + '\n'
        #     string += 'export PATH=$PATH:/usr/local/trilinos/bin \n'
        #     string += 'mpirun --allow-run-as-root -np ' + self.tasks + ' /peridigm/build/src/Peridigm ' + self.filename + '.' + self.filetype + '& echo $! > pid.txt \n'
        #     string += 'pid=`cat pid.txt` \n'
        #     string += 'tail --pid=$pid -f /dev/null \n'
        #     for out in self.outputDict:
        #         string += 'python /peridigm/scripts/MergeFiles.py ' + out.Name + ' ' + str(self.tasks) + '\n'
        #     string += 'rm *.e.' + str(self.tasks) + '.*' +  '\n'
        #     string += 'rm Output*.log \n'
        #     string += 'rm pid.txt \n'

        return string
