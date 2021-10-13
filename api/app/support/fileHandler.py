import shutil
import filecmp
import paramiko
import os

class fileHandler(object):
    
    def getRemoteModelPath(Cluster, username, ModelName):
        
        if Cluster=='FA-Cluster':
            remotepath = './PeridigmJobs/apiModels/' + os.path.join(username, ModelName)
        
        elif Cluster=='Cara':
            remotepath = './PeridigmJobs/apiModels/' + os.path.join(username, ModelName)
        return remotepath

    def getRemoteUserPath(Cluster, username):
        
        if Cluster=='FA-Cluster':
            remotepath = './PeridigmJobs/apiModels/' + username
        
        elif Cluster=='Cara':
            remotepath = './PeridigmJobs/apiModels/' + username
        return remotepath

    def getUserPath(Cluster, username, ModelName):
        
        if Cluster=='FA-Cluster':
            userpath = './PeridigmJobs/apiModels/' + username
        
        elif Cluster=='Cara':
            userpath = './PeridigmJobs/apiModels/' + username
        return userpath

    def copyModelToCluster(username, ModelName, Cluster):
        
        if Cluster=='None':
            localpath = './Output/' + os.path.join(username, ModelName)
            remotepath = './peridigmJobs/' + os.path.join(username, ModelName)
            if not os.path.exists(remotepath):
                os.makedirs(remotepath)
                # os.chown(remotepath, 'test')
            if not os.path.exists(localpath):
                return ModelName + ' has not been created yet'
            for root, dirs, files in os.walk(localpath):
                if len(files)==0:
                    return ModelName + ' has not been created yet'

                for name in files:
                    shutil.copy(os.path.join(root,name), os.path.join(remotepath,name))
                    # os.chmod(os.path.join(remotepath,name), 0o0777)
                    # os.chown(os.path.join(remotepath,name), 'test')
            return ModelName + ' has been copied'

        else:       
            
            remotepath = fileHandler.getRemoteModelPath(Cluster, username, ModelName)
            userpath = fileHandler.getUserPath(Cluster, username, ModelName) 
            ssh, sftp = fileHandler.sftpToCluster(Cluster, username)

            try:
                sftp.chdir(remotepath)  # Test if remote_path exists
            except IOError:
                sftp.mkdir(userpath)
                sftp.mkdir(remotepath)  # Create remote_path
                sftp.chdir(remotepath)
            for root, dirs, files in os.walk('./Output/' + os.path.join(username, ModelName)):
                if len(files)==0:
                    return ModelName + ' has not been created yet'
                for name in files:
                    sftp.put(os.path.join(root,name), name)

            sftp.close()
            ssh.close()
            
            return ModelName + ' has been copied to Cluster'

    def copyResultsFromCluster(username, ModelName, Cluster, allData):
        resultpath = './Results/' + os.path.join(username, ModelName)
        if not os.path.exists(resultpath):
            os.makedirs(resultpath)

        if Cluster=='None':
            remotepath = './peridigmJobs/' + os.path.join(username, ModelName)
            for root, dirs, files in os.walk(remotepath):
                if len(files)==0:
                    return ModelName + ' results can be found on ' + Cluster

                for filename in files:
                    if(allData or '.e' in filename):
                        if(os.path.exists(os.path.join(resultpath, filename))):
                            if(filecmp.cmp(os.path.join(remotepath, filename),os.path.join(resultpath, filename))==False):
                                shutil.copy(os.path.join(remotepath, filename), os.path.join(resultpath,filename))
                        else:
                            shutil.copy(os.path.join(remotepath, filename), os.path.join(resultpath,filename))
                    # os.chmod(os.path.join(remotepath,name), 0o0777)
                    # os.chown(os.path.join(remotepath,name), 'test')
            # return ModelName + ' has been copied'

        else:
            remotepath = fileHandler.getRemoteModelPath(Cluster, username, ModelName)
            ssh, sftp = fileHandler.sftpToCluster(Cluster, username)
            try:
                for filename in sftp.listdir(remotepath):
                    if(allData or '.e' in filename):
                        if(os.path.exists(os.path.join(resultpath, filename))):
                            remoteInfo = sftp.stat(os.path.join(remotepath, filename))
                            remoteSize = remoteInfo.st_size
                            localSize = os.path.getsize(os.path.join(resultpath, filename))
                            print('compare ' + filename + ' remoteSize: ' + remoteSize + ', localsize: ' + localSize)
                            if(abs(remoteSize-localSize)>5):
                                sftp.get(os.path.join(remotepath, filename), os.path.join(resultpath, filename))
                        else:
                            sftp.get(os.path.join(remotepath, filename), os.path.join(resultpath, filename))
            except:
                return ModelName + ' results can be found on ' + Cluster
            sftp.close()
            ssh.close()
            # return ModelName + ' has been copied'

    def sftpToCluster(Cluster, username):
        
        if Cluster=='FA-Cluster':
            username='f_peridi'
            server='129.247.54.37'
            keypath = 'id_rsa_cluster'
        
        elif Cluster=='Cara':
            username='hess_ja'
            server='cara.dlr.de'
            keypath = 'id_rsa_cara'

        ssh = paramiko.SSHClient() 
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(server, username=username, allow_agent=False, key_filename=keypath)
        sftp = ssh.open_sftp()
        return ssh, sftp

    def sshToCluster(Cluster, username):
        
        if Cluster=='FA-Cluster':
            username='f_peridi'
            server='129.247.54.37'
            keypath = 'id_rsa_cluster'
        
        elif Cluster=='Cara':
            username='hess_ja'
            server='cara.dlr.de'
            keypath = 'id_rsa_cara'

        ssh = paramiko.SSHClient() 
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(server, username=username, allow_agent=False, key_filename=keypath)
        return ssh
