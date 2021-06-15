import os, sys
#insert smetana sorce folder to path
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))

from modelGeneratorControl import ModelControl

if __name__ == '__main__':
    
    kwargs = {}
    mycontrol = ModelControl()
    try:        
       
        mycontrol.run(**kwargs)

    except:
        mycontrol.endRunOnError()
 
        raise
    


    

