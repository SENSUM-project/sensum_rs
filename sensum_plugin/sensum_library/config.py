import sys
import os
#import site

if os.name == "posix":
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    os.environ['ITK_AUTOLOAD_PATH'] = "/usr/local/lib/otb/applications"
    #os.chdir('/home/gale/Izmir')
else:
    if os.path.isdir("C:/OSGeo4W64"):
        sys.path.append("C:/OSGeo4W64/apps/Python27/Lib/site-packages")
        sys.path.append("C:/OSGeo4W64/apps/orfeotoolbox/python")
        os.environ["PATH"] = os.environ["PATH"] + ";C:/OSGeo4W64/bin"
    elif os.path.isdir("C:/OSGeo4W"):
        sys.path.append("C:/OSGeo4W64/apps/orfeotoolbox/python")
        sys.path.append("C:/OSGeo4W/apps/Python27/Lib/site-packages")
        os.environ["PATH"] = os.environ["PATH"] + ";C:/OSGeo4W/bin"
    #print site.getsitepackages()[0]
