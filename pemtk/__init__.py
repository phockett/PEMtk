__version__ = '0.0.1'


# For module dev testing, include path to module here, otherwise use global installation
import sys

# Import epsproc if installed, otherwise try local copy
try:
    import epsproc as ep

except ImportError:
    print("*** ePSproc not found, setting for local copy.")

    if sys.platform == "win32":
        modPath = r'D:\code\github\ePSproc'  # Win test machine
        winFlag = True
    else:
        modPath = r'/home/femtolab/github/ePSproc/'  # Linux test machine
        winFlag = False

    sys.path.append(modPath)

    import epsproc as ep
