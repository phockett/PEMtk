"""
Env util routines for PEMtk

"""

def isnotebook():
    """
    Check if code is running in Jupyter Notebook.

    Taken verbatim from https://exceptionshub.com/how-can-i-check-if-code-is-executed-in-the-ipython-notebook.html
    Might be a better/more robust way to do this?

    """

    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter
