from __future__ import print_function
try:
    import numpy as np
    import scipy.io as sio
    print('<CALLRESULT>')
    print('TRUE')
    exit()    
except ImportError:
    print('<CALLRESULT>')
    print('FALSE')
    exit()    
