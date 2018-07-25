import ctypes
hlib = ctypes.cdll.LoadLibrary("src/libpyhuddingec.so")
hlib.pyhuddinge_distance.argtypes = [ctypes.c_char_p,ctypes.c_char_p]

def huddinge_distance(x,y):
    "Computes Huddinge distance between two kmers x and y"

    import six
    if six.PY3:
        x=x.encode('ascii')
        y=y.encode('ascii')

    h_dist = hlib.pyhuddinge_distance(x,y)

    return h_dist


def all_pairs_huddinge_distance(kmers):
    """Compute Huddinge distance between all pairs of kmers. 
    Return value is compressed distance matrix which can be 
    expanded with scipy.spatial.squareform"""

    import six
    if six.PY3:
        kmers = [ x.encode('ascii') for x in kmers]
    import numpy as np
    import numpy.ctypeslib as npct
    array_1d_int = npct.ndpointer(dtype=np.int8, ndim=1, flags='CONTIGUOUS')

    N = len(kmers)
    d_count = int((N*(N-1))/2)
    #print(N,type(N),d_count,type(d_count))
    pair_dists = np.empty(d_count,dtype=np.int8)
    
    array_1d_char_p = ctypes.c_char_p * N
    arr = array_1d_char_p()
    arr[:] = kmers

    hlib.pyhuddinge_all_pairs_distance.argtypes = [array_1d_int,ctypes.c_int,array_1d_char_p]
    hlib.pyhuddinge_all_pairs_distance(pair_dists, 
            N, arr)

    return pair_dists

def squareform(D,kmers):
    "Convert all pairs distance matrix to a dataframe with index and column names"
    from scipy.spatial.distance import squareform
    import pandas as pd

    if len(kmers)>0:
        D = pd.DataFrame(squareform(D),index=kmers,columns=kmers)
    else:
        D = pd.DataFrame()

    return D
