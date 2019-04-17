import ctypes

hlib = ctypes.cdll.LoadLibrary("libpyhuddingec.so")
hlib.pyhuddinge_distance.argtypes = [ctypes.c_char_p,ctypes.c_char_p,ctypes.c_bool]
ALPHABET = set("ACGTacgt")

def huddinge_distance(x,y,reverse_complements=False):
    """Computes Huddinge distance between two kmers x and y 
    and minimum with their reverse compelemnt if reverse_complements==True.
    Gaps are denoted with lower case 'n'"""

    import six
    
    x = x.upper().strip("N")
    y = y.upper().strip("N")

    if not set(x+y).issubset(ALPHABET):
        raise ValueError("Unrecognized character {} in input.".format(set(x+y)-ALPHABET))

    if len(x)==0 or len(y)==0: # Special case which would crash distance calculation
        return 0

    if six.PY3:
        x=x.encode('ascii')
        y=y.encode('ascii')

    h_dist = hlib.pyhuddinge_distance(x,y,reverse_complements)

    return h_dist


def all_pairs_huddinge_distance(kmers,reverse_complements=False):
    """Compute Huddinge distance between all pairs of kmers. 
    Return value is compressed distance matrix which can be 
    expanded with scipy.spatial.squareform"""

    kmers = [x.upper().strip("N") for x in kmers]

    all_chars = set()
    all_chars.update(*[set(x) for x in kmers])
    if not all_chars.issubset(ALPHABET):
        raise ValueError("Unrecognized character {} in input.".format(all_chars-ALPHABET))
    if min(len(x) for x in kmers)==0:
        raise ValueError("All input kmers must have at least one defined character.")

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

    hlib.pyhuddinge_all_pairs_distance.argtypes = [array_1d_int,ctypes.c_int,array_1d_char_p, ctypes.c_bool]
    hlib.pyhuddinge_all_pairs_distance(pair_dists, 
            N, arr,reverse_complements)

    return pair_dists




def huddinge_neighbours(kmer,max_dist=1,min_kmer_len=None,max_kmer_len=None,max_gapped_len=None):
    """Compute Huddinge distance between all pairs of kmers. 
    Return value is compressed distance matrix which can be 
    expanded with scipy.spatial.squareform"""

    import six
    kmer = kmer.upper()


    if six.PY3:
        kmer = kmer.encode('ascii')

    assert max_dist>0
    
    if min_kmer_len is None:
        min_kmer_len = len(kmer)
    if max_kmer_len is None:
        max_kmer_len = min_kmer_len

    if max_gapped_len is None:
        max_gapped_len = max_kmer_len

    assert min_kmer_len>0
    assert max_kmer_len>0
    assert min_kmer_len <= max_kmer_len
    assert max_gapped_len>=max_kmer_len

    hlib.pyhuddinge_neighbourhood.argtypes = [ctypes.c_char_p,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int]

    hlib.pyhuddinge_neighbourhood.restype = ctypes.py_object

    neighbours = hlib.pyhuddinge_neighbourhood(kmer,max_dist,min_kmer_len,max_kmer_len,max_gapped_len)

    return neighbours


def squareform(D,kmers):
    "Convert all pairs distance matrix to a dataframe with index and column names"
    from scipy.spatial.distance import squareform
    import pandas as pd

    if len(kmers)>0:
        D = pd.DataFrame(squareform(D),index=kmers,columns=kmers)
    else:
        D = pd.DataFrame()

    return D
