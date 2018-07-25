from hypothesis import given
from hypothesis.strategies import text,lists,booleans

ALPHABET=list("ACGT")

def reverse_complement(x):
    wc = dict(zip(ALPHABET,ALPHABET[::-1]))
    return "".join(map(wc.get,x[::-1]))

def test_huddinge_import():
    import pyhuddinge


def test_huddingepair1():
    import pyhuddinge as ph

    x,y = "AAAAAAA","AAATAAA"

    h_dist = ph.huddinge_distance(x,y) 

    assert h_dist==1


@given(text(ALPHABET,min_size=1,max_size=15),text(ALPHABET,min_size=1,max_size=15), booleans() )
def test_huddingepair_symmetry(x,y,reverse_complements):
    import pyhuddinge as ph

    
    h_dist = ph.huddinge_distance(x,y,reverse_complements)
    
    if reverse_complements:
        assert (h_dist==0 and (x==y or x==reverse_complement(y)) ) or (h_dist>0 and x!=y and x!=reverse_complement(y))
    else:
        assert (h_dist==0 and x==y) or (h_dist>0 and x!=y)
    
    assert h_dist <= max(len(x),len(y)),h_dist

    h_dist_v = ph.huddinge_distance(y,x,reverse_complements)

    assert h_dist_v == h_dist



@given(text(ALPHABET,min_size=1,max_size=15),text(ALPHABET,min_size=1,max_size=15) )
def test_huddingepair_rev_comp(x,y):
    import pyhuddinge as ph

    
    h_dist = ph.huddinge_distance(x,y)
    h_dist_rc = ph.huddinge_distance(x,y,True)

    h_dist_rcx = ph.huddinge_distance(reverse_complement(x),y)
    h_dist_rcy = ph.huddinge_distance(x,reverse_complement(y))
    h_dist_rcxy = ph.huddinge_distance(reverse_complement(x),reverse_complement(y))


    assert h_dist == h_dist_rcxy
    assert h_dist_rc <= h_dist
    assert h_dist_rc == min(h_dist_rcx, h_dist)
    assert h_dist_rcx == h_dist_rcy


@given(text(ALPHABET,min_size=1,max_size=15),text(ALPHABET,min_size=1,max_size=15) )
def test_huddingepair_rc_closer(x,y):
    import pyhuddinge as ph

    
    h_dist = ph.huddinge_distance(x,y,False)
    h_dist_rc = ph.huddinge_distance(y,x,True)

    assert (h_dist>=h_dist_rc)
    assert (h_dist==0 and x==y ) or h_dist>0
    assert h_dist <= max(len(x),len(y)),h_dist

    





def test_all_huddinge_pairs1():
    import pyhuddinge as ph

 
    kmers = ["AAGGGAA","AAAGCAA","AAACCAA","AAACAAA"]
    import itertools as it

    kmers = ["".join(x) for x in it.product("ACGT",repeat=6)]

    D = ph.all_pairs_huddinge_distance(kmers)

    assert len(D) == len(kmers)*(len(kmers)-1)/2
    assert min(D)>0


@given(lists(text(ALPHABET,min_size=1,max_size=15),max_size=30,unique=True),booleans())
def test_all_huddinge_pairs_concordance(kmers,reverse_complements):
    import pyhuddinge as ph

    D = ph.all_pairs_huddinge_distance(kmers,reverse_complements)

    assert len(D) == len(kmers)*(len(kmers)-1)/2
    assert len(D)==0 or min(D)>0 or reverse_complements

    Dsq = ph.squareform(D,kmers)

    for i in Dsq.index:
        for j in Dsq.columns:
            if i<j:
                assert Dsq.loc[i,j] == ph.huddinge_distance(i,j,reverse_complements)

    
 

@given(lists(text(ALPHABET,min_size=1,max_size=15),max_size=30,unique=True))
def test_all_huddinge_pairs_shuffle(kmers):
    import pyhuddinge as ph

 
    D1 = ph.all_pairs_huddinge_distance(kmers)
    D2 = ph.all_pairs_huddinge_distance(kmers[::-1])


    Dsq1 = ph.squareform(D1,kmers)
    Dsq2 = ph.squareform(D2,kmers[::-1])
    reordered2 = Dsq2[Dsq1.columns].loc[Dsq1.index]
    
    print(Dsq1)
    print(reordered2)
    assert (Dsq1==reordered2).all().all()
    
 