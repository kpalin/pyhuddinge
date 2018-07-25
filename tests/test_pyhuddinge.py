from hypothesis import given
from hypothesis.strategies import text,lists

ALPHABET=list("ACGT")

def test_huddinge_import():
    import pyhuddinge


def test_huddingepair1():
    import pyhuddinge as ph

    x,y = "AAAAAAA","AAATAAA"

    h_dist = ph.huddinge_distance(x,y) 

    assert h_dist==1


@given(text(ALPHABET,min_size=1,max_size=15),text(ALPHABET,min_size=1,max_size=15))
def test_huddingepair2(x,y):
    import pyhuddinge as ph

    
    h_dist = ph.huddinge_distance(x,y)

    assert h_dist>0 or x==y
    assert h_dist <= max(len(x),len(y)),h_dist


def test_all_huddinge_pairs1():
    import pyhuddinge as ph

 
    kmers = ["AAGGGAA","AAAGCAA","AAACCAA","AAACAAA"]
    import itertools as it

    kmers = ["".join(x) for x in it.product("ACGT",repeat=6)]

    D = ph.all_pairs_huddinge_distance(kmers)

    assert len(D) == len(kmers)*(len(kmers)-1)/2
    assert min(D)>0


@given(lists(text(ALPHABET,min_size=1,max_size=15),max_size=60,unique=True))
def test_all_huddinge_pairs_concordance(kmers):
    import pyhuddinge as ph

 
    print(kmers)
    D = ph.all_pairs_huddinge_distance(kmers)
    print(D)
    assert len(D) == len(kmers)*(len(kmers)-1)/2
    assert len(D)==0 or min(D)>0

    Dsq = ph.squareform(D,kmers)

    for i in Dsq.index:
        for j in Dsq.columns:
            if i<j:
                assert Dsq.loc[i,j] == ph.huddinge_distance(i,j)

    
 

@given(lists(text(ALPHABET,min_size=1,max_size=15),max_size=60,unique=True))
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
    
 
