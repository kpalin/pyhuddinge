from hypothesis import given,settings
from hypothesis.strategies import text,lists,booleans,integers
import pytest

ALPHABET=list("acgt")

def reverse_complement(x):
    

    wc = dict(zip("acgtACGT","tgcaTGCA"))
    return "".join(wc.get(y,y) for y in x[::-1])


@given(text(ALPHABET,min_size=1,max_size=15))
def test_rev_complement_symmetry(x):
    assert x == reverse_complement(reverse_complement(x))


@given(text([x.upper() for x in ALPHABET],min_size=1,max_size=15))
def test_rev_complement_symmetry_upper(x):
    assert x == reverse_complement(reverse_complement(x))

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
        assert (h_dist==0 and (x.strip("n")==y.strip("n") or x.strip("n")==reverse_complement(y).strip("n")) or min(len(x.strip("n")),len(y.strip("n")))==0) or (h_dist>0 and x!=y and x!=reverse_complement(y))
    else:
        assert (h_dist==0 and x.strip("n")==y.strip("n") or min(len(x.strip("n")),len(y.strip("n")))==0) or (h_dist>0 and x.strip("n")!=y.strip("n"))
    
    assert h_dist <= max(len(x),len(y)),h_dist

    h_dist_v = ph.huddinge_distance(y,x,reverse_complements)

    assert h_dist_v == h_dist

@given(text([x for x in ALPHABET if x.islower()],min_size=1,max_size=15),
        text([x for x in ALPHABET if x.islower()],min_size=1,max_size=15) )
def test_huddinge_upperlower(x,y):
    import pyhuddinge as ph

    
    h_dist_ll = ph.huddinge_distance(x,y)
    h_dist_lu = ph.huddinge_distance(x,y.upper())
    h_dist_uu = ph.huddinge_distance(x.upper(),y.upper())
    h_dist_ul = ph.huddinge_distance(x.upper(),y)
    assert h_dist_ll==h_dist_lu
    assert h_dist_ll==h_dist_ul
    assert h_dist_ll==h_dist_uu


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

@given(text(min_size=1,max_size=15).filter(lambda x:~(set(x).issubset(ALPHABET))),
        text(ALPHABET,min_size=1,max_size=15).filter(lambda x:~(set(x).issubset(ALPHABET))) )
def test_huddingepair_nonbase(x,y):
    import pyhuddinge as ph

    with pytest.raises(ValueError):
        h_dist = ph.huddinge_distance(x,y)

 
@settings(deadline=1000)
@given(lists(text(min_size=1,max_size=15).filter(lambda x:~(set(x).issubset(ALPHABET))),min_size=1,
    max_size=20,unique=True))
def test_all_huddinge_pairs_nonbase(kmers):
    import pyhuddinge as ph

 
    with pytest.raises(ValueError):
        D1 = ph.all_pairs_huddinge_distance(kmers)


def test_all_huddinge_pairs_empty():
    import pyhuddinge as ph
    kmers=[]
    with pytest.raises(ValueError):
        D1 = ph.all_pairs_huddinge_distance(kmers)
   
def test_all_huddinge_pairs_puregap():
    import pyhuddinge as ph
    kmers=['nnnnnnnnnnnn','nnn']
    with pytest.raises(ValueError):
        D1 = ph.all_pairs_huddinge_distance(kmers)
  



@given(text(ALPHABET,min_size=1,max_size=15),text(ALPHABET,min_size=1,max_size=15) )
def test_huddingepair_rc_closer(x,y):
    import pyhuddinge as ph

    
    h_dist = ph.huddinge_distance(x,y,False)
    h_dist_rc = ph.huddinge_distance(y,x,True)

    assert (h_dist>=h_dist_rc)
    
    if min(len(x.strip("n")),len(y.strip("n"))) >0:
        assert (h_dist==0 and x.strip("n")==y.strip("n") ) or h_dist>0
        assert h_dist <= max(len(x.strip("n")),len(y.strip("n"))),h_dist

    


@pytest.mark.skip("MODER version does not work with gaps.")
def test_gapkmers_PMC4362205_figS1B_ACnGT():
    import pyhuddinge as ph
    A="ACnGT"
    ex = {"ACnCT":1,
        "ACTG":1,
        "ACnnTA":1,
        "ACnG":1,
        "ACnnT":1,
        "ACnGTA":1,
        "ACAGT":1}
    for B,d in ex.items():
        print(A,B)
        assert ph.huddinge_distance(A,B) == d,(A,B,d)

def test_ungapkmers_PMC4362205_figS1B_ACnGT():
    import pyhuddinge as ph
    A="ACGT"
    ex = {"ACCT":1,
        "CGTT":1}
    for B,d in ex.items():
        print(A,B)
        assert ph.huddinge_distance(A,B) == d,(A,B,d)



def test_all_huddinge_pairs1():
    import pyhuddinge as ph

 
    kmers = ["AAGGGAA","AAAGCAA","AAACCAA","AAACAAA"]
    import itertools as it

    kmers = ["".join(x) for x in it.product("ACGT",repeat=6)]

    D = ph.all_pairs_huddinge_distance(kmers)

    assert len(D) == len(kmers)*(len(kmers)-1)/2
    assert min(D)>0

@settings(deadline=None)
@given(lists(text(ALPHABET,min_size=1,max_size=15).filter(lambda kmer:len(kmer.strip("n"))>0),min_size=1,max_size=20,unique=True),booleans())
def test_all_huddinge_pairs_concordance(kmers,reverse_complements):
    import pyhuddinge as ph
    print(repr(kmers),file=open("dump.py","w"))
    
    D = ph.all_pairs_huddinge_distance(kmers,reverse_complements)

    assert len(D) == len(kmers)*(len(kmers)-1)/2

    Dsq = ph.squareform(D,kmers)

    for i in Dsq.index:
        for j in Dsq.columns:
            if i<j:
                assert Dsq.loc[i,j] == ph.huddinge_distance(i,j,reverse_complements)

    
 
@settings(deadline=None)
@given(lists(text(ALPHABET,min_size=1,max_size=15).filter(lambda kmer:len(kmer.strip("n"))>0),min_size=1,
    max_size=20,unique=True))
def test_all_huddinge_pairs_shuffle(kmers):
    import pyhuddinge as ph

 
    D1 = ph.all_pairs_huddinge_distance(kmers)
    D2 = ph.all_pairs_huddinge_distance(kmers[::-1])


    Dsq1 = ph.squareform(D1,kmers)
    Dsq2 = ph.squareform(D2,kmers[::-1])
    reordered2 = Dsq2[Dsq1.columns].loc[Dsq1.index]
    
    #print(Dsq1)
    #print(reordered2)
    assert (Dsq1==reordered2).all().all()
    
 
@given(text(ALPHABET,min_size=1,max_size=15))
def test_huddingeneighbor(x):
    import pyhuddinge as ph
    x=x.upper()
    
    h_neigh = ph.huddinge_neighbours(x)
    assert x in h_neigh
    assert len(h_neigh) > 2

#Slow test so only few runs.
@given(text(ALPHABET,min_size=1,max_size=5),integers(min_value=1,max_value=4))
def test_huddingeneighbors_more_by_distance(x,d):

    import pyhuddinge as ph
    x = x.upper()
    
    h_neigh = ph.huddinge_neighbours(x,max_dist=d)
    h_neigh_bigger = ph.huddinge_neighbours(x,max_dist=d+1)
    assert set(h_neigh).issubset(set(h_neigh_bigger))

    for k in h_neigh_bigger:
        if k in h_neigh:
                assert h_neigh[k] == h_neigh_bigger[k]
                assert h_neigh[k] <=d
        else:
                assert h_neigh_bigger[k] == d+1


