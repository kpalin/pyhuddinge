
#include "MODER/huddinge.hpp"
#include "MODER/common.hpp"

#include <python3.6m/Python.h>

extern "C" int pyhuddinge_distance(const char *x, const char *y,const bool reverse_complements)
{  // Compute huddinge distance between x and y. If reverse_complements, compute as minimum of the distances
    std::string cppx(x);
    std::string cppy(y);
    int d= -1;

    d=huddinge_distance(cppx,cppy);
    if(reverse_complements) {
        int rc_d = huddinge_distance(cppx,reverse_complement(cppy));
        d = std::min(d,rc_d);
    }

    return d;
}


extern "C" int pyhuddinge_all_pairs_distance(void *distances, 
        const int m,const char **kmers,const bool reverse_complements)
{
    /* Compute huddinge distances for all pairs of m kmers. 
    If reverse_complements, the distance is the minimum between 
    the kmers and their reverse compelemtns. The 
    distances array must be preallocated to hold m*(m-1)/2 ints.

    Returns a condensed distance matrix Y. 
    For each i and j (where i<j<m). The metric 
    dist(u=X[i], v=X[j]) is computed and stored in entry i*j.

    This ordering is the same as produced by 
    scipy.spatial.distance.pdist and can be converted to
    distance matrix with squareform
*/

    char *D = (char*)distances;

    #pragma omp parallel for if(m>50)
    for(int i=0;i<m;i++) {
        int mmi = m-i;
        int ij=m*(m-1)/2  - mmi*(mmi-1)/2;
        //int base_ij = ij;
        //std::cout << ij<< "-" << ij_alt <<"==" 
        //std::cout<<ij-ij_alt<<'\n';
        const std::string kmer_i(kmers[i]);
        const std::string rc_kmer_i(reverse_complement(kmer_i));

        for(int j=i+1;j<m;j++,ij++){
            const std::string kmer_j(kmers[j]);

            //int ij = (base_ij +j - i -1);
            D[ij] = (char)huddinge_distance(kmer_i,kmer_j);

            if(reverse_complements) {
                D[ij] = std::min(D[ij],(char)huddinge_distance(rc_kmer_i,kmer_j));
            }

            //assert( (alt_ij +j - i -1)==ij);
            //assert(0);

            //D[ij] = (char)ij;
            //std::cout<<i<<":"<<kmer_i<<"-"<<j<<":"<<kmer_j<<" == "<<(int)D[ij]<<'\n';            
        }
    }
    return 1;
}

extern "C" PyObject *pyhuddinge_neighbourhood(const char *kmer,
    int max_dist,int min_kmer_len,int max_kmer_len, int max_gapped_len)
{
    /* Compute sequence neighbourhood of kmer up to huddinge distance max_dist.
    The considered kmers are form min_kmer_len to max_kmer_len long. 
    max_gapped_len is something similar to max_kmer_len but I'm not sure.
*/
    const std::string kmer_s(kmer);
    PyObject* data = PyDict_New();

    int h=max_dist, L=max_gapped_len;

    huddinge_neighbourhood neighbour_obj(kmer_s, h,L,min_kmer_len,max_kmer_len);
    std::vector<boost::tuple<std::string,int> > neighbours = neighbour_obj.compute(true);

    for(std::vector<boost::tuple<std::string,int> >::iterator it = neighbours.begin(); it != neighbours.end(); ++it) {
        const char *n_kmer = it->get<0>().c_str();
        const long dist = it->get<1>();
        
        PyMapping_SetItemString(data,n_kmer,
                PyLong_FromLong(dist));
    }
    
    return data;
}

