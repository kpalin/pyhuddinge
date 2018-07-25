
#include "MODER/huddinge.hpp"

extern "C" int pyhuddinge_distance(const char *x, const char *y)
{
    std::string cppx(x);
    std::string cppy(y);
    int d= -1;

    d=huddinge_distance(cppx,cppy);

    return d;
}


extern "C" int pyhuddinge_all_pairs_distance(void *distances, 
        const int m,const char **kmers)
{
    /* Compute huddinge distances for all pairs of m kmers. The 
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
        
        for(int j=i+1;j<m;j++,ij++){
            const std::string kmer_j(kmers[j]);

            //int ij = (base_ij +j - i -1);
            D[ij] = (char)huddinge_distance(kmer_i,kmer_j);
            //assert( (alt_ij +j - i -1)==ij);
            //assert(0);

            //D[ij] = (char)ij;
            //std::cout<<i<<":"<<kmer_i<<"-"<<j<<":"<<kmer_j<<" == "<<(int)D[ij]<<'\n';            
        }
    }
    return 1;
}

