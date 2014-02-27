

#ifndef VIENNA_H
#define VIENNA_H

#include <string>
#include "rnainterval.h"
#include <vector>

//include external lib
extern "C" 
{
    #include "vienna/fold.h"
    #include "vienna/Lfold.h"
    #include "vienna/part_func.h"
}


double * boltzmann_fold (const std::string& rna)
{
    double fe = pf_fold (rna.c_str(), NULL);
    double * bppm = export_bppm();
    for (int i = 0; i < rna.size(); ++i)
        for (int j = 0; j < rna.size(); ++j)
            std::cout << bppm[i + j] << " " << std::endl;
    return bppm;
}


RNAInterval zuker_fold(const std::string& rna, const int l, const int r)
{
    //get the interval string and calculate the score
	int sz = r-l+1;
    char* ss = new char[sz+1];
    double score = fold(rna.substr(l, sz).c_str(), ss);
    RNAInterval ret(l, r, -score, ss);
    delete[] ss;
    return ret;
}

std::vector<RNAInterval> zuker_multi_fold(const std::string& rna, const int windowSize) {
	std::vector<RNAInterval> ints;
	for(int i = 0; i < rna.size()-windowSize; ++i) {
		ints.push_back(zuker_fold(rna, i, i+windowSize));
	}
	return ints;
}

std::vector<RNAInterval> rnal_fold(const std::string& rna, int windowSize)
{
    char* ss = new char[rna.size()+1];
    struct mnode* n = mLfold(rna.c_str(), ss, windowSize);
    delete[] ss;
    std::vector<RNAInterval> tmp;
    std::string s;
    while(n->next != 0) 
    {
    	s = std::string(n->sstruct);

        //global left/right index
    	int ll = n->left - 1;
        int rr = ll + s.size() - 1;

        //trim dangling ends
    	int l = ll, r = rr; //l and r are used to store substring indexes in 's'
    	int i;
    	for(i = 0; s[i] == '.'; ++i)
    		++ll;
    	l = i;
    	for(i = s.size()-1; s[i] == '.'; --i)
    		--rr;
    	r = i;

        tmp.push_back ( RNAInterval ( ll, rr, -(n->fe), s.substr (l, (r - l) + 1) ) );
        struct mnode* prev = n;
        n = n->next;
        free (prev->sstruct);
        free (prev);
    }
    return tmp;
}



#endif