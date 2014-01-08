

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
}



RNAInterval zuker_fold(const std::string& rna, const int l, const int r)
{
    //get the interval string and calculate the score
	int sz = r-l+1;
    char* ss = new char[sz+1];
    double score = fold(rna.substr(l, sz).c_str(), ss);
    return RNAInterval(l, r, -score, ss);
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
    std::vector<RNAInterval> tmp;
    std::string s;
    while(n->next != 0) 
    {
    	s = std::string(n->sstruct);
    	int ll = n->left-1, rr = n->left+s.size()-1;
    	int l = ll, r = rr;
    	int i;
    	for(i = 0; s[i] == '.'; ++i)
    		++ll;
    	l = i;
    	for(i = s.size()-1; s[i] == '.'; --i)
    		--rr;
    	r = i;
        tmp.push_back(RNAInterval(ll, rr, -(n->fe), s.substr(l, r)));
        n = n->next;
    }
    return tmp;
}



#endif