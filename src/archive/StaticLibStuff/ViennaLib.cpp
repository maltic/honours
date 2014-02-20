#include "ViennaLib.h"

SStruct RNAfold(const string& rna)
{
    char* ss = new char[rna.size()+1];
    float score = fold(rna.c_str(), ss);
    return SStruct(score, string(ss)); 
}

vector<SStruct> RNALfold(const string& rna, int windowSize)
{
    char* ss = new char[rna.size()+1];
    struct mnode* n = mLfold(rna.c_str(), ss, windowSize);
    vector<SStruct> tmp;
    while(n->next != 0) 
    {
        tmp.push_back(SStruct(n->fe, string(n->sstruct)));
        n = n->next;
    }
    return tmp;
}

float calcFreeEnergy(const string& rna, const string& sstruct)
{
    return energy_of_structure(rna.c_str(), sstruct.c_str(), 0);
}

float calcMoveEnergy(const string& rna, const string& sstruct, int i, int j)
{
    return energy_of_move(rna.c_str(), sstruct.c_str(), i, j);
}


