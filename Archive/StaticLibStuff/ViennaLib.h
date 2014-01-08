#pragma once
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
extern "C"
{
    #include "Lfold.h"
    #include "fold.h"
}

/*
max@max-ThinkPad-Edge-E330:~/Desktop/Honours/StaticLibStuff$ gcc -c -o Lfold.o Lfold.c
max@max-ThinkPad-Edge-E330:~/Desktop/Honours/StaticLibStuff$ ar  rcs libLfold.a      Lfold.o gquad.o params.o fold.o utils.o energy_par.o ribo.o fold_vars.o aln_util.o alifold.o svm.o svm_utils.o
*/

using namespace std;

struct SStruct
{
    float fe;
    string dbNotation;
    int left, right;
    SStruct(const float f, const string& dbn)
    {
        this->fe = f;
        this->dbNotation = dbn;
        this->left = 0;
        this->right = dbn.size()-1;
    }
    SStruct(const float f, const string& dbn, int l, int r)
    {
        this->fe = f;
        this->dbNotation = dbn;
        this->left = 0;
        this->right = dbn.size()-1;
    }
};

SStruct RNAfold(const string& rna);

vector<SStruct> RNALfold(const string& rna, int windowSize);

float calcFreeEnergy(const string& rna, const string& sstruct);

float calcMoveEnergy(const string& rna, const string& sstruct, int i, int j);